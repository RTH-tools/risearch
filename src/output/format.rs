use std::ops::Range;

use crate::alignment::{fingerprint_symbols, AlignColumn};
use crate::config::OutputFormat;
use crate::search::SearchHit;
use crate::types::Base;

const BINDING_SITE_FLANK_LEN: usize = 20;

struct TsvLine<'a> {
    buf: &'a mut Vec<u8>,
    needs_sep: bool,
}

impl<'a> TsvLine<'a> {
    #[inline(always)]
    fn new(buf: &'a mut Vec<u8>) -> Self {
        Self {
            buf,
            needs_sep: false,
        }
    }

    #[inline(always)]
    fn field(&mut self, bytes: &[u8]) {
        if self.needs_sep {
            self.buf.push(b'\t');
        }
        self.needs_sep = true;
        self.buf.extend_from_slice(bytes);
    }

    #[inline(always)]
    fn field_with(&mut self, f: impl FnOnce(&mut Vec<u8>)) {
        if self.needs_sep {
            self.buf.push(b'\t');
        }
        self.needs_sep = true;
        f(self.buf);
    }

    #[inline(always)]
    fn finish(self) {
        self.buf.push(b'\n');
    }
}

#[inline]
pub fn format_hit_into(
    out: &mut Vec<u8>,
    hit: &SearchHit,
    q_name: &str,
    t_name: &str,
    target: &[Base],
    target_range: Range<usize>,
    format: OutputFormat,
) {
    // Uninitialized 40-byte stack array; there is nothing to amortize by
    // hoisting it to the caller.
    let mut itoa = itoa::Buffer::new();
    let alignment = hit.alignment.as_ref();
    let steps_len = alignment.map(|a| a.len()).unwrap_or(0);
    let flank_reserve = if format == OutputFormat::BindingSite {
        2 * BINDING_SITE_FLANK_LEN
    } else {
        0
    };
    out.reserve(q_name.len() + t_name.len() + (steps_len * 2) + flank_reserve + 96);

    if format == OutputFormat::Detailed {
        write_alignment_prelude(out, hit);
    }

    let mut row = TsvLine::new(out);
    write_base_fields(&mut row, &mut itoa, hit, q_name, t_name);
    write_extended_fields(&mut row, hit, target, target_range, format);
    row.finish();
}

#[inline(always)]
fn write_base_fields(
    row: &mut TsvLine<'_>,
    itoa: &mut itoa::Buffer,
    hit: &SearchHit,
    q_name: &str,
    t_name: &str,
) {
    row.field(q_name.as_bytes());
    row.field(itoa.format(hit.q_start + 1).as_bytes());
    row.field(itoa.format(hit.q_end + 1).as_bytes());
    row.field(t_name.as_bytes());
    row.field(itoa.format(hit.t_start + 1).as_bytes());
    row.field(itoa.format(hit.t_end + 1).as_bytes());
    row.field(&[char::from(hit.strand) as u8]);
    row.field_with(|buf| append_score_2dp(buf, itoa, hit.energy.to_kcal()));
}

#[inline(always)]
fn write_extended_fields(
    row: &mut TsvLine<'_>,
    hit: &SearchHit,
    target: &[Base],
    target_range: Range<usize>,
    format: OutputFormat,
) {
    let alignment = hit.alignment.as_ref();

    if matches!(format, OutputFormat::Cigar | OutputFormat::BindingSite) {
        row.field_with(|buf| {
            if let Some(align) = alignment {
                push_pairing_string(buf, align);
            }
        });
    }

    if format == OutputFormat::BindingSite {
        row.field_with(|buf| {
            if let Some(align) = alignment {
                push_target_row(buf, align);
            }
        });

        // The physical duplex target runs 3'->5': bases after the site are its
        // 5' flank, while the 3' flank runs backward through the preceding bases.
        // Legacy -p3 reports the 5' flank first, then the 3'.
        let before = &target[..target_range.start];
        let after = &target[target_range.end..];
        let flank_5 = &after[..after.len().min(BINDING_SITE_FLANK_LEN)];
        let flank_3 = &before[before.len().saturating_sub(BINDING_SITE_FLANK_LEN)..];
        row.field_with(|buf| buf.extend(flank_5.iter().map(|base| base.to_byte())));
        row.field_with(|buf| buf.extend(flank_3.iter().rev().map(|base| base.to_byte())));
    }
}

fn write_alignment_prelude(out: &mut Vec<u8>, hit: &SearchHit) {
    let Some(align) = hit.alignment.as_ref() else {
        return;
    };
    for column in align {
        out.push(column.query().to_byte());
    }
    out.push(b'\n');

    for column in align {
        out.push(column.class().alignment_symbol() as u8);
    }
    out.push(b'\n');

    push_target_row(out, align);
    out.push(b'\n');
}

#[inline]
fn append_score_2dp(buf: &mut Vec<u8>, itoa: &mut itoa::Buffer, score: f64) {
    let scaled = (score * 100.0).round_ties_even() as i64;
    if scaled < 0 {
        buf.push(b'-');
    }
    let v = scaled.unsigned_abs();
    buf.extend_from_slice(itoa.format(v / 100).as_bytes());
    buf.push(b'.');
    let frac = (v % 100) as u8;
    buf.push(b'0' + (frac / 10));
    buf.push(b'0' + (frac % 10));
}

fn push_pairing_string(buf: &mut Vec<u8>, alignment: &[AlignColumn]) {
    buf.extend(fingerprint_symbols(alignment).map(|symbol| symbol as u8));
}

/// The target row, shared by the detailed block and the binding-site column.
fn push_target_row(buf: &mut Vec<u8>, alignment: &[AlignColumn]) {
    for column in alignment {
        buf.push(column.target().to_byte());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Energy, Strand};

    #[test]
    fn detailed_rows_preserve_asymmetry_gaps_wobble_and_trailing_space() {
        let alignment = [
            AlignColumn::new(Base::A, Base::U),
            AlignColumn::new(Base::C, Base::G),
            AlignColumn::new(Base::G, Base::U),
            AlignColumn::new(Base::U, Base::U),
            AlignColumn::new(Base::G, Base::Gap),
            AlignColumn::new(Base::Gap, Base::A),
            AlignColumn::new(Base::C, Base::G),
            AlignColumn::new(Base::A, Base::A),
        ];
        let hit = SearchHit {
            query_idx: 0,
            target_idx: 0,
            q_start: 2,
            q_end: 8,
            t_start: 10,
            t_end: 16,
            strand: Strand::Reverse,
            energy: Energy::from_kcal(-12.34),
            alignment: Some(Box::new(alignment)),
        };
        let mut out = Vec::new();
        format_hit_into(&mut out, &hit, "q", "t", &[], 0..0, OutputFormat::Detailed);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "acgug-ca\n||:   | \nuguu-aga\nq\t3\t9\tt\t11\t17\t-\t-12.34\n"
        );
    }

    #[test]
    fn binding_site_writes_physical_five_prime_then_three_prime_flank() {
        let hit = SearchHit {
            query_idx: 0,
            target_idx: 0,
            q_start: 0,
            q_end: 1,
            t_start: 2,
            t_end: 3,
            strand: Strand::Reverse,
            energy: Energy::from_kcal(0.0),
            alignment: None,
        };
        let target = [Base::A, Base::C, Base::G, Base::U, Base::A, Base::C];
        let mut out = Vec::new();

        format_hit_into(
            &mut out,
            &hit,
            "q",
            "t",
            &target,
            2..4,
            OutputFormat::BindingSite,
        );

        let text = String::from_utf8(out).unwrap();
        let fields = text.trim_end().split('\t').collect::<Vec<_>>();
        // A zero energy is unsigned; the sign branch must test strictly below zero.
        assert_eq!(fields[7], "0.00", "zero energy must not be signed");
        assert_eq!(fields[10], "ac", "5' flank must be reported first");
        assert_eq!(fields[11], "ca", "3' flank must be reported second");
    }

    #[test]
    fn cigar_rows_carry_the_pairing_symbol_string() {
        let alignment = [
            AlignColumn::new(Base::A, Base::U),
            AlignColumn::new(Base::C, Base::G),
            AlignColumn::new(Base::G, Base::U),
            AlignColumn::new(Base::U, Base::U),
            AlignColumn::new(Base::G, Base::Gap),
            AlignColumn::new(Base::Gap, Base::A),
            AlignColumn::new(Base::C, Base::G),
            AlignColumn::new(Base::A, Base::A),
        ];
        let hit = SearchHit {
            query_idx: 0,
            target_idx: 0,
            q_start: 2,
            q_end: 8,
            t_start: 10,
            t_end: 16,
            strand: Strand::Reverse,
            energy: Energy::from_kcal(-12.34),
            alignment: Some(Box::new(alignment)),
        };
        let mut out = Vec::new();

        format_hit_into(&mut out, &hit, "q", "t", &[], 0..0, OutputFormat::Cigar);

        let text = String::from_utf8(out).unwrap();
        let fields = text.trim_end_matches('\n').split('\t').collect::<Vec<_>>();
        // Class codes, not the detailed block's glyphs: canonical, canonical,
        // wobble, mismatch, query bulge, target bulge, canonical, mismatch.
        assert_eq!(
            fields[8], "PPWUQTPU",
            "pairing column must carry one class code per alignment column"
        );
    }
}
