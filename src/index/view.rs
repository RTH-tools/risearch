//! Borrowed target-index views and target-coordinate mapping.

use std::ops::Range;

use crate::index::sa::SuffixIndexView;
use crate::types::{Base, Strand};

/// Zero-copy target index used by seed search.
///
/// The paired suffix index covers the concatenation
/// `R(T) + Gap + C(T) + Gap` for every target `T`. The offset directory maps
/// positions in that global sequence back to a target, physical strand, and
/// strand-local coordinate. A block ends at the next target offset, or at the
/// end of the indexed sequence for the final target.
///
/// All fields are borrowed from the same validated [`crate::TargetRegistry`].
#[derive(Clone, Copy)]
pub struct TargetView<'a> {
    suffixes: SuffixIndexView<'a>,
    /// Start offset of each target block in the indexed sequence.
    offsets: &'a [usize],
}

impl<'a> TargetView<'a> {
    /// Construct a target view from a suffix index and its target directory.
    ///
    /// `offsets[i]` must point to the start of target `i`'s
    /// `R(T) + Gap + C(T) + Gap` block. The offsets must be strictly increasing
    /// and collectively cover the indexed sequence.
    pub(crate) fn new(suffixes: SuffixIndexView<'a>, offsets: &'a [usize]) -> Self {
        Self { suffixes, offsets }
    }

    /// Return the paired sequence and suffix positions used for traversal.
    #[inline]
    pub(crate) fn suffixes(self) -> SuffixIndexView<'a> {
        self.suffixes
    }

    /// Return the selected physical target strand in duplex-column order.
    ///
    /// For input `T` written 5' to 3', Forward is `R(T)` and Reverse is `C(T)`.
    ///
    /// # Panics
    ///
    /// Panics if `target_idx` is outside the target directory.
    #[inline]
    pub fn target(&self, target_idx: usize, strand: Strand) -> &'a [Base] {
        let seq_len = self.seq_len(target_idx);
        let offset = self.offsets[target_idx]
            + match strand {
                Strand::Forward => 0,
                Strand::Reverse => seq_len + 1,
            };
        &self.suffixes.sequence()[offset..offset + seq_len]
    }

    /// Convert a half-open target range between physical duplex-view and
    /// original FASTA coordinates.
    ///
    /// This is an involution: Forward `R(T)` mirrors; Reverse `C(T)` preserves.
    #[inline(always)]
    pub(crate) fn map_target_range(
        self,
        target_idx: usize,
        strand: Strand,
        range: Range<usize>,
    ) -> Range<usize> {
        let target_len = self.seq_len(target_idx);
        debug_assert!(range.start <= range.end && range.end <= target_len);
        match strand {
            Strand::Forward => target_len - range.end..target_len - range.start,
            Strand::Reverse => range,
        }
    }

    /// Map a global sequence position to a physical target coordinate.
    ///
    /// The returned start is relative to the physical strand selected by
    /// [`Strand`]. Returns `None` when `global_pos` is outside a physical strand,
    /// points at a separator [`Base::Gap`], or leaves fewer than `seed_len`
    /// bases before the next separator.
    #[inline]
    pub(crate) fn map_seed_pos(
        self,
        global_pos: usize,
        seed_len: usize,
    ) -> Option<(usize, Strand, usize)> {
        let target_idx = self
            .offsets
            .partition_point(|&offset| offset <= global_pos)
            .checked_sub(1)?;
        let local_pos = global_pos - self.offsets[target_idx];
        let seq_len = self.seq_len(target_idx);

        let (strand, start) = if local_pos < seq_len {
            (Strand::Forward, local_pos)
        } else {
            (Strand::Reverse, local_pos.checked_sub(seq_len + 1)?)
        };
        (start < seq_len && seed_len <= seq_len - start).then_some((target_idx, strand, start))
    }

    /// Derive one physical strand's length from its enclosing target block.
    #[inline]
    fn seq_len(self, target_idx: usize) -> usize {
        let block_start = self.offsets[target_idx];
        let block_end = self
            .offsets
            .get(target_idx + 1)
            .copied()
            .unwrap_or(self.suffixes.sequence().len());
        (block_end - block_start - 2) / 2
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Range;

    use crate::index::store::TargetRegistry;
    use crate::types::Strand;
    use crate::Sequence;

    /// Independent oracle for duplex-frame → FASTA coordinate conversion.
    ///
    /// Derived from [`crate::seed::reference_tests::target_in_duplex_order`]:
    /// Forward = `reverse(FASTA)`, so duplex pos `d` = FASTA pos `len-1-d`.
    /// Reverse = `complement(FASTA)`, same positions.
    fn naive_duplex_to_fasta(
        target_len: usize,
        strand: Strand,
        duplex_range: Range<usize>,
    ) -> Range<usize> {
        match strand {
            Strand::Forward => (target_len - duplex_range.end)..(target_len - duplex_range.start),
            Strand::Reverse => duplex_range,
        }
    }

    #[test]
    fn map_target_range_matches_independent_coordinate_oracle() {
        let target_text = "ACGUACGUACGU";
        let target_len = target_text.len();
        let (target, _) = Sequence::normalize("t", target_text.as_bytes()).unwrap();
        let targets = TargetRegistry::build(vec![("t".into(), target)], None).unwrap();
        let tview = targets.view();

        for &strand in &[Strand::Forward, Strand::Reverse] {
            let ranges: &[Range<usize>] = &[
                0..1,
                0..target_len,
                0..target_len / 2,
                target_len / 2..target_len,
                3..7,
                target_len - 1..target_len,
                5..5,
            ];
            for range in ranges {
                let production = tview.map_target_range(0, strand, range.clone());
                let oracle = naive_duplex_to_fasta(target_len, strand, range.clone());
                assert_eq!(
                    production, oracle,
                    "strand={strand} duplex_range={range:?} target_len={target_len}"
                );
            }
        }
    }
}

#[cfg(kani)]
mod kani_proofs {
    use super::*;

    // One valid four-base target block: R(T) + Gap + C(T) + Gap for
    // T = AAGC. Distinct strand contents expose an offset/orientation mixup.
    // The suffix array is empty because these proofs exercise only
    // target geometry; no traversal method reads it.
    static SEQUENCE: [u8; 10] = [2, 3, 1, 1, 0, 5, 5, 2, 3, 0];
    static SUFFIX_ARRAY: [u64; 0] = [];
    static OFFSETS: [usize; 1] = [0];

    fn fixture() -> TargetView<'static> {
        // SAFETY: every byte is a valid Base discriminant, and the static
        // suffix array and sequence share the same immutable lifetime.
        let suffixes = unsafe { SuffixIndexView::from_bytes_unchecked(&SEQUENCE, &SUFFIX_ARRAY) };
        TargetView::new(suffixes, &OFFSETS)
    }

    /// Check physical strand contents and seed containment against a literal
    /// archive layout, including positions on/beyond its separators.
    #[kani::proof]
    fn target_view_preserves_bounded_ranges_and_span_lengths() {
        let view = fixture();
        assert_eq!(
            view.target(0, Strand::Forward),
            &[Base::C, Base::G, Base::A, Base::A]
        );
        assert_eq!(
            view.target(0, Strand::Reverse),
            &[Base::U, Base::U, Base::C, Base::G]
        );

        let start: usize = kani::any();
        let end: usize = kani::any();
        kani::assume(start <= end);
        kani::assume(end <= 4);
        let range = start..end;

        let forward = view.map_target_range(0, Strand::Forward, range.clone());
        assert_eq!(forward, 4 - end..4 - start);

        let position: usize = kani::any();
        let seed_len: usize = kani::any();
        kani::assume(position <= 10);
        kani::assume(seed_len >= 1 && seed_len <= 5);
        let expected = if position < 4 && seed_len <= 4 - position {
            Some((0, Strand::Forward, position))
        } else if (5..9).contains(&position) && seed_len <= 9 - position {
            Some((0, Strand::Reverse, position - 5))
        } else {
            None
        };
        assert_eq!(view.map_seed_pos(position, seed_len), expected);
    }
}
