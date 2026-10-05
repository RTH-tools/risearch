//! Parse a user-provided long-form DSM TSV into an initiation offset and table.

use super::{DsmTable, DSM_HEADER};
use crate::error::{Error, Result};
use crate::types::{Base, Energy, BASE_COUNT};

const UNOBSERVED_KCAL: f64 = 20.0;

/// Parse a `q1 q2 t1 t2 delta_g_kcal_per_mol` TSV into `(initiation, table)`.
///
/// Absent cells and cells with ΔG >= 20 kcal/mol are unobserved. Input ΔG is
/// taken at offset 0; the initiation offset is derived from the initiation
/// cells and folded into them; a table without initiation cells is rejected.
pub(crate) fn parse_tsv(input: &str) -> Result<(Energy, DsmTable)> {
    let input = input.strip_prefix('\u{feff}').unwrap_or(input);
    let mut lines = input.lines().enumerate();
    let (_, header) = lines
        .next()
        .ok_or_else(|| Error::Dsm("empty energy parameter file".into()))?;
    if header.split('\t').ne(DSM_HEADER) {
        return Err(Error::Dsm(format!(
            "bad energy parameter file header '{header}', expected '{}'",
            DSM_HEADER.join("\t")
        )));
    }

    let mut cells = [[[[None::<f64>; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
    for (n, line) in lines {
        let row = n + 1;
        let fields: Vec<&str> = line.split('\t').collect();
        let &[q1, q2, t1, t2, dg] = fields.as_slice() else {
            return Err(Error::Dsm(format!(
                "line {row}: expected 5 tab-separated fields, found {}",
                fields.len()
            )));
        };
        let [q1, q2, t1, t2] = [q1, q2, t1, t2].map(|f| base(f, row));
        let (q1, q2, t1, t2) = (q1?, q2?, t1?, t2?);
        let dg: f64 = dg
            .parse()
            .map_err(|e| Error::Dsm(format!("line {row}: invalid energy '{dg}': {e}")))?;
        let cell = &mut cells[q1.as_usize()][q2.as_usize()][t1.as_usize()][t2.as_usize()];
        if cell.is_some() {
            return Err(Error::Dsm(format!(
                "line {row}: duplicate cell {q1:?}{q2:?}/{t1:?}{t2:?}"
            )));
        }
        *cell = Some(dg);
    }

    let mut peak = f64::NEG_INFINITY;
    for_each_cell(|q1, q2, t1, t2| {
        if let Some(v) = cells[q1][q2][t1][t2] {
            if is_init(q1, q2, t1, t2) && v < UNOBSERVED_KCAL {
                peak = peak.max(v);
            }
        }
    });
    if !peak.is_finite() {
        return Err(Error::Dsm(
            "energy parameter file has no initiation rows (X - X - or - X - X)".into(),
        ));
    }
    let offset = peak * 2.0 + 1.0;

    let mut table = [[[[0i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
    let mut err = None;
    for_each_cell(|q1, q2, t1, t2| {
        let val = match cells[q1][q2][t1][t2] {
            Some(v) if v < UNOBSERVED_KCAL && is_init(q1, q2, t1, t2) => v - offset / 2.0,
            Some(v) if v < UNOBSERVED_KCAL => v,
            _ => UNOBSERVED_KCAL,
        };
        match Energy::try_from(-val) {
            Ok(e) => table[q1][q2][t1][t2] = e.0,
            Err(e) => {
                err.get_or_insert(e);
            }
        }
    });
    if let Some(e) = err {
        return Err(Error::Dsm(e));
    }
    Ok((Energy::try_from(offset).map_err(Error::Dsm)?, table))
}

fn base(field: &str, row: usize) -> Result<Base> {
    let mut chars = field.chars();
    match (chars.next(), chars.next()) {
        (Some(c), None) => Base::try_from(c).map_err(|e| Error::Dsm(format!("line {row}: {e}"))),
        _ => Err(Error::Dsm(format!("line {row}: invalid base '{field}'"))),
    }
}

/// Helix initiation: both dinucleotides are `X-` or both are `-X`.
fn is_init(q1: usize, q2: usize, t1: usize, t2: usize) -> bool {
    let gap = Base::Gap.as_usize();
    let lead = q1 != gap && q2 == gap && t1 != gap && t2 == gap;
    let trail = q1 == gap && q2 != gap && t1 == gap && t2 != gap;
    lead || trail
}

fn for_each_cell(mut f: impl FnMut(usize, usize, usize, usize)) {
    for q1 in 0..BASE_COUNT {
        for q2 in 0..BASE_COUNT {
            for t1 in 0..BASE_COUNT {
                for t2 in 0..BASE_COUNT {
                    f(q1, q2, t1, t2);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "q1\tq2\tt1\tt2\tdelta_g_kcal_per_mol\n";

    #[test]
    fn parses_cells_and_derives_initiation_offset() {
        let tsv = format!(
            "{HEADER}A\tC\tU\tG\t-2.1805\nA\t-\tU\t-\t3.0\n-\tA\t-\tU\t2.0\nG\tG\tC\tC\t25.0\n"
        );
        let (init, table) = parse_tsv(&tsv).unwrap();
        assert_eq!(init, Energy(70000));
        assert_eq!(table[1][2][5][3], 21805);
        assert_eq!(table[1][0][5][0], 5000);
        assert_eq!(table[0][1][0][5], 15000);
        assert_eq!(table[3][3][2][2], -200000);
        assert_eq!(table[0][0][0][0], -200000);
        assert_eq!(parse_tsv(&format!("\u{feff}{tsv}")).unwrap(), (init, table));
    }

    #[test]
    fn rejects_malformed_input() {
        let bad = [
            "",
            "q1\tq2\tt1\tt2\tdG\nA\tC\tU\tG\t-1.0\n",
            &format!("{HEADER}A\tC\tU\tG\t-1.0\n"),
            &format!("{HEADER}A\tC\tU\t-1.0\n"),
            &format!("{HEADER}A\tC\tU\tX\t-1.0\n"),
            &format!("{HEADER}AC\tC\tU\tG\t-1.0\n"),
            &format!("{HEADER}A\tC\tU\tG\tabc\n"),
            &format!("{HEADER}A\tC\tU\tG\t-1.0\nA\tC\tU\tG\t-2.0\n"),
        ];
        for tsv in bad {
            assert!(matches!(parse_tsv(tsv), Err(Error::Dsm(_))), "{tsv:?}");
        }
    }
}
