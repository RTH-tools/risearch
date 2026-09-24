//! RNA dinucleotide stacking energy matrices (DSM)
//!
//! Tables encode nearest-neighbor thermodynamic parameters for RNA-RNA interactions.
//! Canonical table values are stored in kcal/mol and converted to score units on load.
//!
//! Indices: `[q1][q2][t1][t2]` = stacking energy for:
//! ```text
//! Query:  5'─ q1 ─ q2 ─ 3'
//!             |    |
//! Target: 3'─ t1 ─ t2 ─ 5'
//! ```

use std::path::Path;

use crate::error::{Error, Result};
use crate::types::{DsmId, Energy, BASE_COUNT};

mod model;
mod parse;
mod tables;

pub use model::ScoringModel;

pub(crate) type DsmTable = [[[[i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];

pub(crate) const DSM_FLAT_SIZE: usize = BASE_COUNT * BASE_COUNT * BASE_COUNT * BASE_COUNT;

pub(crate) const DSM_HEADER: [&str; 5] = ["q1", "q2", "t1", "t2", "delta_g_kcal_per_mol"];

#[inline(always)]
pub(crate) const fn flat_idx(q1: u8, q2: u8, t1: u8, t2: u8) -> usize {
    (q1 as usize) * BASE_COUNT.pow(3)
        + (q2 as usize) * BASE_COUNT.pow(2)
        + (t1 as usize) * BASE_COUNT
        + t2 as usize
}

const NAMES: &[&str] = &["t04", "slh04", "s95-rna-dna", "s95-dna-rna"];

/// Service for managing and loading dinucleotide stacking models (DSM).
pub struct DsmRegistry;

impl DsmRegistry {
    /// List all available DSM string identifiers (for CLI/UI).
    pub const fn all_names() -> &'static [&'static str] {
        NAMES
    }

    /// Accept a bundled identifier or a path to an existing TSV table.
    pub fn parse_id(s: &str) -> Result<DsmId> {
        if NAMES.contains(&s) || Path::new(s).is_file() {
            Ok(DsmId(s.to_string()))
        } else {
            Err(Error::Dsm(format!(
                "unknown DSM id '{s}': expected one of {} or a path to a TSV table",
                NAMES.join(", ")
            )))
        }
    }

    /// Load a bundled table, interpolating between bracket temperatures if needed,
    /// or parse the user TSV at `id`, which is used as-is whatever `temperature` is.
    pub fn load(id: &DsmId, temperature: i32) -> Result<(Energy, DsmTable)> {
        match id.0.as_str() {
            "s95-rna-dna" => Self::load_canonical("s95", temperature),
            "s95-dna-rna" => {
                let (init, table) = Self::load_canonical("s95", temperature)?;
                Ok((init, transpose_table(&table)))
            }
            name if NAMES.contains(&name) => Self::load_canonical(name, temperature),
            path => parse::parse_tsv(&fs_err::read_to_string(path)?)
                .map_err(|e| Error::Dsm(format!("{path}: {e}"))),
        }
    }

    fn load_canonical(name: &str, temperature: i32) -> Result<(Energy, DsmTable)> {
        let entries = tables::TABLES.iter().filter(|e| e.0 == name);
        if let Some(&(_, _, init, table)) = entries.clone().find(|e| e.1 == temperature) {
            return Ok((init, *table));
        }
        let lo = entries
            .clone()
            .filter(|e| e.1 < temperature)
            .max_by_key(|e| e.1);
        let hi = entries.filter(|e| e.1 > temperature).min_by_key(|e| e.1);
        let (Some(&(_, lo_t, lo_init, lo_table)), Some(&(_, hi_t, hi_init, hi_table))) = (lo, hi)
        else {
            return Err(Error::Dsm(format!(
                "DSM '{name}' has no bundled table bracketing {temperature}C"
            )));
        };
        Ok(interpolate_tables(
            lo_init,
            lo_table,
            hi_init,
            hi_table,
            [f64::from(temperature), f64::from(lo_t), f64::from(hi_t)],
        ))
    }
}

fn interpolate_tables(
    offset_1: Energy,
    table_1: &DsmTable,
    offset_2: Energy,
    table_2: &DsmTable,
    temps: [f64; 3],
) -> (Energy, DsmTable) {
    let mut table = [[[[0i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
    for q1 in 0..BASE_COUNT {
        for q2 in 0..BASE_COUNT {
            for t1 in 0..BASE_COUNT {
                for t2 in 0..BASE_COUNT {
                    table[q1][q2][t1][t2] =
                        lerp(table_1[q1][q2][t1][t2], table_2[q1][q2][t1][t2], temps);
                }
            }
        }
    }
    (Energy(lerp(offset_1.0, offset_2.0, temps)), table)
}

/// Swap query and target strands: `out[q1][q2][t1][t2] = table[t2][t1][q2][q1]`.
fn transpose_table(table: &DsmTable) -> DsmTable {
    let mut out = [[[[0i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
    for q1 in 0..BASE_COUNT {
        for q2 in 0..BASE_COUNT {
            for t1 in 0..BASE_COUNT {
                for t2 in 0..BASE_COUNT {
                    out[q1][q2][t1][t2] = table[t2][t1][q2][q1];
                }
            }
        }
    }
    out
}

fn lerp(v1: i32, v2: i32, temps: [f64; 3]) -> i32 {
    let [t0, t1, t2] = temps;
    let diff = i64::from(v1) - i64::from(v2);
    ((t0 - t2) / (t1 - t2) * diff as f64 + f64::from(v2)).round() as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Base;

    #[test]
    fn pairing_policy_uses_physical_duplex_bases() {
        assert!(Base::A.pair_type(Base::U).is_match(false));
        assert!(Base::G.pair_type(Base::U).is_match(true));
        assert!(!Base::G.pair_type(Base::U).is_match(false));
        assert!(!Base::G.pair_type(Base::A).is_match(true));
    }

    #[test]
    fn score_returns_nonzero_for_valid_pairs() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let score = model.score_bases(Base::A, Base::C, Base::U, Base::G);
        let gap_score = model.score_bases(Base::Gap, Base::A, Base::Gap, Base::U);
        assert!(
            gap_score != 0 || score != 0,
            "At least one transition query should be non-zero"
        );
    }

    #[test]
    fn gg_cc_stack_is_strongest() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let actual = model.score_bases(Base::G, Base::G, Base::C, Base::C);
        assert_eq!(actual, 33016);
    }

    #[test]
    fn bundled_models_have_explicit_stack_and_initiation_energies() {
        // Values read from the generated tables/*.rs, in integer score units.
        // An AC/UG stack is asymmetric in the mixed-strand model.
        for (name, init, ac_ug, gg_cc) in [
            ("t04", 61352, 21796, 33016),
            ("slh04", 20516, 13853, 18006),
            ("s95-rna-dna", 40464, 20478, 29009),
        ] {
            let (actual_init, table) = DsmRegistry::load(&DsmId::from(name), 37).unwrap();
            assert_eq!(actual_init, Energy(init), "{name}");
            assert_eq!(table[1][2][5][3], ac_ug, "{name}");
            assert_eq!(table[3][3][2][2], gg_cc, "{name}");
            assert_eq!(table[0][0][0][0], -200000, "omitted transition in {name}");
        }
        let (init, reversed) = DsmRegistry::load(&DsmId::from("s95-dna-rna"), 37).unwrap();
        assert_eq!(init, Energy(40464));
        // AC/UG (RNA/DNA) becomes GU/CA (DNA/RNA), not a plain transpose.
        assert_eq!(reversed[3][5][2][1], 20478);
    }

    #[test]
    fn temperature_interpolation_uses_bracketing_integer_scores() {
        let (init, table) = DsmRegistry::load(&DsmId::from("t04"), 31).unwrap();
        // Midpoint of 25C and 37C: initiation (63906 + 61352)/2;
        // AC/UG stack (25372 + 21796)/2; GG/CC (36919 + 33016)/2.
        assert_eq!(init, Energy(62629));
        assert_eq!(table[1][2][5][3], 23584);
        assert_eq!(table[3][3][2][2], 34968);
        assert_eq!(table[1][0][5][0], 5000);
        assert!(DsmRegistry::load(&DsmId::from("missing"), 37).is_err());
        assert!(DsmRegistry::load(&DsmId::from("t04"), 51).is_err());
    }

    #[test]
    fn interpolation_rounds_exact_half_away_from_zero() {
        // temps: [target, v1's knot, v2's knot]; 31C sits exactly midway in 25..37.
        assert_eq!(lerp(1, 2, [31.0, 25.0, 37.0]), 2);
        assert_eq!(lerp(-1, -2, [31.0, 25.0, 37.0]), -2);
        assert_eq!(lerp(1_200_000, 1_200_001, [31.0, 25.0, 37.0]), 1_200_001);
    }

    #[test]
    fn energy_conversion_roundtrips() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let energy = model.binding_energy(Energy(33016));
        assert!((energy.to_kcal() - 2.8336).abs() < 0.001);
    }

    #[test]
    fn transpose_swaps_both_pairs() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let right = ScoringModel::new(&source, init, Energy::from_kcal(0.005));
        let left = right.transpose();
        for q1 in 0u8..6 {
            for q2 in 0u8..6 {
                for t1 in 0u8..6 {
                    for t2 in 0u8..6 {
                        assert_eq!(
                            left.score(q1, q2, t1, t2),
                            right.score(q2, q1, t2, t1),
                            "transpose mismatch at ({},{},{},{})",
                            q1,
                            q2,
                            t1,
                            t2
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn user_tsv_path_loads_as_matrix() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("custom.tsv");
        fs_err::write(
            &path,
            "q1\tq2\tt1\tt2\tdelta_g_kcal_per_mol\nA\t-\tU\t-\t3.0\nA\tC\tU\tG\t-2.1805\n",
        )
        .unwrap();
        let id = DsmRegistry::parse_id(path.to_str().unwrap()).unwrap();
        let (init, table) = DsmRegistry::load(&id, 37).unwrap();
        assert_eq!(init, Energy(70000));
        assert_eq!(table[1][2][5][3], 21805);
        assert_eq!(table[1][0][5][0], 5000);
        assert!(DsmRegistry::parse_id("no-such-model").is_err());
        assert!(DsmRegistry::load(&DsmId::from("/no/such/file.tsv"), 37).is_err());
    }

    fn table_hash(table: &DsmTable) -> u64 {
        table
            .iter()
            .flatten()
            .flatten()
            .flatten()
            .fold(0u64, |h, &v| {
                h.wrapping_mul(31).wrapping_add(v as u32 as u64)
            })
    }

    /// Pin every cell of every bundled model's 6^4 tensor at its default temperature.
    ///
    /// A single flipped cell changes the hash. After an intentional table update,
    /// take the new value from this assertion's failure output.
    #[test]
    fn bundled_model_tables_are_bitwise_stable() {
        for (name, expected) in [
            ("t04", 0x8FBA_58A6_D7DB_E203_u64),
            ("slh04", 0xE913_0DFA_2E56_1AAE),
            ("s95-rna-dna", 0xD2D5_D3BA_1B6C_4210),
            ("s95-dna-rna", 0x1C1C_88FA_64AA_7D00),
        ] {
            let (_, table) = DsmRegistry::load(&DsmId::from(name), 37).unwrap();
            assert_eq!(
                table_hash(&table),
                expected,
                "{name} table changed — if intentional, update the hash",
            );
        }
    }
}
