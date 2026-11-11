use crate::error::Result;

use crate::dp::MAX_TRANSITION_SCORE;
use crate::types::{Base, DsmId, Energy, BASE_COUNT};

use super::{flat_idx, DsmRegistry, DsmTable, DSM_FLAT_SIZE};

const BASES: [Base; BASE_COUNT] = [Base::Gap, Base::A, Base::C, Base::G, Base::N, Base::U];

/// Penalty-adjusted flat scoring table for one orientation.
///
/// External inputs and outputs are kcal/mol.  Internally uses integer score
/// units: `score = -delta_g * 10000`.
#[derive(Clone, Debug)]
pub struct ScoringModel {
    table: [i32; DSM_FLAT_SIZE],
    initiation: Energy,
}

impl ScoringModel {
    #[inline]
    fn ext_penalty_mult(q1: Base, q2: Base, t1: Base, t2: Base) -> i32 {
        if q2 == Base::Gap && t2 == Base::Gap {
            0
        } else if q2 == Base::Gap || t2 == Base::Gap {
            1
        } else if q1 == Base::Gap && t1 == Base::Gap {
            if q2.pair_type(t2).is_match(true) {
                2
            } else {
                0
            }
        } else {
            2
        }
    }

    /// Load the bundled table for `id` at `temperature` and apply `penalty`.
    pub(crate) fn load(id: &DsmId, temperature: i32, penalty: Energy) -> Result<Self> {
        let (initiation, source_table) = DsmRegistry::load(id, temperature)?;
        Ok(Self::new(&source_table, initiation, penalty))
    }

    /// Create a new scoring model from a source table and initiation energy.
    pub fn new(source_table: &DsmTable, initiation: Energy, penalty: Energy) -> Self {
        let mut table = [0i32; DSM_FLAT_SIZE];
        for q1 in BASES {
            for q2 in BASES {
                for t1 in BASES {
                    for t2 in BASES {
                        let idx = flat_idx(q1.as_u8(), q2.as_u8(), t1.as_u8(), t2.as_u8());
                        let ext_penalty_mult = Self::ext_penalty_mult(q1, q2, t1, t2);

                        table[idx] = source_table[q1.as_usize()][q2.as_usize()][t1.as_usize()]
                            [t2.as_usize()]
                            - penalty.0 * ext_penalty_mult;
                    }
                }
            }
        }

        for (i, &val) in table.iter().enumerate() {
            assert!(
                (val as i64).unsigned_abs() <= MAX_TRANSITION_SCORE as u64,
                "table entry {} = {} exceeds MAX_TRANSITION_SCORE bound {}",
                i,
                val,
                MAX_TRANSITION_SCORE,
            );
        }

        Self { table, initiation }
    }

    /// Convert a total stacking score to binding free energy.
    pub(crate) fn binding_energy(&self, stacking_score: Energy) -> Energy {
        self.initiation - stacking_score
    }

    #[inline(always)]
    pub(crate) fn table_ptr(&self) -> *const i32 {
        self.table.as_ptr()
    }

    /// Produce a left-canonical (transposed) copy.
    pub fn transpose(&self) -> ScoringModel {
        let mut table = [0i32; DSM_FLAT_SIZE];
        for q1 in 0..BASE_COUNT as u8 {
            for q2 in 0..BASE_COUNT as u8 {
                for t1 in 0..BASE_COUNT as u8 {
                    for t2 in 0..BASE_COUNT as u8 {
                        let dst = flat_idx(q1, q2, t1, t2);
                        let src = flat_idx(q2, q1, t2, t1);
                        table[dst] = self.table[src];
                    }
                }
            }
        }
        Self {
            table,
            initiation: self.initiation,
        }
    }

    /// Point query for one stacking score in the canonical 4-base tensor.
    #[inline(always)]
    pub(crate) fn score(&self, q1: u8, q2: u8, t1: u8, t2: u8) -> i32 {
        debug_assert!(q1 < 6 && q2 < 6 && t1 < 6 && t2 < 6);
        let idx = flat_idx(q1, q2, t1, t2);
        // SAFETY: All args are in 0..6. Max idx = 1295 < 1296.
        unsafe { *self.table.get_unchecked(idx) }
    }

    /// Point query for one stacking score using semantic bases.
    #[inline(always)]
    pub(crate) fn score_bases(&self, q1: Base, q2: Base, t1: Base, t2: Base) -> i32 {
        self.score(q1.as_u8(), q2.as_u8(), t1.as_u8(), t2.as_u8())
    }

    /// Calculate the thermodynamic score of a continuous, ungapped anti-parallel duplex.
    ///
    /// Query and target are both in duplex-column order: query 5'→3' and the
    /// physical target 3'→5', so paired positions advance together.
    pub(crate) fn ungapped_duplex_score(&self, query: &[Base], target: &[Base]) -> Energy {
        assert_eq!(query.len(), target.len());
        Energy(
            query
                .windows(2)
                .zip(target.windows(2))
                .map(|(q, t)| self.score_bases(q[0], q[1], t[0], t[1]))
                .sum(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_source() -> DsmTable {
        [[[[0; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]
    }

    #[test]
    fn runtime_table_and_ungapped_score_use_physical_target_order() {
        let mut source = empty_source();
        source[Base::A.as_usize()][Base::C.as_usize()][Base::U.as_usize()][Base::G.as_usize()] = 11;
        source[Base::C.as_usize()][Base::G.as_usize()][Base::G.as_usize()][Base::C.as_usize()] = 17;
        let model = ScoringModel::new(&source, Energy::from_kcal(0.0), Energy::from_kcal(0.0));

        assert_eq!(model.score_bases(Base::A, Base::C, Base::U, Base::G), 11);
        assert_eq!(model.score_bases(Base::A, Base::C, Base::A, Base::C), 0);
        assert_eq!(
            model.ungapped_duplex_score(&[Base::A, Base::C, Base::G], &[Base::U, Base::G, Base::C]),
            Energy(28)
        );
    }

    #[test]
    fn anchor_penalty_uses_physical_wc_wobble_and_mismatch_pairs() {
        let source = empty_source();
        let unit = Energy::from_kcal(0.005);
        let charged = ScoringModel::new(&source, Energy::from_kcal(0.0), unit);

        for query in BASES {
            for target in BASES {
                let expected_mult = match (query == Base::Gap, target == Base::Gap) {
                    (true, true) => 0,
                    (true, false) | (false, true) => 1,
                    (false, false) => 2 * i32::from(query.pair_type(target).is_match(true)),
                };
                assert_eq!(
                    charged.score_bases(Base::Gap, query, Base::Gap, target),
                    -unit.0 * expected_mult,
                    "anchor penalty for physical pair {query:?}-{target:?}"
                );
            }
        }
    }
}
