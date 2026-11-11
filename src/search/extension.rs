use std::ops::Range;

use log::trace;

use crate::adapter::gotoh::GotohModel;
use crate::alignment::{self, AlignColumn};
use crate::dp::gotoh::Gotoh;
use crate::dp::{DpGrid, TraceOp, MAX_EXT};
use crate::dsm::ScoringModel;
use crate::seed::SeedHit;
use crate::types::{Base, Energy};

/// Target range is in the physical duplex frame; FASTA conversion is the caller's job.
pub(super) struct SeedExtension {
    pub(super) q_range: Range<usize>,
    pub(super) t_range: Range<usize>,
    pub(super) energy: Energy,
}

/// Per-worker extension engine. `MAX_EXT` bounds the buffers; raising it means
/// changing `MAX_EXTENSION`, the clap range on `-l`, and `check_unlimited_fits`.
pub(super) struct ExtensionEngine {
    grid: DpGrid,
    q_buf: [u8; MAX_EXT],
    t_buf: [u8; MAX_EXT],
    engine_left: GotohModel,
    engine_right: GotohModel,
    steps: Option<Vec<TraceOp>>,
    max_window: Option<usize>,
}

impl ExtensionEngine {
    pub(super) fn new(
        max_window: Option<usize>,
        build_alignment: bool,
        model: &ScoringModel,
    ) -> Self {
        Self {
            grid: DpGrid::new(max_window.unwrap_or(MAX_EXT)),
            q_buf: [0u8; MAX_EXT],
            t_buf: [0u8; MAX_EXT],
            engine_left: Gotoh::new(&model.transpose()),
            engine_right: Gotoh::new(model),
            steps: build_alignment.then(Vec::new),
            max_window: max_window.map(|n| n.max(1)),
        }
    }

    /// Extend both flanks, score the whole duplex. Left flank reversed +
    /// transposed model; right flank forward + canonical model.
    /// Left = lower duplex index: query 5' side, target 3' side.
    pub(super) fn extend_seed(
        &mut self,
        query: &[Base],
        target: &[Base],
        seed: &SeedHit,
    ) -> SeedExtension {
        let qr = seed.query_range();
        let tr = seed.target_range();

        let q_left = &query[..=qr.start];
        let t_left = &target[..=tr.start];
        let cap = self.max_window.unwrap_or(q_left.len());
        let q_win = fill(&mut self.q_buf, q_left.iter().rev().copied(), cap);
        let t_win = fill(&mut self.t_buf, t_left.iter().rev().copied(), cap);
        trace!("[left] q_len={} t_len={}", q_win.len(), t_win.len());
        let left = self.engine_left.extend(q_win, t_win, &mut self.grid);
        if let Some(steps) = &mut self.steps {
            steps.clear();
            self.engine_left
                .traceback(q_win, t_win, &self.grid, left, |step| steps.push(step));
            steps.resize(steps.len() + qr.len(), TraceOp::Match);
        }

        let q_right = &query[qr.end - 1..];
        let t_right = &target[tr.end - 1..];
        let cap = self.max_window.unwrap_or(q_right.len());
        let q_win = fill(&mut self.q_buf, q_right.iter().copied(), cap);
        let t_win = fill(&mut self.t_buf, t_right.iter().copied(), cap);
        trace!("[right] q_len={} t_len={}", q_win.len(), t_win.len());
        let right = self.engine_right.extend(q_win, t_win, &mut self.grid);
        if let Some(steps) = &mut self.steps {
            let seed_end = steps.len();
            self.engine_right
                .traceback(q_win, t_win, &self.grid, right, |step| steps.push(step));
            // Left steps arrive 5'->3'; right steps arrive 3'->5'.
            steps[seed_end..].reverse();
        }

        debug_assert!(left.q_idx <= qr.start);
        debug_assert!(left.t_idx <= tr.start);
        let q_range = qr.start - left.q_idx..qr.end + right.q_idx;
        let t_range = tr.start - left.t_idx..tr.end + right.t_idx;
        let model = &self.engine_right.scoring;
        let score = model.ungapped_duplex_score(&query[qr], &target[tr])
            + Energy(left.score)
            + Energy(right.score);
        let energy = model.binding_energy(score);

        SeedExtension {
            q_range,
            t_range,
            energy,
        }
    }

    /// Step buffer is overwritten by the next `extend_seed`; call before continuing.
    pub(super) fn materialize_alignment(
        &self,
        query: &[Base],
        target: &[Base],
        ext: &SeedExtension,
    ) -> Option<Box<[AlignColumn]>> {
        let steps = self.steps.as_deref()?;
        let query = &query[ext.q_range.clone()];
        let target = &target[ext.t_range.clone()];
        Some(alignment::resolve(steps, query, target))
    }
}

/// Caller controls orientation via the iterator (`.rev()` for left flank).
fn fill(dst: &mut [u8], src: impl Iterator<Item = Base>, cap: usize) -> &[u8] {
    let mut len = 0;
    for (slot, base) in dst.iter_mut().zip(src.take(cap)) {
        *slot = base.as_u8();
        len += 1;
    }
    &dst[..len]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::DsmId;

    fn test_model() -> ScoringModel {
        ScoringModel::load(&DsmId::from("t04"), 37, Energy::from_kcal(0.0)).unwrap()
    }

    #[test]
    fn penalty_changes_the_selected_endpoint_on_both_flanks() {
        for left in [false, true] {
            let (q, t, start) = if left {
                (
                    vec![Base::A, Base::G, Base::C],
                    vec![Base::U, Base::C, Base::G],
                    1,
                )
            } else {
                (
                    vec![Base::G, Base::C, Base::A],
                    vec![Base::C, Base::G, Base::U],
                    0,
                )
            };
            let mut table = [[[[-10000; 6]; 6]; 6]; 6];
            for (&qb, &tb) in q.iter().zip(&t) {
                table[0][qb.as_usize()][0][tb.as_usize()] = 0;
                table[qb.as_usize()][0][tb.as_usize()][0] = 0;
            }
            for i in 0..2 {
                table[q[i].as_usize()][q[i + 1].as_usize()][t[i].as_usize()][t[i + 1].as_usize()] =
                    if i == start { 1000 } else { 100 };
            }
            let seed = SeedHit::new(start, 0, start, 2, crate::Strand::Forward);
            for penalty in [49, 50, 51] {
                let model = ScoringModel::new(&table, Energy(1000), Energy(penalty));
                let mut engine = ExtensionEngine::new(Some(3), true, &model);
                let hit = engine.extend_seed(&q, &t, &seed);
                let extended = hit.q_range.len() == 3;
                if penalty != 50 {
                    assert_eq!(extended, penalty == 49);
                }
                assert_eq!(
                    hit.energy,
                    Energy(if extended {
                        6 * penalty - 100
                    } else {
                        4 * penalty
                    })
                );
                assert_eq!(hit.q_range, hit.t_range);
                assert_eq!(
                    engine.materialize_alignment(&q, &t, &hit).unwrap().len(),
                    hit.q_range.len()
                );
            }
        }
    }

    #[test]
    fn anchored_duplexes_have_explicit_columns_and_energy() {
        for (q_row, t_row, seed_column) in [
            ("ACGUAC", "UGCAUG", 0),
            ("ACGUAC", "UGCAUG", 4),
            ("ACGUAC", "UGCAUG", 2),
            ("ACGUAC", "UGAAUG", 0),
            ("ACAGUAC", "UG-CAUG", 0),
            ("ACAGUAC", "UG-CAUG", 5),
            ("AC-GUAC", "UGACAUG", 0),
            ("AC-GUAC", "UGACAUG", 5),
        ] {
            let qcols: Vec<_> = q_row.chars().map(|b| Base::try_from(b).unwrap()).collect();
            let tcols: Vec<_> = t_row.chars().map(|b| Base::try_from(b).unwrap()).collect();
            let query: Vec<_> = qcols.iter().copied().filter(|b| *b != Base::Gap).collect();
            let target: Vec<_> = tcols.iter().copied().filter(|b| *b != Base::Gap).collect();
            let mut table = [[[[-100000; 6]; 6]; 6]; 6];
            for (&q, &t) in qcols.iter().zip(&tcols) {
                if q != Base::Gap && t != Base::Gap {
                    table[0][q.as_usize()][0][t.as_usize()] = 11;
                    table[q.as_usize()][0][t.as_usize()][0] = 17;
                }
            }
            for (q, t) in qcols.windows(2).zip(tcols.windows(2)) {
                table[q[0].as_usize()][q[1].as_usize()][t[0].as_usize()][t[1].as_usize()] =
                    100 + q[0].as_usize() as i32 * 10 + t[1].as_usize() as i32;
            }
            let score: i32 = qcols
                .windows(2)
                .zip(tcols.windows(2))
                .map(|(q, t)| 100 + q[0].as_usize() as i32 * 10 + t[1].as_usize() as i32)
                .sum();
            let qs = qcols[..seed_column]
                .iter()
                .filter(|b| **b != Base::Gap)
                .count();
            let ts = tcols[..seed_column]
                .iter()
                .filter(|b| **b != Base::Gap)
                .count();
            assert!(qcols[seed_column..seed_column + 2]
                .iter()
                .all(|b| *b != Base::Gap));
            assert!(tcols[seed_column..seed_column + 2]
                .iter()
                .all(|b| *b != Base::Gap));
            let seed = SeedHit::new(qs, 0, ts, 2, crate::Strand::Forward);
            for penalty in [0, 7] {
                let model = ScoringModel::new(&table, Energy(1000), Energy(penalty));
                let mut engine = ExtensionEngine::new(Some(20), true, &model);
                let hit = engine.extend_seed(&query, &target, &seed);
                assert_eq!(
                    hit.q_range,
                    0..query.len(),
                    "{q_row}/{t_row} seed={seed_column}"
                );
                assert_eq!(hit.t_range, 0..target.len());
                let extension_penalty = penalty * (query.len() + target.len()) as i32;
                assert_eq!(
                    hit.energy,
                    Energy(1000 - score - 11 - 17 + extension_penalty)
                );
                let alignment = engine.materialize_alignment(&query, &target, &hit).unwrap();
                assert_eq!(
                    alignment.iter().map(|c| c.query()).collect::<Vec<_>>(),
                    qcols
                );
                assert_eq!(
                    alignment.iter().map(|c| c.target()).collect::<Vec<_>>(),
                    tcols
                );

                let mut seed_only = ExtensionEngine::new(Some(0), true, &model);
                let hit = seed_only.extend_seed(&query, &target, &seed);
                assert_eq!(hit.q_range, qs..qs + 2);
                assert_eq!(hit.t_range, ts..ts + 2);
                let seed_score =
                    100 + query[qs].as_usize() as i32 * 10 + target[ts + 1].as_usize() as i32;
                assert_eq!(
                    hit.energy,
                    Energy(1000 - seed_score - 11 - 17 + 4 * penalty)
                );
            }
        }
    }

    #[test]
    fn fill_walks_query_and_physical_target_together() {
        let query = [Base::Gap, Base::A, Base::C, Base::G, Base::U];
        let target = [Base::N, Base::U, Base::G, Base::C, Base::A, Base::Gap];
        let expected = [Base::G, Base::C, Base::A].map(Base::as_u8);

        let mut buf = [u8::MAX; 3];
        assert_eq!(
            fill(&mut buf, query[..4].iter().rev().copied(), 3),
            expected
        );
        assert_eq!(
            fill(&mut buf, target[..5].iter().rev().copied(), 3),
            [Base::A, Base::C, Base::G].map(Base::as_u8)
        );
        assert_eq!(
            fill(&mut buf, query[1..].iter().copied(), 3),
            [Base::A, Base::C, Base::G].map(Base::as_u8)
        );
        assert_eq!(fill(&mut buf, target[2..].iter().copied(), 3), expected);
    }

    #[test]
    fn unlimited_window_follows_the_query_not_the_target() {
        let query = [Base::A; 5];
        let long = [Base::A; MAX_EXT + 10];
        let mut buf = [0u8; MAX_EXT];

        // Unlimited: cap follows the query length
        let cap = query.len();
        assert_eq!(fill(&mut buf, query.iter().copied(), cap).len(), 5);
        assert_eq!(fill(&mut buf, long.iter().copied(), cap).len(), 5);

        // Fixed: cap is the window, independent of available sequence
        let cap = 20;
        assert_eq!(fill(&mut buf, query.iter().copied(), cap).len(), 5);
        assert_eq!(fill(&mut buf, long.iter().copied(), cap).len(), 20);

        // Buffer is the ceiling even when cap exceeds it
        for cap in [usize::MAX, MAX_EXT + 1, 1] {
            assert!(fill(&mut buf, long.iter().copied(), cap).len() <= MAX_EXT);
        }
    }

    #[test]
    fn unlimited_extension_matches_a_window_fixed_to_the_query_flank() {
        let model = test_model();
        let query = [Base::A, Base::U, Base::G, Base::C];
        let target = [Base::U; 200];
        let seed = SeedHit::new(1, 0, 1, 2, crate::Strand::Forward);

        let mut unlimited = ExtensionEngine::new(None, true, &model);
        let mut fixed = ExtensionEngine::new(Some(query.len()), true, &model);

        let free = unlimited.extend_seed(&query, &target, &seed);
        let pinned = fixed.extend_seed(&query, &target, &seed);
        assert_eq!(
            (free.energy, free.q_range.clone(), free.t_range.clone()),
            (
                pinned.energy,
                pinned.q_range.clone(),
                pinned.t_range.clone()
            ),
            "unlimited must equal a window fixed to the query flank"
        );
        assert!(
            free.t_range.len() <= query.len(),
            "target extent {} escaped the {}-symbol query window",
            free.t_range.len(),
            query.len()
        );
        let free_cols = unlimited
            .materialize_alignment(&query, &target, &free)
            .unwrap();
        let pinned_cols = fixed
            .materialize_alignment(&query, &target, &pinned)
            .unwrap();
        assert_eq!(free_cols.len(), pinned_cols.len());
    }
}
