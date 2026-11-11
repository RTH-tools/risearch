//! Gotoh 3-state DP recurrence.
//!
//! This module implements the core DP engine for directional extension and
//! obtains every transition score through [`GotohScoring`]. Its three states
//! differ only in coordinate consumption:
//! - `M`: consume one symbol from `q` and `t`
//! - `GapQ`: consume one symbol from `q`
//! - `GapT`: consume one symbol from `t`
//!
//! `M` covers every diagonal transition. The caller interprets the symbols and
//! the recovered state path.
//!
//! # Index safety
//!
//! All index arguments to point-lookup and slice-accessor methods must be in
//! the valid symbol index range of the scoring source. Callers hand
//! `Gotoh::extend()` the window already materialized as dense symbol indices.

use crate::dp::scoring::GotohScoring;
use log::trace;

use super::{is_valid_score, BestScore, DpGrid, TraceOp, MAX_EXT, NEG_INF};

#[inline(always)]
fn is_transition(val: i32, pred: i32, score: i32) -> bool {
    is_valid_score(pred) && val == pred + score
}

/// A Gotoh 3-state DP engine.
///
/// Direction-agnostic: each instance corresponds to one orientation of the
/// underlying scoring source.
pub struct Gotoh<S: GotohScoring> {
    pub(crate) scoring: S,
}

impl<S: GotohScoring> Gotoh<S> {
    /// Create a new Gotoh engine for the given scoring source.
    pub fn new(scoring: &S) -> Self
    where
        S: Clone,
    {
        Self {
            scoring: scoring.clone(),
        }
    }

    #[cfg_attr(feature = "prof", inline(never))]
    /// DP forward pass over one extension window, reusing the caller's grid.
    ///
    /// `q` and `t` hold rank-indexed symbols in DP order: position 0 is the
    /// anchor column at the seed boundary. Neither may exceed [`MAX_EXT`] — the
    /// NEG_INF drift proof is only valid within that bound.
    pub fn extend(&self, q: &[u8], t: &[u8], grid: &mut DpGrid) -> BestScore {
        let (q_len, t_len) = (q.len(), t.len());

        debug_assert!(
            q_len <= MAX_EXT && t_len <= MAX_EXT,
            "DP window exceeds MAX_EXT: {q_len}x{t_len}"
        );

        // A one-position window scores only its anchor, and an empty one has no
        // anchor at all. Taking that branch first means the raw reads below are
        // covered by a check the hot path already pays for.
        if q_len <= 1 || t_len <= 1 {
            return match (q.first(), t.first()) {
                (Some(&q0), Some(&t0)) => BestScore::new(self.scoring.boundary(q0, t0)),
                _ => BestScore::new(NEG_INF),
            };
        }

        let q_ptr = q.as_ptr();
        let t_ptr = t.as_ptr();

        // SAFETY: both windows hold at least 2 symbols per the check above.
        let q0 = unsafe { *q_ptr };
        let t0 = unsafe { *t_ptr };
        let mut best = BestScore::new(self.scoring.boundary(q0, t0));

        grid.resize(t_len + 1, q_len + 1);

        let has_main_region = self.init_frontier(q_ptr, t_ptr, grid, q_len, t_len, &mut best);
        if has_main_region {
            self.main_loop(q_ptr, t_ptr, grid, q_len, t_len, &mut best);
        }

        trace!(
            "dp result: score={} q_idx={} t_idx={}",
            best.score,
            best.q_idx,
            best.t_idx,
        );

        best
    }

    /// Walk the filled grid backward from `best` to recover the state-transition
    /// path, emitting each step to the closure, endpoint first.
    pub(crate) fn traceback(
        &self,
        q: &[u8],
        t: &[u8],
        grid: &DpGrid,
        best: BestScore,
        mut emit: impl FnMut(TraceOp),
    ) {
        let (mut i, mut j) = (best.q_idx, best.t_idx);
        let mut state = TraceOp::Match;

        while i > 0 || j > 0 {
            let next = match state {
                TraceOp::Match if i > 0 && j > 0 => {
                    let c = grid.get(i, j);
                    let diag = grid.get(i - 1, j - 1);

                    let qi_prev = q[i - 1];
                    let qi = q[i];
                    let tj_prev = t[j - 1];
                    let tj = t[j];

                    emit(TraceOp::Match);

                    let r#match = self.scoring.r#match(qi_prev, qi, tj_prev, tj);
                    let close_query_gap = self.scoring.close_query_gap(qi_prev, qi, tj);
                    let close_target_gap = self.scoring.close_target_gap(qi, tj_prev, tj);

                    i -= 1;
                    j -= 1;

                    if is_transition(c.m, diag.m, r#match) {
                        TraceOp::Match
                    } else if is_transition(c.m, diag.gap_q, close_query_gap) {
                        TraceOp::GapQ
                    } else if is_transition(c.m, diag.gap_t, close_target_gap) {
                        TraceOp::GapT
                    } else {
                        debug_assert!(
                            false,
                            "traceback: no valid predecessor at ({i},{j}) in Match state"
                        );
                        break;
                    }
                }
                TraceOp::GapQ if i > 0 => {
                    let c = grid.get(i, j);
                    let up = grid.get(i - 1, j);

                    let qi_prev = q[i - 1];
                    let qi = q[i];
                    let tj = t[j];

                    emit(TraceOp::GapQ);

                    let open_query_gap = self.scoring.open_query_gap(qi_prev, qi, tj);
                    let extend_query_gap = self.scoring.extend_query_gap(qi_prev, qi);

                    i -= 1;

                    if is_transition(c.gap_q, up.m, open_query_gap) {
                        TraceOp::Match
                    } else if is_transition(c.gap_q, up.gap_q, extend_query_gap) {
                        TraceOp::GapQ
                    } else {
                        debug_assert!(
                            false,
                            "traceback: no valid predecessor at ({i},{j}) in GapQ state"
                        );
                        break;
                    }
                }
                TraceOp::GapT if j > 0 => {
                    let c = grid.get(i, j);
                    let left = grid.get(i, j - 1);

                    let qi = q[i];
                    let tj_prev = t[j - 1];
                    let tj = t[j];

                    emit(TraceOp::GapT);

                    let open_target_gap = self.scoring.open_target_gap(qi, tj_prev, tj);
                    let extend_target_gap = self.scoring.extend_target_gap(tj_prev, tj);

                    j -= 1;

                    if is_transition(c.gap_t, left.m, open_target_gap) {
                        TraceOp::Match
                    } else if is_transition(c.gap_t, left.gap_t, extend_target_gap) {
                        TraceOp::GapT
                    } else {
                        debug_assert!(
                            false,
                            "traceback: no valid predecessor at ({i},{j}) in GapT state"
                        );
                        break;
                    }
                }
                _ => break,
            };
            state = next;
        }
    }
}

#[cfg(test)]
mod tests {
    // End-to-end tests for the DP algorithm are located in the /tests directory.

    use super::is_transition;
    use crate::dp::NEG_INF;

    /// `is_transition` must reject a predecessor whose state is unreachable
    /// (NEG_INF sentinel), even when the arithmetic identity coincidentally
    /// holds. Without the `is_valid_score(pred)` guard, traceback could
    /// follow a phantom predecessor through the NEG_INF region.
    ///
    /// Catches the `&& → ||` mutation surfaced by cargo-mutants on this fn.
    #[test]
    fn is_transition_rejects_invalid_pred_with_coincidental_arithmetic() {
        let score = 100;
        let val = NEG_INF + score;
        // pred is NEG_INF (invalid). val happens to equal pred + score.
        // Under `&& `: returns false (invalid pred fails the guard).
        // Under `||`: would return true (arithmetic check matches).
        assert!(
            !is_transition(val, NEG_INF, score),
            "is_transition must reject invalid predecessor regardless of arithmetic"
        );
    }

    #[test]
    fn is_transition_accepts_valid_pred_with_matching_arithmetic() {
        // Normal happy path: pred reachable, arithmetic matches → true.
        assert!(is_transition(150, 50, 100));
    }

    #[test]
    fn is_transition_rejects_valid_pred_with_mismatched_arithmetic() {
        // Pred reachable but arithmetic doesn't match → false.
        // Catches mutations that drop the `==` check.
        assert!(!is_transition(151, 50, 100));
    }

    #[test]
    fn is_transition_rejects_invalid_pred_with_mismatched_arithmetic() {
        // Both conditions fail → false. Baseline correctness check.
        assert!(!is_transition(0, NEG_INF, 100));
    }

    /// Small safety-focused domain for the interpreter: early returns, both tiny
    /// frontiers, the 3x3 corner, the pointer-based interior, and shrink/reuse.
    /// Execution alone is the check; these shapes are compared against the
    /// reference by the exhaustive and generated tests.
    #[cfg(miri)]
    #[test]
    fn frontier_and_grid_reuse_safety() {
        use crate::dp::{gotoh::Gotoh, DpGrid};
        use crate::dsm::{DsmTable, ScoringModel};
        use crate::types::Energy;

        let table: DsmTable = std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                std::array::from_fn(|c| {
                    std::array::from_fn(|d| {
                        ((a * 71 + b * 47 + c * 19 + d * 11 + 7) % 29) as i32 - 18
                    })
                })
            })
        });
        let engine = Gotoh::new(&ScoringModel::new(&table, Energy(0), Energy(0)));
        let query = [1, 3, 2, 5, 4, 1];
        let target = [5, 2, 3, 1, 4, 3];
        let mut grid = DpGrid::new(1);
        for (qlen, tlen) in [
            (0, 0),
            (0, 1),
            (1, 0),
            (1, 1),
            (2, 6),
            (6, 2),
            (3, 3),
            (3, 6),
            (6, 3),
            (6, 6),
            (2, 2),
            (4, 5),
            (1, 6),
        ] {
            let (q, t) = (&query[..qlen], &target[..tlen]);
            let best = engine.extend(q, t, &mut grid);
            if !q.is_empty() && !t.is_empty() {
                engine.traceback(q, t, &grid, best, |_, _, _| {});
            }
        }
    }
}
