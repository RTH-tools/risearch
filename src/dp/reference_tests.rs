//! Independent checks of anchored alignment scores and traceback consumption.
//!
//! The production DP is written for speed. Scores come from a row profile of raw
//! pointers, the matrix is split into a hand-written frontier and a pointer-walking
//! interior, the grid is reused without clearing, and unreachable cells carry a
//! numeric `NEG_INF` rather than an option. A wrong index, a stale cell or a
//! sentinel that drifted into valid range all yield numbers that still look
//! reasonable, so the output has to be checked against something else.
//!
//! `naive_gotoh` is the single reference implementation: three explicit scalar
//! recurrences over a fully initialized matrix, predecessor tracking, and
//! pointer-chase traceback. Exhaustive sequence pairs compare production with it.
//! It does not use production scoring accessors, frontier initialization, pointer
//! arithmetic, or sentinels. The only shared input is the raw `DsmTable` tensor.
//!
//! Inputs are canonical ranks, matching the DP boundary. A result tuple is
//! `(score_including_terminal, query_offset, target_offset)`: offsets are zero-based
//! from the paired anchor, so they also count consumed bases beyond that anchor.
//! Ranks are `Base` discriminants: 0 is a gap, 1 to 5 are A, C, G, N, U. Rank 0
//! never occurs in an input window, only as a tensor coordinate naming the gapped
//! side of a column. Scores are raw tensor units, with no energy conversion.
//!
//! Each sequence pair checks three things: every grid cell production leaves
//! behind, state by state; the reported best score and the endpoint it names; and
//! the trace, compared directly against the reference. Tie-breaking matches
//! production exactly: predecessor priority PAIRED → QUERY_ONLY → TARGET_ONLY,
//! and first-seen optimal endpoint in row-major order.

use super::{gotoh::Gotoh, DpGrid, TraceOp, MAX_EXT};
use crate::dsm::{DsmTable, ScoringModel};
use crate::types::{Base, Energy};

// State indices into the reference's `[Option<i64>; 3]` cell. The order follows
// `DpCell`'s fields (`m`, `gap_q`, `gap_t`), which is what lets
// `assert_grid_matches_reference` index both representations with one counter.
const PAIRED: usize = 0;
const QUERY_ONLY: usize = 1;
const TARGET_ONLY: usize = 2;

/// Transition score between two `(query_rank, target_rank)` columns.
/// Rank 0 is a gap; `(0, 0)` as `next` gives the terminal `table[q][0][t][0]`.
fn score_column_transition(
    table: &DsmTable,
    previous: (usize, usize),
    next: (usize, usize),
) -> i64 {
    let (previous_query, previous_target) = previous;
    let (next_query, next_target) = next;
    i64::from(table[previous_query][next_query][previous_target][next_target])
}

/// Add a transition score to a reachable predecessor; `None` stays `None`.
/// Production uses a `NEG_INF` sentinel instead — a disagreement here catches
/// sentinel drift.
fn add_transition_score(predecessor: Option<i64>, transition: i32) -> Option<i64> {
    predecessor.map(|score| score + i64::from(transition))
}

/// Naive, anchored Gotoh alignment over the raw transition tensor.
///
/// Three recurrences, evaluated once at each `(i, j)`:
///
/// - `PAIRED`: consume query and target; arrive diagonally from any state.
/// - `QUERY_ONLY`: consume query; arrive from above, opening or extending a gap.
/// - `TARGET_ONLY`: consume target; arrive from the left, opening or extending a gap.
///
/// `(i, j)` names the last consumed position on each input, not prefix lengths.
/// The inputs include an already paired anchor at `(0, 0)`, whose score is zero.
/// A gap keeps the other input's position fixed. Opposite gap states cannot
/// follow each other directly: a paired column must separate them.
///
/// The tensor is indexed `[previous_query][next_query][previous_target][next_target]`.
/// Inputs contain non-gap ranks 1 through 5 (A/C/G/N/U), in extension order.
/// Rank 0 is a gap. Each recurrence below spells out these four column symbols;
/// it uses neither the production scoring adapter nor production DP helpers.
///
/// Every cell starts unreachable. There is no local-alignment reset to zero:
/// even a negative-scoring path must remain connected to the anchor. Empty
/// inputs have no anchor. At the top/left edges, the missing predecessors simply
/// leave states unreachable; there is no separate frontier implementation.
/// Production hand-codes rows 0-2 and columns 0-2, with another path for windows
/// shorter than three, and reaches its pointer loop only at `(3, 3)`. Here those
/// cells come out of the same two `if` guards as every other cell.
///
/// Predecessor pointers are recorded during the forward pass. Traceback is a
/// pointer chase from the best endpoint — no score recomputation, no transition
/// lookups. Tie-breaking matches production: candidates are checked in priority
/// order (PAIRED → QUERY_ONLY → TARGET_ONLY) and the first to beat the current
/// best wins (`>`, not `>=`). The exhaustive test proves this across all short
/// inputs.
///
/// Returns three facts: the full grid of `Option<i64>` scores, the best
/// terminal-inclusive score with its coordinates, and the traceback path. The
/// wider integers keep expected arithmetic independent of production's i32 scores
/// and negative-infinity sentinel.
///
/// Grid scores exclude terminals. The returned best score includes one: for
/// each reachable paired cell the terminal `table[q][0][t][0]` is added, and
/// the maximum over all such cells is kept. Empty inputs have no anchor and
/// yield `(i64::MIN, 0, 0)` with an empty trace.
type NaiveGrid = Vec<Vec<[Option<i64>; 3]>>;

#[allow(clippy::type_complexity)]
fn naive_gotoh(
    query: &[u8],
    target: &[u8],
    table: &DsmTable,
) -> (NaiveGrid, (i64, usize, usize), Vec<TraceOp>) {
    let mut scores = vec![vec![[None; 3]; target.len()]; query.len()];
    let mut preds = vec![vec![[None::<usize>; 3]; target.len()]; query.len()];

    if query.is_empty() || target.is_empty() {
        return (scores, (i64::MIN, 0, 0), Vec::new());
    }
    scores[0][0][PAIRED] = Some(0);

    for i in 0..query.len() {
        for j in 0..target.len() {
            let q = query[i] as usize;
            let t = target[j] as usize;

            // A paired column consumes both inputs. The previous column may
            // have been paired, query-only, or target-only, respectively.
            //
            // All three predecessors sit at (i-1, j-1) and differ only in what
            // their own column looked like, which the tensor's "previous"
            // coordinates spell out:
            //   from PAIRED       (q_prev, t_prev) -> (q, t)
            //   from QUERY_ONLY   (q_prev, gap)    -> (q, t)   target side 0
            //   from TARGET_ONLY  (gap,    t_prev) -> (q, t)   query  side 0
            // Both gap states can close into a pair, so a paired column is the
            // only place the two gap sides meet.
            if i > 0 && j > 0 {
                let q_prev = query[i - 1] as usize;
                let t_prev = target[j - 1] as usize;
                let diag = scores[i - 1][j - 1];
                let candidates = [
                    (
                        PAIRED,
                        add_transition_score(diag[PAIRED], table[q_prev][q][t_prev][t]),
                    ),
                    (
                        QUERY_ONLY,
                        add_transition_score(diag[QUERY_ONLY], table[q_prev][q][0][t]),
                    ),
                    (
                        TARGET_ONLY,
                        add_transition_score(diag[TARGET_ONLY], table[0][q][t_prev][t]),
                    ),
                ];
                for &(src, val) in &candidates {
                    if val > scores[i][j][PAIRED] {
                        scores[i][j][PAIRED] = val;
                        preds[i][j][PAIRED] = Some(src);
                    }
                }
            }

            // A query-only column is (q, gap). Open from a paired column
            // (q_prev, t), or extend (q_prev, gap). Never arrive from TARGET_ONLY.
            //
            // `j` does not move, so the predecessor's target symbol is the same `t`
            // this column skips over: open reads table[q_prev][q][t][0], not
            // table[q_prev][q][t_prev][0]. Leaving out a TARGET_ONLY predecessor is
            // the gap-adjacency rule: (gap, t) followed by (q, gap) aligns nothing,
            // and is written as one paired column or a longer detour instead.
            if i > 0 {
                let q_prev = query[i - 1] as usize;
                let above = scores[i - 1][j];
                let candidates = [
                    (
                        PAIRED,
                        add_transition_score(above[PAIRED], table[q_prev][q][t][0]),
                    ),
                    (
                        QUERY_ONLY,
                        add_transition_score(above[QUERY_ONLY], table[q_prev][q][0][0]),
                    ),
                ];
                for &(src, val) in &candidates {
                    if val > scores[i][j][QUERY_ONLY] {
                        scores[i][j][QUERY_ONLY] = val;
                        preds[i][j][QUERY_ONLY] = Some(src);
                    }
                }
            }

            // A target-only column is (gap, t). Open from a paired column
            // (q, t_prev), or extend (gap, t_prev). Never arrive from QUERY_ONLY.
            //
            // Mirror of the block above: `i` is fixed, so the predecessor's query
            // symbol is this cell's own `q`, and extension carries gap on the query
            // side of both columns.
            if j > 0 {
                let t_prev = target[j - 1] as usize;
                let left = scores[i][j - 1];
                let candidates = [
                    (
                        PAIRED,
                        add_transition_score(left[PAIRED], table[q][0][t_prev][t]),
                    ),
                    (
                        TARGET_ONLY,
                        add_transition_score(left[TARGET_ONLY], table[0][0][t_prev][t]),
                    ),
                ];
                for &(src, val) in &candidates {
                    if val > scores[i][j][TARGET_ONLY] {
                        scores[i][j][TARGET_ONLY] = val;
                        preds[i][j][TARGET_ONLY] = Some(src);
                    }
                }
            }
        }
    }

    // Best terminal-inclusive score over all reachable paired cells.
    let mut best = (i64::MIN, 0, 0);
    for (i, row) in scores.iter().enumerate() {
        for (j, states) in row.iter().enumerate() {
            if let Some(score) = states[PAIRED] {
                let with_terminal = score
                    + score_column_transition(
                        table,
                        (query[i] as usize, target[j] as usize),
                        (0, 0),
                    );
                if with_terminal > best.0 {
                    best = (with_terminal, i, j);
                }
            }
        }
    }

    // Traceback: follow predecessor pointers from the best endpoint.
    let mut trace = Vec::new();
    let (mut i, mut j) = (best.1, best.2);
    let mut state = PAIRED;
    while i > 0 || j > 0 {
        let src = match preds[i][j][state] {
            Some(s) => s,
            None => break,
        };
        match state {
            PAIRED => {
                trace.push(TraceOp::Match);
                i -= 1;
                j -= 1;
            }
            QUERY_ONLY => {
                trace.push(TraceOp::GapQ);
                i -= 1;
            }
            TARGET_ONLY => {
                trace.push(TraceOp::GapT);
                j -= 1;
            }
            _ => unreachable!(),
        }
        state = src;
    }
    debug_assert!(
        (best.1 > 0 && best.2 > 0) || trace.is_empty(),
        "edge endpoint must produce empty trace"
    );

    (scores, best, trace)
}

/// Every canonical input word, including empty, ordered by length then symbols.
/// Gap is an alignment operation, not an input symbol; N is an input symbol.
///
/// Built breadth-first: `previous_length` holds the index range of the words of the
/// current length, and each round appends every one-symbol extension of them. The
/// empty word is there because `extend` has to handle a zero-length window. N is a
/// real rank with its own row in the tensor, so index arithmetic that works only for
/// the four unambiguous bases would pass without it.
///
/// The result has `5^0 + 5^1 + ... + 5^max_length` entries. The caller asserts that
/// count, so a change in enumeration cannot quietly shrink coverage.
fn all_sequence_windows(max_length: usize) -> Vec<Vec<u8>> {
    let alphabet = [Base::A, Base::C, Base::G, Base::N, Base::U];
    let mut windows = vec![Vec::new()];
    let mut previous_length = 0..1;
    for _ in 0..max_length {
        let start = windows.len();
        for index in previous_length {
            for base in alphabet {
                let mut next = windows[index].clone();
                next.push(base.as_u8());
                windows.push(next);
            }
        }
        previous_length = start..windows.len();
    }
    windows
}

/// Synthetic tensor with asymmetric mixed-sign entries (`-18..=10`).
/// Different odd multipliers per axis — swapping any two coordinates changes
/// almost every cell, catching transpositions a symmetric table would hide.
fn synthetic_table() -> DsmTable {
    std::array::from_fn(|a| {
        std::array::from_fn(|b| {
            std::array::from_fn(|c| {
                std::array::from_fn(|d| ((a * 71 + b * 47 + c * 19 + d * 11 + 7) % 29) as i32 - 18)
            })
        })
    })
}

/// Length 3 (debug): 24K pairs, reaches the first interior cell `(3, 3)`.
/// Length 5 (release): 15M pairs. Both cross the frontier boundary.
const MAX_LEN: usize = if cfg!(debug_assertions) { 3 } else { 5 };

/// Every word pair up to `MAX_LEN` against the naive reference.
/// Covers every frontier code path and asymmetric length combinations.
/// `Energy(0)` makes `ScoringModel` a pass-through of the raw tensor.
#[test]
fn exhaustive_sequences_match_the_naive_reference() {
    let windows = all_sequence_windows(MAX_LEN);
    assert_eq!(windows.len(), (5usize.pow(MAX_LEN as u32 + 1) - 1) / 4);
    let table = synthetic_table();
    let engine = Gotoh::new(&ScoringModel::new(&table, Energy(0), Energy(0)));
    let mut grid = DpGrid::new(1);
    for query in &windows {
        for target in &windows {
            // The oracle allocates a fresh, fully initialized matrix. Only
            // production reuses its grid, so stale cells cannot be shared.
            let (reference, expected, trace) = naive_gotoh(query, target, &table);
            assert_extension_matches(
                query, target, &engine, &mut grid, &reference, expected, &trace,
            );
        }
    }
}

/// Cell-by-cell grid comparison. Skips gap states on the last row/column —
/// production omits those (dead-end gaps that can never close into a pair),
/// and the reused grid may hold stale data there.
///
/// Unreachable cells are not pinned to a specific sentinel; the assertion only
/// checks they stayed below `2 * MAX_EXT * MAX_TRANSITION_SCORE`, the largest
/// score any bounded path could reach.
fn assert_grid_matches_reference(
    query: &[u8],
    target: &[u8],
    reference: &[Vec<[Option<i64>; 3]>],
    grid: &DpGrid,
) {
    if query.len() < 2 || target.len() < 2 {
        return;
    }
    for (i, row) in reference.iter().enumerate() {
        for (j, states) in row.iter().enumerate() {
            let actual = grid.get(i, j);
            // `DpCell` lists its fields in PAIRED/QUERY_ONLY/TARGET_ONLY order, so
            // `state` indexes both representations.
            for (state, value) in [actual.m, actual.gap_q, actual.gap_t]
                .into_iter()
                .enumerate()
            {
                if state != PAIRED && (i + 1 == query.len() || j + 1 == target.len()) {
                    continue;
                }
                if let Some(expected) = states[state] {
                    assert_eq!(
                        i64::from(value),
                        expected,
                        "state {state} at ({i},{j}), q={query:?} t={target:?}"
                    );
                } else {
                    // Do not prescribe the numeric sentinel or its drift. It must
                    // remain below every score of a reachable bounded path.
                    assert!(
                        i64::from(value) < -(2 * MAX_EXT as i64 * super::MAX_TRANSITION_SCORE),
                        "unreachable state {state} at ({i},{j}) acquired score {value}"
                    );
                }
            }
        }
    }
}

/// Assert that production's best score and endpoint match the naive reference.
fn assert_best_score_matches_reference(
    best: &super::BestScore,
    expected: (i64, usize, usize),
    query: &[u8],
    target: &[u8],
) {
    let (expected_score, expected_qi, expected_tj) = expected;
    assert_eq!(
        i64::from(best.score),
        expected_score,
        "score: q={query:?} t={target:?}"
    );
    assert_eq!(
        (best.q_idx, best.t_idx),
        (expected_qi, expected_tj),
        "endpoint: q={query:?} t={target:?}"
    );
}

/// Assert that production's traceback matches the naive reference.
fn assert_trace_matches_reference(
    engine: &Gotoh<ScoringModel>,
    query: &[u8],
    target: &[u8],
    grid: &DpGrid,
    best: &super::BestScore,
    expected_trace: &[TraceOp],
) {
    let mut prod_trace = Vec::new();
    engine.traceback(query, target, grid, *best, |op| prod_trace.push(op));
    assert_eq!(
        prod_trace.as_slice(),
        expected_trace,
        "trace: q={query:?} t={target:?}"
    );
}

/// Run production and compare all three facts against the naive reference:
/// grid, best score, and trace.
fn assert_extension_matches(
    query: &[u8],
    target: &[u8],
    engine: &Gotoh<ScoringModel>,
    grid: &mut DpGrid,
    reference: &[Vec<[Option<i64>; 3]>],
    expected: (i64, usize, usize),
    expected_trace: &[TraceOp],
) -> super::BestScore {
    let best = engine.extend(query, target, grid);
    assert_grid_matches_reference(query, target, reference, grid);
    if query.is_empty() || target.is_empty() {
        assert_eq!(best.score, super::NEG_INF);
        return best;
    }
    assert_best_score_matches_reference(&best, expected, query, target);
    assert_trace_matches_reference(engine, query, target, grid, &best, expected_trace);
    best
}

/// Domain boundaries: `MAX_EXT`-length windows and extreme transition scores.
///
/// Lengths cover `MAX_EXT`, one below it, the three small windows with their own
/// code paths, and zero. Two scales bracket the arithmetic: `scale == 1` for
/// readable off-by-one diagnostics, and `scale == MAX_TRANSITION_SCORE` where a
/// 256-long window accumulates close to `NEG_INF`'s headroom limit.
#[test]
fn boundary_dimensions_and_extreme_scores_match_wide_reference() {
    let mut grid = DpGrid::new(1);
    for scale in [1, super::MAX_TRANSITION_SCORE as i32] {
        let mut table = [[[[0; 6]; 6]; 6]; 6];
        for (a, slab) in table.iter_mut().enumerate() {
            for (b, plane) in slab.iter_mut().enumerate() {
                for (c, row) in plane.iter_mut().enumerate() {
                    for (d, cell) in row.iter_mut().enumerate() {
                        *cell = ((a * 13 + b * 7 + c * 5 + d) % 3) as i32 * scale - scale;
                    }
                }
            }
        }
        for qlen in [256, 1, 3, 0, 255, 2, 4] {
            for tlen in [0, 256, 2, 255, 1, 4, 3] {
                let query: Vec<_> = (0..qlen).map(|i| 1 + (i % 5) as u8).collect();
                let target: Vec<_> = (0..tlen).map(|i| 5 - (i % 5) as u8).collect();
                let (reference, expected, trace) = naive_gotoh(&query, &target, &table);
                let engine = Gotoh::new(&ScoringModel::new(&table, Energy(0), Energy(0)));
                assert_extension_matches(
                    &query, &target, &engine, &mut grid, &reference, expected, &trace,
                );
            }
        }
    }
}

// The fixed tables above are closed-form — a wrong row-profile coordinate can
// stay consistent across every cell it reaches. Drawing all 1_296 entries (6^4)
// independently breaks that: a wrong lane has no relation to the right one.
// The distribution is 4 parts zero, 4 parts small, 1 part +MAX, 1 part −MAX,
// keeping ties common while mixing in extreme accumulation.
//
// Windows arrive in batches through a single reused grid, exercising reuse
// against shapes proptest picks rather than the fixed ordering above.
// Excluded under Miri; `gotoh::tests::frontier_and_grid_reuse_safety` covers
// the same paths there.
#[cfg(not(miri))]
mod generated {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn generated_tables_and_reused_grids(
            entries in prop::collection::vec(prop_oneof![
                4 => Just(0i32),
                4 => -1000i32..=1000,
                1 => Just(super::super::MAX_TRANSITION_SCORE as i32),
                1 => Just(-(super::super::MAX_TRANSITION_SCORE as i32)),
            ], 1296),
            windows in prop::collection::vec((prop::collection::vec(1u8..=5, 0..20),
                prop::collection::vec(1u8..=5, 0..20)), 1..8),
        ) {
            let mut table = [[[[0; 6]; 6]; 6]; 6];
            for (cell, value) in table.iter_mut().flatten().flatten().flatten().zip(entries) { *cell = value; }
            let mut grid = DpGrid::new(1);
            for (query, target) in windows {
                let (reference, expected, trace) = naive_gotoh(&query, &target, &table);
                let engine = Gotoh::new(&ScoringModel::new(&table, Energy(0), Energy(0)));
                assert_extension_matches(&query, &target, &engine, &mut grid, &reference, expected, &trace);
            }
        }
    }
}
