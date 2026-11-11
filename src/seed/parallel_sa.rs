//! Parallel suffix-array traversal kernel for seed search.

use std::ops::Range;

use crate::index::sa::SuffixIndexView;
use crate::types::{Base, BASE_COUNT};

const LINEAR_PARTITION_CUTOFF: usize = 1024;
const MATCHABLE_BUCKETS: [(usize, Base); 4] = [
    (Base::A.as_usize() - 1, Base::A),
    (Base::C.as_usize() - 1, Base::C),
    (Base::G.as_usize() - 1, Base::G),
    (Base::U.as_usize() - 1, Base::U),
];

/// A raw seed candidate emitted by the parallel SA traversal.
#[derive(Debug, Clone)]
pub(super) struct SeedMatch {
    pub query_interval: Range<usize>,
    pub target_interval: Range<usize>,
    pub seed_len: usize,
}

// Cohesive traversal parameters: query/target views plus the seed-length and
// mismatch bounds that drive the SA descent; grouping would only obscure them.
#[allow(clippy::too_many_arguments)]
pub(super) fn traverse<const WOBBLE: bool, F: FnMut(SeedMatch)>(
    query: SuffixIndexView<'_>,
    target: SuffixIndexView<'_>,
    min_len: usize,
    max_len: usize,
    max_mm: usize,
    min_prefix: usize,
    min_suffix: usize,
    on_match: &mut F,
) {
    if min_len > max_len || query.is_empty() || target.is_empty() {
        return;
    }

    let mut traverser = ParallelSaTraverser {
        query,
        target,
        min_len,
        max_len,
        max_mm,
        min_prefix,
        min_suffix,
        on_match,
    };
    let q_range = 0..traverser.query.len();
    let t_range = 0..traverser.target.len();

    traverser.recurse::<WOBBLE>(&q_range, &t_range, 0, 0, 0);
}

struct ParallelSaTraverser<'a, F: FnMut(SeedMatch)> {
    query: SuffixIndexView<'a>,
    target: SuffixIndexView<'a>,
    min_len: usize,
    max_len: usize,
    max_mm: usize,
    min_prefix: usize,
    min_suffix: usize,
    on_match: &'a mut F,
}

impl<F: FnMut(SeedMatch)> ParallelSaTraverser<'_, F> {
    /// Partition a sorted SA interval by base character at `depth`.
    ///
    /// Returns boundary positions `[A, C, G, N, U, end]` in suffix-array order.
    /// The leading `Gap` bucket and the `N` slot are omitted by traversal.
    #[inline(always)]
    fn partition(
        suffixes: SuffixIndexView<'_>,
        range: &Range<usize>,
        depth: usize,
    ) -> [usize; BASE_COUNT] {
        let (start, end) = (range.start, range.end);
        if start >= end {
            return [start; BASE_COUNT];
        }

        let mut out = [0usize; BASE_COUNT];
        if end - start <= LINEAR_PARTITION_CUTOFF {
            let mut i = start;
            for (slot, out_slot) in out.iter_mut().enumerate().take(BASE_COUNT - 1) {
                let target = (slot + 1) as u8;
                // SAFETY: `i` is inside this SA interval.
                while i < end && unsafe { suffixes.base_unchecked(i, depth).as_u8() } < target {
                    i += 1;
                }
                *out_slot = i;
            }
        } else {
            for (slot, out_slot) in out.iter_mut().enumerate().take(BASE_COUNT - 1) {
                let target = (slot + 1) as u8;
                *out_slot = Self::lower_bound(suffixes, start..end, depth, target);
            }
        }
        out[BASE_COUNT - 1] = end;
        out
    }

    /// Leftmost position where the character at `depth` is at least `target`.
    #[inline(always)]
    fn lower_bound(
        suffixes: SuffixIndexView<'_>,
        mut range: Range<usize>,
        depth: usize,
        target: u8,
    ) -> usize {
        let mut half = (range.end - range.start) >> 1;
        while range.start < range.end {
            let mid = range.start + half;
            // SAFETY: `mid` is inside this SA interval.
            if unsafe { suffixes.base_unchecked(mid, depth).as_u8() } >= target {
                range.end = range.start + half;
            } else {
                range.start += if half != 0 { half } else { 1 };
            }
            half >>= 1;
        }
        range.start
    }

    #[inline(always)]
    fn should_emit(&self, depth: usize, match_streak: usize, mm_count: usize) -> bool {
        // Interrupted seeds must not contain a complete perfect sub-seed.
        // can_mismatch_next prevents crossing such a run internally; this
        // condition also excludes it at the trailing end.
        depth >= self.min_len
            && depth <= self.max_len
            && (mm_count == 0 || (match_streak >= self.min_suffix && match_streak < self.min_len))
    }

    #[inline(always)]
    fn can_mismatch_next(&self, depth: usize, match_streak: usize, mm_count: usize) -> bool {
        let next_depth = depth + 1;
        self.max_mm > 0
            && mm_count < self.max_mm
            && next_depth > self.min_prefix
            && match_streak < self.min_len
            && self.max_len - next_depth >= self.min_suffix
    }
    #[inline(always)]
    fn recurse<const WOBBLE: bool>(
        &mut self,
        query_sa: &Range<usize>,
        target_sa: &Range<usize>,
        depth: usize,
        match_streak: usize,
        mm_count: usize,
    ) {
        if self.should_emit(depth, match_streak, mm_count) {
            (self.on_match)(SeedMatch {
                query_interval: query_sa.start..query_sa.end,
                target_interval: target_sa.start..target_sa.end,
                seed_len: depth,
            });
        }

        if depth >= self.max_len {
            return;
        }

        if query_sa.len() == 1 && target_sa.len() == 1 {
            self.recurse_singleton::<WOBBLE>(
                query_sa.start,
                target_sa.start,
                depth,
                match_streak,
                mm_count,
            );
            return;
        }
        if query_sa.len() == 1 {
            self.recurse_half_singleton::<WOBBLE, true>(
                query_sa.start,
                target_sa,
                depth,
                match_streak,
                mm_count,
            );
            return;
        }
        if target_sa.len() == 1 {
            self.recurse_half_singleton::<WOBBLE, false>(
                target_sa.start,
                query_sa,
                depth,
                match_streak,
                mm_count,
            );
            return;
        }

        let query_partitions = Self::partition(self.query, query_sa, depth);
        let target_partitions = Self::partition(self.target, target_sa, depth);

        if query_partitions[0] == query_sa.end || target_partitions[0] == target_sa.end {
            return;
        }

        let next_depth = depth + 1;
        let ms = match_streak + 1;
        let can_mm = self.can_mismatch_next(depth, match_streak, mm_count);

        for (query_slot, query_base) in MATCHABLE_BUCKETS {
            if query_partitions[query_slot] >= query_partitions[query_slot + 1] {
                continue;
            }

            for (target_slot, target_base) in MATCHABLE_BUCKETS {
                if target_partitions[target_slot] >= target_partitions[target_slot + 1] {
                    continue;
                }

                let q = query_partitions[query_slot]..query_partitions[query_slot + 1];
                let t = target_partitions[target_slot]..target_partitions[target_slot + 1];
                if query_base.pair_type(target_base).is_match(WOBBLE) {
                    self.recurse::<WOBBLE>(&q, &t, next_depth, ms, mm_count);
                } else if can_mm {
                    self.recurse::<WOBBLE>(&q, &t, next_depth, 0, mm_count + 1);
                }
            }
        }
    }

    /// Half-singleton: one SA interval has a single entry, the other has multiple.
    /// `Q_SINGLETON=true` -> query is the single entry, partition target.
    /// `Q_SINGLETON=false` -> target is the single entry, partition query.
    #[inline(always)]
    fn recurse_half_singleton<const WOBBLE: bool, const Q_SINGLETON: bool>(
        &mut self,
        singleton_sa_idx: usize,
        multi_sa: &Range<usize>,
        depth: usize,
        match_streak: usize,
        mm_count: usize,
    ) {
        let (singleton, multi) = if Q_SINGLETON {
            (self.query, self.target)
        } else {
            (self.target, self.query)
        };

        // SAFETY: traversal intervals only contain valid SA indices.
        let singleton_base = unsafe { singleton.base_unchecked(singleton_sa_idx, depth) };
        if !singleton_base.is_matchable() {
            return;
        }

        let multi_partitions = Self::partition(multi, multi_sa, depth);
        let d1 = depth + 1;
        let can_mm = self.can_mismatch_next(depth, match_streak, mm_count);
        let singleton_range = singleton_sa_idx..singleton_sa_idx + 1;
        let ms = match_streak + 1;
        let mm1 = mm_count + 1;

        for (multi_slot, multi_base) in MATCHABLE_BUCKETS {
            if multi_partitions[multi_slot] >= multi_partitions[multi_slot + 1] {
                continue;
            }
            let multi_range = multi_partitions[multi_slot]..multi_partitions[multi_slot + 1];
            let (query_sa, target_sa) = if Q_SINGLETON {
                (&singleton_range, &multi_range)
            } else {
                (&multi_range, &singleton_range)
            };
            if singleton_base.pair_type(multi_base).is_match(WOBBLE) {
                self.recurse::<WOBBLE>(query_sa, target_sa, d1, ms, mm_count);
            } else if can_mm {
                self.recurse::<WOBBLE>(query_sa, target_sa, d1, 0, mm1);
            }
        }
    }

    /// Both-singleton fast path: tight linear scan when both SA intervals have one entry.
    ///
    /// Assumes the caller already handled emission at the current depth, then advances
    /// linearly and emits only for newly reached depths.
    #[inline(always)]
    fn recurse_singleton<const WOBBLE: bool>(
        &mut self,
        query_sa_idx: usize,
        target_sa_idx: usize,
        mut depth: usize,
        mut match_streak: usize,
        mut mm_count: usize,
    ) {
        loop {
            if depth >= self.max_len {
                return;
            }

            let d1 = depth + 1;
            let can_mm = self.can_mismatch_next(depth, match_streak, mm_count);

            // SAFETY: traversal intervals only contain valid SA indices.
            let query_base = unsafe { self.query.base_unchecked(query_sa_idx, depth) };
            let target_base = unsafe { self.target.base_unchecked(target_sa_idx, depth) };
            if !query_base.is_matchable() || !target_base.is_matchable() {
                return;
            }

            if query_base.pair_type(target_base).is_match(WOBBLE) {
                depth = d1;
                match_streak += 1;
            } else if can_mm {
                depth = d1;
                mm_count += 1;
                match_streak = 0;
            } else {
                return;
            }

            if self.should_emit(depth, match_streak, mm_count) {
                (self.on_match)(SeedMatch {
                    query_interval: query_sa_idx..query_sa_idx + 1,
                    target_interval: target_sa_idx..target_sa_idx + 1,
                    seed_len: depth,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::config::SeedConfig;
    use crate::index::sa::SuffixIndex;
    use crate::seq::Sequence;
    use crate::types::Base;

    use super::*;

    #[test]
    fn pruning_predicates_respect_the_suffix_and_mismatch_bounds() {
        let seq = Sequence::from(vec![Base::A, Base::C]);
        let prepared = SuffixIndex::build_for_seed(&seq, 1).unwrap();
        let view = prepared.view();
        let mut sink = |_: SeedMatch| {};
        let traverser = ParallelSaTraverser {
            query: view,
            target: view,
            min_len: 4,
            max_len: 8,
            max_mm: 2,
            min_prefix: 1,
            min_suffix: 3,
            on_match: &mut sink,
        };

        assert!(traverser.can_mismatch_next(2, 0, 0));
        // max_len - next_depth = 1, short of min_suffix.
        assert!(!traverser.can_mismatch_next(6, 0, 0));
    }

    #[test]
    fn partition_splits_by_base() {
        let seq_bases = vec![Base::A, Base::G, Base::C, Base::U];
        let seq = Sequence::from(seq_bases);
        let prepared = SuffixIndex::build_for_seed(&seq, 1).unwrap();
        let suffixes = prepared.view();
        let parts =
            ParallelSaTraverser::<fn(SeedMatch)>::partition(suffixes, &(0..suffixes.len()), 0);

        assert_eq!(parts[1] - parts[0], 1); // A
        assert_eq!(parts[2] - parts[1], 1); // C
        assert_eq!(parts[3] - parts[2], 1); // G
        assert_eq!(parts[4] - parts[3], 0); // N
        assert_eq!(parts[5] - parts[4], 1); // U
    }

    #[test]
    fn lower_bound_matches_a_linear_scan_above_the_cutoff() {
        let bases = [Base::A, Base::C, Base::G, Base::U];
        let seq: Vec<Base> = (0..LINEAR_PARTITION_CUTOFF * 2)
            .map(|i| bases[i % bases.len()])
            .collect();
        let prepared = SuffixIndex::build_for_seed(&seq, 1).unwrap();
        let suffixes = prepared.view();
        assert!(suffixes.len() > LINEAR_PARTITION_CUTOFF);

        for depth in [0, 1] {
            for rank in 1..BASE_COUNT as u8 {
                let expected = (0..suffixes.len())
                    .find(|&i| unsafe { suffixes.base_unchecked(i, depth).as_u8() } >= rank)
                    .unwrap_or(suffixes.len());
                assert_eq!(
                    ParallelSaTraverser::<fn(SeedMatch)>::lower_bound(
                        suffixes,
                        0..suffixes.len(),
                        depth,
                        rank,
                    ),
                    expected,
                    "rank {rank} at depth {depth}"
                );
            }
        }
    }

    #[test]
    fn terminal_gap_stops_descent_at_sequence_end() {
        let len = 40;
        let query = vec![Base::A; len];
        let target = vec![Base::U; len];
        let prepared_query = SuffixIndex::build_for_seed(&query, 1).unwrap();
        let prepared_target = SuffixIndex::build_for_seed(&target, 1).unwrap();
        let query_suffixes = prepared_query.view();
        let target_suffixes = prepared_target.view();
        let mut emitted_lengths = Vec::new();

        // Every suffix is retained, so each depth partitions a multi-suffix
        // interval in which one suffix has just reached the terminal Gap.
        assert_eq!(query_suffixes.sequence().len(), len + 1);
        assert_eq!(query_suffixes.len(), len);
        traverse::<false, _>(
            query_suffixes,
            target_suffixes,
            len,
            len + 1,
            0,
            1,
            0,
            &mut |seed| emitted_lengths.push(seed.seed_len),
        );

        assert_eq!(emitted_lengths, [len]);
    }

    #[test]
    fn singleton_handoff_does_not_double_emit() {
        // Distinct first symbols force singleton intervals at depth=1 in both
        // query and target branches, exercising the singleton fast-path handoff.
        let q_seq = Sequence::from(vec![Base::A, Base::C]);
        let t_seq = Sequence::from(vec![Base::U, Base::G]);
        let prepared_query = SuffixIndex::build_for_seed(&q_seq, 1).unwrap();
        let prepared_target = SuffixIndex::build_for_seed(&t_seq, 1).unwrap();
        let query = prepared_query.view();
        let target = prepared_target.view();

        let cfg = SeedConfig {
            seed_start: None,
            seed_end: None,
            seed_length: Some(1),
            seed_wobble: false,
            no_max_prune: false,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        };
        let mut seen = std::collections::HashSet::new();
        let mut ctx = ParallelSaTraverser {
            query,
            target,
            min_len: 1,
            max_len: 2,
            max_mm: cfg.max_mismatches,
            min_prefix: cfg.min_prefix_matches,
            min_suffix: cfg.min_suffix_matches,
            on_match: &mut |m| {
                assert!(
                    seen.insert((
                        m.query_interval.start,
                        m.query_interval.end,
                        m.target_interval.start,
                        m.target_interval.end,
                        m.seed_len
                    )),
                    "duplicate seed match emitted: {:?}",
                    m
                );
            },
        };

        let q_range = 0..ctx.query.len();
        let t_range = 0..ctx.target.len();
        ctx.recurse::<false>(&q_range, &t_range, 0, 0, 0);
    }
}
