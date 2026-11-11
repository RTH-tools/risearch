mod engine;
mod parallel_sa;

#[cfg(test)]
pub(crate) mod reference_tests;

pub use engine::SeedingEngine;

use std::ops::Range;

use crate::types::{Base, Strand};

/// A candidate seed match found during suffix array search.
///
/// Query coordinates address the full query sequence. Target coordinates
/// address the strand-selected physical duplex view returned by
/// `TargetRegistry::target(target_idx, strand)`. Both ranges are half-open and
/// always have the same non-zero length.
#[derive(Debug, Clone)]
pub struct SeedHit {
    query_start: usize,
    target_idx: usize,
    target_start: usize,
    len: usize,
    strand: Strand,
}

impl SeedHit {
    pub(crate) fn new(
        query_start: usize,
        target_idx: usize,
        target_start: usize,
        len: usize,
        strand: Strand,
    ) -> Self {
        debug_assert!(len > 0, "seed matches must be non-empty");
        Self {
            query_start,
            target_idx,
            target_start,
            len,
            strand,
        }
    }

    /// Half-open range in the full query sequence.
    #[inline]
    pub(crate) fn query_range(&self) -> Range<usize> {
        self.query_start..self.query_start + self.len
    }

    /// Index of the target in its registry.
    #[inline]
    pub(crate) fn target_idx(&self) -> usize {
        self.target_idx
    }

    /// Half-open range in the strand-selected physical target view.
    #[inline]
    pub(crate) fn target_range(&self) -> Range<usize> {
        self.target_start..self.target_start + self.len
    }

    /// Shared length of the query and target ranges.
    #[cfg(test)]
    #[inline]
    pub(crate) fn seed_len(&self) -> usize {
        self.len
    }

    /// Strand used to select the physical target view containing the match.
    #[inline]
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Whether the match cannot be grown by one more pairing column on either
    /// side. A non-maximal seed is a shorter copy of a longer one, so extending
    /// it only rediscovers the same duplex.
    ///
    /// `query` and `target` are the same views seeding used; `seed_interval`
    /// bounds the query region seeds were drawn from.
    pub(crate) fn is_maximal(
        &self,
        query: &[Base],
        target: &[Base],
        seed_interval: &Range<usize>,
        seed_wobble: bool,
    ) -> bool {
        let pairs = |qi: usize, ti: usize| query[qi].pair_type(target[ti]).is_match(seed_wobble);
        let query_range = self.query_range();
        let target_range = self.target_range();

        if query_range.start > seed_interval.start
            && target_range.start > 0
            && pairs(query_range.start - 1, target_range.start - 1)
        {
            return false;
        }

        if query_range.end < seed_interval.end
            && target_range.end < target.len()
            && pairs(query_range.end, target_range.end)
        {
            return false;
        }

        true
    }
}
