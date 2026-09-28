//! Loaded queries and the name/index bookkeeping shared with the target store.
//!
//! [`QueryRegistry`] is built once from FASTA against a [`SeedConfig`], so each
//! [`Query`] arrives with its seed interval, length bounds, and N-prefix already
//! resolved.

use crate::error::{Error, Result};
use rayon::prelude::*;
use std::collections::HashSet;
use std::ops::{Index, Range};

use crate::config::SeedConfig;
use crate::seq::Sequence;
use crate::types::Base;

#[doc(hidden)]
pub struct Registry<T> {
    entries: Vec<T>,
}

impl<T> Registry<T> {
    pub fn new(entries: Vec<T>) -> Self {
        Self { entries }
    }

    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (usize, &T)> {
        self.entries.iter().enumerate()
    }
}

impl<T> Index<usize> for Registry<T> {
    type Output = T;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.entries[idx]
    }
}

/// Query prepared once at load time for the current seed configuration.
pub struct Query {
    /// Sequence identifier
    name: String,
    /// Full forward sequence used for extension and output
    sequence: Sequence,
    /// Seed-interval slice; the per-query seeding suffix array is built from this
    seed_sequence: Sequence,
    /// Pre-computed normalized seed interval bounds on the full query
    pub(crate) seed_interval: Range<usize>,
    /// Minimum admissible seed length for this query after normalization.
    pub(crate) min_seed_len: usize,
    /// Maximum admissible seed length for this query after normalization.
    pub(crate) max_seed_len: usize,
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: Vec<u32>,
    /// Fast path when query has no Ns
    has_n_any: bool,
}

impl Query {
    fn from_parts(id: String, sequence: Sequence, config: &SeedConfig) -> Result<Self> {
        let q_len = sequence.len();

        // Compute N-prefix for O(1) N-checking
        let mut n_prefix = Vec::with_capacity(q_len + 1);
        n_prefix.push(0);
        let mut n_total = 0;
        for &base in sequence.iter() {
            if base == Base::N {
                n_total += 1;
            }
            n_prefix.push(n_total);
        }
        let has_n_any = n_total != 0;

        // Compute seed interval once and fail early at boundary if invalid.
        let (start1, end1, min_seed_len) = config
            .resolve(q_len)
            .map_err(|err| Error::Config(format!("Invalid seed spec for query '{id}': {err}")))?;
        let seed_interval = (start1 - 1)..end1;
        let max_seed_len = seed_interval.end.saturating_sub(seed_interval.start);
        let seed_sequence =
            Sequence::from(sequence[seed_interval.start..seed_interval.end].to_vec());

        Ok(Self {
            name: id,
            sequence,
            seed_sequence,
            seed_interval,
            min_seed_len,
            max_seed_len,
            n_prefix,
            has_n_any,
        })
    }

    #[inline]
    /// Sequence identifier from the FASTA header.
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline(always)]
    /// Full forward sequence, as used for extension and output.
    pub fn sequence(&self) -> &[Base] {
        &self.sequence
    }

    #[inline]
    /// Prefix sums of `N` counts, for O(1) N-free span checks.
    pub fn n_prefix(&self) -> &[u32] {
        &self.n_prefix
    }

    #[inline]
    /// Whether the sequence contains any `N`.
    pub fn has_n_any(&self) -> bool {
        self.has_n_any
    }

    #[inline]
    /// Normalized seed-interval bounds on the full sequence.
    pub fn seed_interval(&self) -> Range<usize> {
        self.seed_interval.clone()
    }

    #[inline(always)]
    /// The seed interval as its own slice; the per-query seeding SA is built from this.
    pub fn seed_sequence(&self) -> &[Base] {
        &self.seed_sequence
    }

    /// Map a local position within the seed view to a global query start coordinate.
    ///
    /// Returns `None` if the seed overflows the query's seed-search interval,
    /// violates the query-specific length constraints, or contains an 'N' base.
    #[inline]
    pub fn map_seed_pos(&self, local_pos: usize, seed_len: usize) -> Option<usize> {
        if local_pos + seed_len > self.seed_sequence.len() {
            return None;
        }

        if seed_len < self.min_seed_len || seed_len > self.max_seed_len {
            return None;
        }

        let query_start = self.seed_interval.start + local_pos;

        if self.has_n_any && self.n_prefix[query_start + seed_len] != self.n_prefix[query_start] {
            return None;
        }

        Some(query_start)
    }
}

/// Registry of queries. Each `Query` carries its own metadata (sequence, seed
/// interval, N-prefix); the per-query suffix array used for seeding is built
/// on demand inside the seeding worker (see `seed::engine::SeedingEngine::seed_query`),
/// so the registry holds no combined query SA.
pub struct QueryRegistry {
    inner: Registry<Query>,
}

impl QueryRegistry {
    /// Number of loaded queries.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Whether no queries were loaded.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate `(index, query)` pairs in load order.
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Query)> {
        self.inner.iter()
    }
}

impl Index<usize> for QueryRegistry {
    type Output = Query;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.inner[idx]
    }
}

impl QueryRegistry {
    /// Prepare queries from normalized records against `config`; this is what
    /// [`run_search`](crate::run_search) does with its `queries`.
    ///
    /// Record order is kept and duplicate IDs are rejected. Query SA
    /// construction is parallelised via rayon.
    pub fn build(queries: Vec<(String, Sequence)>, config: &SeedConfig) -> Result<Self> {
        if queries.is_empty() {
            return Err(Error::Input("No query sequences provided".into()));
        }

        let mut seen = HashSet::with_capacity(queries.len());
        for (id, _) in &queries {
            if !seen.insert(id.as_str()) {
                return Err(Error::Input(format!(
                    "Duplicate query id '{id}' across inputs"
                )));
            }
        }

        let entries = queries
            .into_par_iter()
            .map(|(id, sequence)| Query::from_parts(id, sequence, config))
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            inner: Registry::new(entries),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastx::read_sequences_from;

    fn make_query_data(
        sequence: Sequence,
        seed_start: Option<i64>,
        seed_end: Option<i64>,
        seed_length: Option<i64>,
    ) -> Query {
        let cfg = SeedConfig {
            seed_start,
            seed_end,
            seed_length,
            seed_wobble: false,
            no_max_prune: false,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        };
        Query::from_parts("q".into(), sequence, &cfg).expect("query data")
    }

    #[test]
    fn query_registry_exposes_the_loaded_queries() {
        let empty = QueryRegistry {
            inner: Registry::new(Vec::new()),
        };
        assert!(empty.is_empty());

        let query = make_query_data(
            Sequence::from(vec![Base::A, Base::C, Base::G]),
            None,
            None,
            None,
        );
        let loaded = QueryRegistry {
            inner: Registry::new(vec![query]),
        };
        assert!(!loaded.is_empty());
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].name(), "q");
    }

    #[test]
    fn has_n_any_is_false_without_an_n() {
        let clean = make_query_data(
            Sequence::from(vec![Base::A, Base::C, Base::G]),
            None,
            None,
            None,
        );
        assert!(!clean.has_n_any());
    }

    #[test]
    fn seed_sequence_matches_interval_slice() {
        let sequence = Sequence::from(vec![Base::A, Base::U, Base::G, Base::C, Base::A, Base::U]);
        let query = make_query_data(sequence, Some(2), Some(5), Some(2));

        assert_eq!(query.seed_interval, 1..5);
        assert_eq!(query.min_seed_len, 2);
        assert_eq!(query.max_seed_len, 4);
        assert_eq!(query.seed_sequence(), &query.sequence()[1..5]);
    }

    #[test]
    fn seed_sequence_tracks_only_valid_interval_starts() {
        let sequence = Sequence::from(vec![Base::A, Base::G, Base::C, Base::U, Base::A]);
        let query = make_query_data(sequence, None, None, Some(3));

        assert_eq!(query.seed_interval, 0..5);
        assert_eq!(query.min_seed_len, 3);
        assert_eq!(query.max_seed_len, 5);
        assert_eq!(query.seed_sequence().len(), 5);
        assert_eq!(query.seed_sequence(), query.sequence());
    }

    fn seed_config() -> SeedConfig {
        SeedConfig {
            seed_start: None,
            seed_end: None,
            seed_length: Some(4),
            seed_wobble: true,
            no_max_prune: false,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        }
    }

    fn records(content: &str) -> Vec<(String, Sequence)> {
        read_sequences_from(content.as_bytes(), "inline").unwrap()
    }

    #[test]
    fn build_preserves_record_order() {
        let queries = records(">alpha\nACGUACGU\n>beta\nUUUUAAAA\n>gamma\nGGGGCCCC\n");

        let registry = QueryRegistry::build(queries, &seed_config()).unwrap();

        assert_eq!(registry.len(), 3);
        assert_eq!(registry[0].name(), "alpha");
        assert_eq!(registry[1].name(), "beta");
        assert_eq!(registry[2].name(), "gamma");
    }

    #[test]
    fn build_rejects_duplicate_ids_across_inputs() {
        let mut queries = records(">seq1\nACGUACGU\n");
        queries.extend(records(">seq1\nUUUUAAAA\n"));

        let msg = QueryRegistry::build(queries, &seed_config())
            .err()
            .unwrap()
            .to_string();
        assert!(
            msg.contains("Duplicate"),
            "expected duplicate error, got: {msg}"
        );
    }

    #[test]
    fn build_rejects_no_queries() {
        assert!(QueryRegistry::build(Vec::new(), &seed_config()).is_err());
    }

    #[test]
    fn seeds_spanning_an_n_are_rejected() {
        let sequence = Sequence::from(vec![Base::A, Base::N, Base::G, Base::C]);
        let query = make_query_data(sequence, None, None, Some(2));

        assert!(query.has_n_any());
        assert_eq!(query.n_prefix(), &[0, 0, 1, 1, 1]);
        assert_eq!(query.map_seed_pos(0, 2), None);
        assert_eq!(query.map_seed_pos(2, 2), Some(2));
        assert_eq!(query.map_seed_pos(0, 1), None);
    }
}
