//! Queries prepared for search.
//!
//! [`QueryRegistry`] is built once from normalized records against a
//! [`SeedConfig`], so each [`Query`] arrives with its seed interval and length
//! bounds already resolved.

use crate::error::{Error, Result};
use rayon::prelude::*;
use std::collections::HashSet;
use std::ops::{Index, Range};

use crate::config::SeedConfig;
use crate::seq::Sequence;
use crate::types::Base;

/// Query prepared once at load time for the current seed configuration.
pub struct Query {
    /// Sequence identifier
    name: String,
    /// Full forward sequence used for extension and output
    sequence: Sequence,
    /// Pre-computed normalized seed interval bounds on the full query
    pub(crate) seed_interval: Range<usize>,
    /// Minimum admissible seed length for this query after normalization.
    pub(crate) min_seed_len: usize,
}

impl Query {
    fn from_parts(id: String, sequence: Sequence, config: &SeedConfig) -> Result<Self> {
        let q_len = sequence.len();

        // Compute seed interval once and fail early at boundary if invalid.
        let (start1, end1, min_seed_len) = config
            .resolve(q_len)
            .map_err(|err| Error::Config(format!("Invalid seed spec for query '{id}': {err}")))?;
        let seed_interval = (start1 - 1)..end1;

        Ok(Self {
            name: id,
            sequence,
            seed_interval,
            min_seed_len,
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
    /// Normalized seed-interval bounds on the full sequence.
    pub fn seed_interval(&self) -> Range<usize> {
        self.seed_interval.clone()
    }
}

/// Registry of queries. Each `Query` carries its own metadata (sequence, seed
/// interval); the per-query suffix array used for seeding is built
/// on demand inside the seeding worker (see `seed::engine::SeedingEngine::seed_query`),
/// so the registry holds no combined query SA.
pub struct QueryRegistry {
    inner: Vec<Query>,
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
        self.inner.iter().enumerate()
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
    /// Record order is kept and duplicate IDs are rejected. Seed intervals are
    /// resolved in parallel via rayon.
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

        Ok(Self { inner: entries })
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
        let empty = QueryRegistry { inner: Vec::new() };
        assert!(empty.is_empty());

        let query = make_query_data(
            Sequence::from(vec![Base::A, Base::C, Base::G]),
            None,
            None,
            None,
        );
        let loaded = QueryRegistry { inner: vec![query] };
        assert!(!loaded.is_empty());
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].name(), "q");
    }

    #[test]
    fn explicit_interval_resolves_to_half_open_bounds() {
        let sequence = Sequence::from(vec![Base::A, Base::U, Base::G, Base::C, Base::A, Base::U]);
        let query = make_query_data(sequence, Some(2), Some(5), Some(2));

        assert_eq!(query.seed_interval, 1..5);
        assert_eq!(query.min_seed_len, 2);
    }

    #[test]
    fn default_interval_spans_the_whole_query() {
        let sequence = Sequence::from(vec![Base::A, Base::G, Base::C, Base::U, Base::A]);
        let query = make_query_data(sequence, None, None, Some(3));

        assert_eq!(query.seed_interval, 0..5);
        assert_eq!(query.min_seed_len, 3);
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
}
