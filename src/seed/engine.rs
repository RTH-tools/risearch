use crate::config::SeedConfig;
use crate::error::{Error, Result};
use crate::index::sa::{SuffixIndex, SuffixIndexView};
use crate::index::store::TargetRegistry;
use crate::index::TargetView;
use crate::registry::{Query, QueryRegistry};

use super::parallel_sa::{traverse, SeedMatch};
use super::SeedHit;

pub struct SeedingEngine<'a> {
    queries: &'a QueryRegistry,
    tview: TargetView<'a>,
}

impl<'a> SeedingEngine<'a> {
    pub fn new(queries: &'a QueryRegistry, targets: &'a TargetRegistry) -> Self {
        Self {
            queries,
            tview: targets.view(),
        }
    }

    /// Collect every seed grouped by query, materializing them all in memory.
    /// Retained for benchmarks and unit tests; the search paths stream one query
    /// at a time via [`seed_query`](Self::seed_query).
    pub fn run(&self, config: &SeedConfig) -> Result<Vec<(usize, Vec<SeedHit>)>> {
        let mut groups = Vec::new();
        for qi in 0..self.queries.len() {
            let mut seeds = Vec::new();
            self.seed_query(qi, config, |seed| seeds.push(seed))?;
            if !seeds.is_empty() {
                groups.push((qi, seeds));
            }
        }
        Ok(groups)
    }

    /// Seed a single query against the shared target SA: build the query's own
    /// (tiny) suffix array, traverse, emit each hit to `on_seed`, then drop it.
    /// Every seed for query `qi` is produced here, so a caller can own one query
    /// end-to-end (per-query parallelism + local dedup).
    pub fn seed_query<F: FnMut(SeedHit)>(
        &self,
        qi: usize,
        config: &SeedConfig,
        mut on_seed: F,
    ) -> Result<()> {
        // `traverse` is monomorphized on whether G-U wobble pairs are allowed.
        let query = &self.queries[qi];
        if config.seed_wobble {
            self.seed_one::<true, _>(qi, query, config, &mut on_seed)
        } else {
            self.seed_one::<false, _>(qi, query, config, &mut on_seed)
        }
    }

    fn seed_one<const WOBBLE: bool, F: FnMut(SeedHit)>(
        &self,
        qi: usize,
        query: &Query,
        config: &SeedConfig,
        on_seed: &mut F,
    ) -> Result<()> {
        // Returning here also skips the per-query suffix array build.
        if query.min_seed_len > query.max_seed_len {
            return Ok(());
        }
        let tview = self.tview;

        let seed = query.seed_sequence();
        let prepared_query =
            SuffixIndex::build_for_seed(seed, query.min_seed_len).map_err(|err| {
                Error::Index(format!(
                    "building suffix array for query '{}': {err}",
                    query.name()
                ))
            })?;
        let query_suffixes = prepared_query.view();

        traverse::<WOBBLE, _>(
            query_suffixes,
            tview.suffixes(),
            query.min_seed_len,
            query.max_seed_len,
            config.max_mismatches,
            config.min_prefix_matches,
            config.min_suffix_matches,
            &mut |m| {
                emit_seed_match::<WOBBLE, _>(
                    qi,
                    query,
                    query_suffixes,
                    tview,
                    config.no_max_prune,
                    m,
                    &mut *on_seed,
                )
            },
        );
        Ok(())
    }
}

/// Emit every (query position × target position) seed for one `SeedMatch` of a
/// single query. Query positions come straight out of the single-query SA and
/// map back through [`Query::map_seed_pos`]; the target view maps its global
/// suffix positions through the index directory.
///
/// Non-maximal seeds are dropped here unless `no_max_prune`: they are shorter
/// copies of a longer match, so extending them only rediscovers the same duplex.
fn emit_seed_match<const WOBBLE: bool, F: FnMut(SeedHit)>(
    qi: usize,
    query: &Query,
    query_suffixes: SuffixIndexView<'_>,
    tview: TargetView<'_>,
    no_max_prune: bool,
    raw_match: SeedMatch,
    on_seed: &mut F,
) {
    let seed_len = raw_match.seed_len;
    let query_bases = query.sequence();
    let seed_interval = query.seed_interval();
    for &query_sa_pos in query_suffixes.suffix_positions(raw_match.query_interval.clone()) {
        let Some(query_start) = query.map_seed_pos(query_sa_pos as usize, seed_len) else {
            continue;
        };

        for &target_sa_pos in tview.suffix_positions(raw_match.target_interval.clone()) {
            let Some((target_idx, strand, target_start)) =
                tview.map_seed_pos(target_sa_pos as usize, seed_len)
            else {
                continue;
            };

            let hit = SeedHit::new(qi, query_start, target_idx, target_start, seed_len, strand);
            if no_max_prune
                || hit.is_maximal(
                    query_bases,
                    tview.target(target_idx, strand),
                    &seed_interval,
                    WOBBLE,
                )
            {
                on_seed(hit);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use crate::config::SeedConfig;
    use crate::fastx::read_sequences;
    use crate::index::store::TargetRegistry;
    use crate::registry::QueryRegistry;
    use crate::types::Strand;

    use super::*;

    fn build_store(fasta: &str) -> (TargetRegistry, tempfile::TempDir) {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        fs_err::write(&fasta_path, fasta).unwrap();
        let targets = read_sequences(&fasta_path).unwrap();
        (TargetRegistry::build(targets, None).unwrap(), dir)
    }

    fn build_queries(fasta: &str, config: &SeedConfig) -> QueryRegistry {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(fasta.as_bytes()).unwrap();
        QueryRegistry::build(read_sequences(file.path()).unwrap(), config).unwrap()
    }

    #[test]
    fn collect_maps_seed_interval_back_to_full_query_coordinates() {
        let config = SeedConfig {
            seed_start: Some(3),
            seed_end: Some(4),
            seed_length: Some(2),
            seed_wobble: false,
            no_max_prune: true,
            ..Default::default()
        };
        let queries = build_queries(">q1\nGGAC\n", &config);
        let (targets, _dir) = build_store(">t1\nGU\n");

        let groups = SeedingEngine::new(&queries, &targets).run(&config).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0, 0);
        assert_eq!(groups[0].1.len(), 1);

        let seed = &groups[0].1[0];
        assert_eq!(seed.query_idx(), 0);
        assert_eq!(seed.query_range(), 2..4);
        assert_eq!(seed.target_idx(), 0);
        assert_eq!(seed.target_range(), 0..2);
        assert_eq!(seed.seed_len(), 2);
        assert_eq!(seed.strand(), Strand::Forward);
    }

    #[test]
    fn collect_preserves_same_block_coordinates_on_both_strands() {
        let config = SeedConfig {
            seed_length: Some(2),
            seed_wobble: false,
            no_max_prune: true,
            ..Default::default()
        };
        let queries = build_queries(">q1\nCG\n", &config);
        let (targets, _dir) = build_store(">t1\nAACGU\n");

        let groups = SeedingEngine::new(&queries, &targets).run(&config).unwrap();
        assert_eq!(groups.len(), 1);
        let seeds = &groups[0].1;
        assert_eq!(seeds.len(), 2);
        assert!(seeds
            .iter()
            .any(|seed| seed.strand() == Strand::Forward && seed.target_range() == (1..3)));
        assert!(seeds
            .iter()
            .any(|seed| seed.strand() == Strand::Reverse && seed.target_range() == (2..4)));
    }
}
