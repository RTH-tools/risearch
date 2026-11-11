//! Search results must be invariant under interface-level settings.
//!
//! Each check runs the production search twice with one setting changed and
//! requires the same named hits.  Search completion order is deliberately
//! unspecified, so a sorted named hit key is the stable contract exercised here.
//! These compare two production runs against each other.
//!
//! The corpus is the shared `tests/data` set: 19 transcripts of one gene, so
//! a query lands in several records at different global offsets, and the same
//! local sequence recurs across records.  Toy records cannot pose that question.

use super::test_support::{
    agreement_config, collision_config, no_dedup_config, search, QUERIES, TARGETS, TARGET_RECORDS,
};
use super::SearchHit;
use crate::types::Strand;
use crate::{SearchConfig, Sequence, TargetRegistry};
use rstest::rstest;
use std::collections::HashSet;
use std::num::NonZeroUsize;

fn alignment_key(hit: &SearchHit) -> Option<String> {
    hit.alignment.as_ref().map(|columns| {
        columns
            .iter()
            .map(|column| {
                format!(
                    "{}{}{}",
                    column.class().symbol(),
                    column.query().to_byte() as char,
                    column.target().to_byte() as char
                )
            })
            .collect()
    })
}

fn hit_key(
    hit: &SearchHit,
    queries: &[(String, Sequence)],
    targets: &TargetRegistry,
) -> impl Ord + std::fmt::Debug {
    (
        queries[hit.query_idx as usize].0.clone(),
        targets.get_name(hit.target_index()).to_string(),
        (hit.q_start, hit.q_end),
        (hit.t_start, hit.t_end),
        hit.strand,
        hit.energy,
        alignment_key(hit),
    )
}

fn hit_keys(
    hits: Vec<SearchHit>,
    queries: &[(String, Sequence)],
    targets: &TargetRegistry,
) -> Vec<impl Ord + std::fmt::Debug> {
    let mut keys: Vec<_> = hits
        .iter()
        .map(|hit| hit_key(hit, queries, targets))
        .collect();
    keys.sort();
    keys
}

/// Records share one concatenated text, so each hit's name and coordinates are
/// read through a per-record offset.  Regrouping the records moves every offset
/// and renumbers every target index; the named hits must not move with them.
/// 19 is prime, so each group size also leaves a ragged final group.
#[rstest]
fn record_grouping_preserves_named_hits(
    #[values(
        SearchConfig::default,
        agreement_config,
        collision_config,
        no_dedup_config
    )]
    preset: fn() -> SearchConfig,
    #[values(1, 4, 7)] group_size: usize,
) {
    let config = preset();
    let combined = TargetRegistry::build(TARGET_RECORDS.clone(), Some(NonZeroUsize::MIN)).unwrap();

    let hits = search(&QUERIES, &combined, &config).unwrap();
    let touched: HashSet<_> = hits.iter().map(|hit| hit.target_index()).collect();
    assert!(
        touched.len() >= 2,
        "corpus must place hits in more than one record, or grouping is vacuous"
    );
    let combined_keys = hit_keys(hits, &QUERIES, &combined);

    let mut grouped = Vec::new();
    for group in TARGET_RECORDS.chunks(group_size) {
        let part = TargetRegistry::build(group.to_vec(), Some(NonZeroUsize::MIN)).unwrap();
        grouped.extend(hit_keys(
            search(&QUERIES, &part, &config).unwrap(),
            &QUERIES,
            &part,
        ));
    }
    grouped.sort();

    assert_eq!(combined_keys, grouped, "record grouping changed named hits");
}

/// Queries are seeded, extended and deduplicated one at a time on a rayon
/// worker, so splitting the query file must not change any query's own hits.
#[rstest]
fn query_partitioning_preserves_named_hits(
    #[values(
        SearchConfig::default,
        agreement_config,
        collision_config,
        no_dedup_config
    )]
    preset: fn() -> SearchConfig,
) {
    let config = preset();
    let combined_keys = hit_keys(
        search(&QUERIES, &TARGETS, &config).unwrap(),
        &QUERIES,
        &TARGETS,
    );
    assert!(!combined_keys.is_empty(), "corpus must produce hits");

    let mut split = Vec::new();
    for one in QUERIES.chunks(1) {
        split.extend(hit_keys(
            search(one, &TARGETS, &config).unwrap(),
            one,
            &TARGETS,
        ));
    }
    split.sort();

    assert_eq!(
        combined_keys, split,
        "query partitioning changed named hits"
    );
}

/// `delta_g` is a post-extension accept/reject, so tightening it must remove
/// exactly the hits above the new bound and nothing else.  The cutoff is an
/// energy the loose run actually reported, at three points of its distribution,
/// so no quantization stands between the two runs.
#[rstest]
fn stricter_cutoff_returns_exactly_the_eligible_hits(
    #[values(
        SearchConfig::default,
        agreement_config,
        collision_config,
        no_dedup_config
    )]
    preset: fn() -> SearchConfig,
    #[values(0.1, 0.5, 0.9)] quantile: f64,
) {
    let config = preset();

    let loose_hits = search(&QUERIES, &TARGETS, &config).unwrap();
    let mut energies: Vec<_> = loose_hits.iter().map(|hit| hit.energy).collect();
    energies.sort();
    energies.dedup();
    assert!(
        energies.len() >= 2,
        "corpus must contain both accepted and rejected energies"
    );

    let mut strict = config.clone();
    strict.filter.delta_g = energies[(quantile * (energies.len() - 1) as f64) as usize];

    let loose_count = loose_hits.len();
    let expected = hit_keys(
        loose_hits
            .into_iter()
            .filter(|hit| hit.energy <= strict.filter.delta_g)
            .collect(),
        &QUERIES,
        &TARGETS,
    );
    assert!(
        !expected.is_empty() && expected.len() < loose_count,
        "strict cutoff should retain and reject corpus hits"
    );

    let strict_keys = hit_keys(
        search(&QUERIES, &TARGETS, &strict).unwrap(),
        &QUERIES,
        &TARGETS,
    );
    assert_eq!(
        strict_keys, expected,
        "cutoff lost eligible hits or retained ineligible hits"
    );
}

/// Queries are spread over rayon workers and reach the sink in completion
/// order, so the worker count must not change which hits are found.
#[rstest]
fn worker_count_preserves_named_hits(
    #[values(
        SearchConfig::default,
        agreement_config,
        collision_config,
        no_dedup_config
    )]
    preset: fn() -> SearchConfig,
) {
    let config = preset();
    let keys_with = |workers: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        let hits = pool.install(|| search(&QUERIES, &TARGETS, &config).unwrap());
        hit_keys(hits, &QUERIES, &TARGETS)
    };

    let single = keys_with(1);
    assert!(!single.is_empty(), "corpus must produce hits");
    assert_eq!(single, keys_with(4), "worker count changed named hits");
}

/// A target and its reverse complement hold the same two physical strands,
/// swapped. Every hit must therefore reappear on the other strand, at FASTA
/// coordinates mirrored within its record.
#[rstest]
fn reverse_complemented_targets_mirror_named_hits(
    #[values(
        SearchConfig::default,
        agreement_config,
        collision_config,
        no_dedup_config
    )]
    preset: fn() -> SearchConfig,
) {
    let config = preset();
    let flipped_records = TARGET_RECORDS
        .iter()
        .map(|(name, sequence)| {
            let flipped: Vec<_> = sequence
                .iter()
                .rev()
                .map(|base| base.complement())
                .collect();
            (name.clone(), Sequence::from(flipped))
        })
        .collect();
    let flipped = TargetRegistry::build(flipped_records, Some(NonZeroUsize::MIN)).unwrap();

    let mirrored = search(&QUERIES, &flipped, &config)
        .unwrap()
        .into_iter()
        .map(|mut hit| {
            let len = TARGET_RECORDS[hit.target_index()].1.len();
            (hit.t_start, hit.t_end) = (len - 1 - hit.t_end, len - 1 - hit.t_start);
            hit.strand = match hit.strand {
                Strand::Forward => Strand::Reverse,
                Strand::Reverse => Strand::Forward,
            };
            hit
        })
        .collect();

    let expected = hit_keys(
        search(&QUERIES, &TARGETS, &config).unwrap(),
        &QUERIES,
        &TARGETS,
    );
    assert!(!expected.is_empty(), "corpus must produce hits");
    assert_eq!(
        hit_keys(mirrored, &QUERIES, &flipped),
        expected,
        "reverse complement did not mirror named hits"
    );
}
