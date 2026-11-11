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

use super::test_support::{config, search, QUERIES, TARGETS, TARGET_RECORDS};
use super::SearchHit;
use crate::{Energy, SearchConfig, Sequence, TargetRegistry};
use std::collections::HashSet;
use std::num::NonZeroUsize;

/// Mirrors the search settings of `tests/cli_library_agreement.rs`. Alignment
/// columns, on by default, let the hit key separate duplexes sharing a
/// bounding box and an energy.
fn agreement_config() -> SearchConfig {
    config(|c| {
        c.seed.seed_length = Some(6);
        c.seed.max_mismatches = 2;
        c.seed.min_prefix_matches = 2;
        c.seed.min_suffix_matches = 2;
        c.extend.max_extension = 30;
        c.filter.delta_g = Energy::from_kcal(-5.0);
    })
}

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

fn hit_key(hit: &SearchHit, queries: &[(String, Sequence)], targets: &TargetRegistry) -> String {
    // Energy is integer-scaled internally.  Re-quantizing the public kcal/mol
    // value keeps this helper independent of the private representation while
    // retaining exact equality at the advertised precision.
    let energy = (hit.energy.to_kcal() * 10_000.0).round() as i64;
    format!(
        "{}\0{}\0{}:{}\0{}:{}\0{}\0{}\0{}",
        queries[hit.query_idx as usize].0,
        targets.get_name(hit.target_index()),
        hit.q_start,
        hit.q_end,
        hit.t_start,
        hit.t_end,
        hit.strand,
        energy,
        alignment_key(hit).unwrap_or_default()
    )
}

fn hit_keys(
    hits: Vec<SearchHit>,
    queries: &[(String, Sequence)],
    targets: &TargetRegistry,
) -> Vec<String> {
    let mut keys: Vec<_> = hits
        .into_iter()
        .map(|hit| hit_key(&hit, queries, targets))
        .collect();
    keys.sort();
    keys
}

/// Records share one concatenated text, so each hit's name and coordinates are
/// read through a per-record offset.  Regrouping the records moves every offset
/// and renumbers every target index; the named hits must not move with them.
/// 19 is prime, so each group size also leaves a ragged final group.
#[rstest::rstest]
fn record_grouping_preserves_named_hits(#[values(1, 4, 7)] group_size: usize) {
    let config = agreement_config();
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
#[test]
fn query_partitioning_preserves_named_hits() {
    let config = agreement_config();
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
#[rstest::rstest]
fn stricter_cutoff_returns_exactly_the_eligible_hits(#[values(0.1, 0.5, 0.9)] quantile: f64) {
    let config = agreement_config();

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
