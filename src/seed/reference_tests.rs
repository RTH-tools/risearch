//! Independent seed eligibility checks over normalized domain values.
//!
//! `enumerate_seed_spans` is the reference algorithm; it never traverses a suffix
//! array. `assert_seeds_match_enumeration` builds production registries and compares
//! their seeds with that reference. `Base`, `Sequence`, and `Strand` remain the
//! shared representations; seed eligibility and strand arrangement are checked
//! explicitly here.

use std::ops::Range;

use crate::fastx::read_sequences;
use crate::registry::QueryRegistry;
use crate::{Base, SeedConfig, Sequence, Strand, TargetRegistry};

/// Arrange a normalized FASTA target in physical duplex-column order.
///
/// Alongside the query's 5' -> 3' order, the forward strand uses reverse(target)
/// and the reverse strand uses complement(target). Constructing those views here
/// checks the target registry's orientation independently; base complementation
/// itself belongs to `Base` and is shared.
pub(crate) fn target_in_duplex_order(target: &[Base], strand: Strand) -> Vec<Base> {
    match strand {
        Strand::Forward => target.iter().rev().copied().collect(),
        Strand::Reverse => target.iter().map(|base| base.complement()).collect(),
    }
}

/// Brute-force seeding oracle: test every pair of starts and every fitting length.
///
/// `query` is in 5' -> 3' order; `target` is in duplex-column order (see
/// [`target_in_duplex_order`]). Results are `(query_start, target_start, length)`,
/// with zero-based offsets into these slices, not original FASTA coordinates.
/// Results follow query-start, target-start, then length order; no deduplication
/// is applied, so callers can compare collections with multiplicity.
///
/// Callers supply a nonempty, half-open `query_interval` within `query` and a
/// positive `minimum_length`. These resolved values take precedence over the
/// corresponding config fields, keeping config resolution out of the oracle.
/// An interval or target shorter than the minimum produces no seeds.
///
/// A candidate must contain no ambiguous bases and stay within the mismatch
/// budget. Only candidates with mismatches must satisfy the matching prefix and
/// suffix lengths and contain no perfect run of at least `minimum_length`.
/// Unless `no_max_prune` is set, a candidate is discarded if a matching pair can
/// extend either end within the query interval and target bounds. This checks
/// one adjacent pair, not whether a larger candidate passes every other rule.
///
/// Enumeration and eligibility are independent of production seeding; only
/// `Base`'s pairing semantics (including optional wobble) are shared.
pub(crate) fn naive_seeding(
    query: &[Base],
    target: &[Base],
    query_interval: Range<usize>,
    minimum_length: usize,
    config: &SeedConfig,
) -> Vec<(usize, usize, usize)> {
    let is_match = |query: Base, target: Base| query.pair_type(target).is_match(config.seed_wobble);
    let mut seeds = Vec::new();

    for query_start in query_interval.clone() {
        for target_start in 0..target.len() {
            let maximum_length =
                (query_interval.end - query_start).min(target.len() - target_start);
            for length in minimum_length..=maximum_length {
                let query_span = &query[query_start..query_start + length];
                let target_span = &target[target_start..target_start + length];
                // Ambiguity disqualifies a span even when mismatches are allowed.
                if query_span.contains(&Base::N) || target_span.contains(&Base::N) {
                    continue;
                }

                let paired: Vec<_> = query_span
                    .iter()
                    .zip(target_span)
                    .map(|(&query, &target)| is_match(query, target))
                    .collect();
                let mismatch_count = paired.iter().filter(|&&is_paired| !is_paired).count();
                if mismatch_count > config.max_mismatches {
                    continue;
                }

                if mismatch_count > 0 {
                    // Interrupted seeds fill gaps between perfect seeds. A full
                    // perfect seed inside the candidate makes it redundant.
                    let contains_perfect_seed = paired
                        .windows(minimum_length)
                        .any(|run| run.iter().all(|&is_paired| is_paired));
                    let matching_prefix = paired.iter().take_while(|&&is_paired| is_paired).count();
                    let matching_suffix = paired
                        .iter()
                        .rev()
                        .take_while(|&&is_paired| is_paired)
                        .count();
                    if contains_perfect_seed
                        || matching_prefix < config.min_prefix_matches
                        || matching_suffix < config.min_suffix_matches
                    {
                        continue;
                    }
                }

                if !config.no_max_prune {
                    let query_end = query_start + length;
                    let target_end = target_start + length;
                    let extends_left = query_start > query_interval.start
                        && target_start > 0
                        && is_match(query[query_start - 1], target[target_start - 1]);
                    let extends_right = query_end < query_interval.end
                        && target_end < target.len()
                        && is_match(query[query_end], target[target_end]);
                    if extends_left || extends_right {
                        continue;
                    }
                }

                seeds.push((query_start, target_start, length));
            }
        }
    }

    seeds
}

/// Compare every production seed with direct enumeration across all records.
///
/// Sorting removes traversal order only; duplicate entries remain observable.
/// Explicit interval expectations also check positive/negative config coordinates.
fn assert_seeds_match_enumeration(
    queries: &[String],
    targets: &[String],
    query_interval: Option<Range<usize>>,
    config: &SeedConfig,
) {
    let fasta = queries
        .iter()
        .enumerate()
        .map(|(i, q)| format!(">q{i}\n{q}\n"))
        .collect::<String>();
    let file = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(file.path(), fasta).unwrap();
    let queries_registry =
        QueryRegistry::build(read_sequences(file.path()).unwrap(), config).unwrap();
    let targets_registry = TargetRegistry::build(
        targets
            .iter()
            .enumerate()
            .map(|(i, t)| {
                let name = format!("t{i}");
                let (seq, _) = Sequence::normalize(&name, t.as_bytes()).unwrap();
                (name, seq)
            })
            .collect(),
        None,
    )
    .unwrap();
    let engine = super::SeedingEngine::new(&queries_registry, &targets_registry);
    for (query_index, query) in queries.iter().enumerate() {
        let (query_bases, _) = Sequence::normalize("query", query.as_bytes()).unwrap();
        let bounds = query_interval.clone().unwrap_or(0..query_bases.len());
        // The public default is six bases, clamped only for shorter queries.
        // An explicitly selected interval instead supplies its whole width.
        // Keep this literal independent of production's DEFAULT_SEED_LEN.
        let minimum_length = match (config.seed_length, query_interval.is_some()) {
            (Some(length), _) if length > 0 => (length as usize).min(query_bases.len()),
            (_, true) => bounds.len(),
            _ => 6.min(query_bases.len()),
        };
        let mut expected = Vec::new();
        for (target_index, target) in targets.iter().enumerate() {
            let (target_bases, _) = Sequence::normalize("target", target.as_bytes()).unwrap();
            for strand in [Strand::Forward, Strand::Reverse] {
                let physical = target_in_duplex_order(&target_bases, strand);
                expected.extend(
                    naive_seeding(
                        &query_bases,
                        &physical,
                        bounds.clone(),
                        minimum_length,
                        config,
                    )
                    .into_iter()
                    .map(|(query_start, target_start, length)| {
                        (
                            target_index,
                            char::from(strand),
                            query_start,
                            target_start,
                            length,
                        )
                    }),
                );
            }
        }
        let mut actual = Vec::new();
        engine
            .seed_query(query_index, config, |s| {
                assert_eq!(s.query_idx(), query_index);
                actual.push((
                    s.target_idx(),
                    char::from(s.strand()),
                    s.query_range().start,
                    s.target_range().start,
                    s.seed_len(),
                ));
            })
            .unwrap();
        expected.sort_unstable();
        actual.sort_unstable();
        assert_eq!(
            actual, expected,
            "query={query} query_interval={bounds:?} config={config:?}"
        );
    }
}

/// All nonempty A/C/G/U sequences up to `max_length`, plus fixed ambiguity cases.
///
/// One corpus exercises singleton and shared suffix intervals together.
/// The canonical corpus grows exponentially: `max_length = 3` gives 84 sequences.
/// The additional `N`, `ANU`, and `NCG` cases are independent of that limit.
fn corpus(max_length: u32) -> Vec<String> {
    let mut corpus = Vec::new();
    for len in 1..=max_length {
        for mut code in 0..4usize.pow(len) {
            corpus.push(
                (0..len)
                    .map(|_| {
                        let b = b"ACGU"[code % 4] as char;
                        code /= 4;
                        b
                    })
                    .collect(),
            );
        }
    }
    corpus.extend(["N", "ANU", "NCG"].map(str::to_owned));
    corpus
}

#[rstest::rstest]
#[case::unrestricted(0, 0)]
#[case::both_flanks(1, 1)]
#[case::suffix_only(0, 1)]
#[case::prefix_only(1, 0)]
fn production_seeding_matches_brute_force_oracle(
    #[values(0, 1, 2, 3, 4)] max_mismatches: usize,
    #[case] min_prefix_matches: usize,
    #[case] min_suffix_matches: usize,
    #[values(1, 2, 3)] seed_length: i64,
    #[values(false, true)] seed_wobble: bool,
    #[values(false, true)] no_max_prune: bool,
) {
    let corpus = corpus(4);
    assert_seeds_match_enumeration(
        &corpus,
        &corpus,
        None,
        &SeedConfig {
            seed_length: Some(seed_length),
            seed_wobble,
            no_max_prune,
            max_mismatches,
            min_prefix_matches,
            min_suffix_matches,
            ..Default::default()
        },
    );
}

#[test]
fn interval_edges_repeats_and_mismatch_runs_match_enumeration() {
    let queries = ["ACGUACGU", "AAAAAAAU", "GUGUGUGU", "ACNUGCAN"].map(str::to_owned);
    let targets = ["ACGUACGU", "AAAAAAAU", "UGUGUGUG", "UUUUUUUU", "ACNUGCAN"].map(str::to_owned);
    for (start, end, bounds) in [(2, 7, 1..7), (-7, -2, 1..7), (1, 8, 0..8)] {
        for length in [None, Some(0), Some(2), Some(3), Some(5)] {
            for wobble in [false, true] {
                for no_max_prune in [false, true] {
                    for (max_mismatches, prefix, suffix) in
                        [(0, 1, 0), (1, 0, 0), (2, 2, 1), (2, 1, 2)]
                    {
                        assert_seeds_match_enumeration(
                            &queries,
                            &targets,
                            Some(bounds.clone()),
                            &SeedConfig {
                                seed_start: Some(start),
                                seed_end: Some(end),
                                seed_length: length,
                                seed_wobble: wobble,
                                no_max_prune,
                                max_mismatches,
                                min_prefix_matches: prefix,
                                min_suffix_matches: suffix,
                            },
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn root_prepruning_keeps_the_reverse_strand_ac_seed() {
    let config = SeedConfig {
        seed_length: Some(2),
        ..Default::default()
    };
    assert_eq!(
        naive_seeding(&[Base::A, Base::C], &[Base::U, Base::G], 0..2, 2, &config),
        [(0, 0, 2)]
    );
    assert!(naive_seeding(&[Base::A, Base::C], &[Base::C, Base::A], 0..2, 2, &config).is_empty());
    assert_seeds_match_enumeration(&["AC".into()], &["AC".into()], None, &config);
}

#[test]
fn mismatch_candidates_exclude_a_complete_perfect_subseed() {
    let config = SeedConfig {
        seed_length: Some(2),
        max_mismatches: 1,
        min_prefix_matches: 0,
        no_max_prune: true,
        ..Default::default()
    };
    // Pairing is ||. (not an interrupted seed) versus |.| (an interrupted seed).
    assert!(!naive_seeding(
        &[Base::A, Base::A, Base::A],
        &[Base::U, Base::U, Base::C],
        0..3,
        2,
        &config
    )
    .contains(&(0, 0, 3)));
    assert!(naive_seeding(
        &[Base::A, Base::A, Base::A],
        &[Base::U, Base::C, Base::U],
        0..3,
        2,
        &config
    )
    .contains(&(0, 0, 3)));
    assert_seeds_match_enumeration(
        &["AAA".into()],
        &["CUU".into(), "UCU".into()],
        None,
        &config,
    );
}

/// Complementary runs, a mixed-base repeat, and an ambiguous-only sequence.
///
/// Random corpora need not contain long complementary runs. Since this corpus
/// supplies both queries and targets, G8/C8 guarantee perfect seeds up to the
/// longest seed length checked here, with repeated starts that exercise
/// maximality pruning. ACGUACGU adds a mixed-base repeat that pairs with its
/// reversal in duplex order. N is ambiguous only and must never seed.
#[rstest::rstest]
#[case::unrestricted(0, 0)]
#[case::both_flanks(1, 1)]
#[case::suffix_only(0, 1)]
#[case::prefix_only(1, 0)]
fn complementary_runs_and_ambiguity_match_enumeration(
    #[values(1, 2, 4, 8)] seed_length: i64,
    #[values(0, 1, 2, 3, 4)] max_mismatches: usize,
    #[case] min_prefix_matches: usize,
    #[case] min_suffix_matches: usize,
    #[values(false, true)] seed_wobble: bool,
    #[values(false, true)] no_max_prune: bool,
) {
    let corpus = ["GGGGGGGG", "CCCCCCCC", "ACGUACGU", "N"].map(str::to_owned);
    assert_seeds_match_enumeration(
        &corpus,
        &corpus,
        None,
        &SeedConfig {
            seed_length: Some(seed_length),
            seed_wobble,
            no_max_prune,
            max_mismatches,
            min_prefix_matches,
            min_suffix_matches,
            ..Default::default()
        },
    );
}

// Proptest reads the current directory to persist regression files, which Miri
// rejects under isolation (proptest-rs/proptest#253).
#[cfg(not(miri))]
mod generated {
    use super::*;
    use proptest::prelude::*;

    /// Generate a corpus, a feasible seed config, and optionally an interval within
    /// the first query, encoded with positive or negative coordinates.
    ///
    /// Resolution rejects a length wider than its interval and coordinates outside
    /// the query, so both are drawn in range rather than folded back into it, which
    /// keeps shrinking inside the feasible region. Only the first query carries an
    /// interval: fixed coordinates need not fit queries of other lengths. Since
    /// `seed_start` and `seed_end` reach production only through `SeedConfig::resolve`,
    /// one encoding per case covers the negative path without a second traversal.
    fn seeding_case(
    ) -> impl Strategy<Value = (Vec<String>, Vec<String>, Option<Range<usize>>, SeedConfig)> {
        let corpus = prop::collection::vec(prop::collection::vec(0usize..5, 1..33), 1..8).prop_map(
            |sequences| {
                sequences
                    .into_iter()
                    .map(|seq| seq.into_iter().map(|b| b"ACGUN"[b] as char).collect())
                    .collect::<Vec<String>>()
            },
        );
        // Mix feasible mismatch flanks with restrictive settings: the latter must
        // still allow perfect seeds, whose eligibility ignores flanks.
        let flanks = (1i64..=8, 0usize..=8, 0usize..=8, any::<bool>()).prop_map(
            |(length, prefix, suffix, constrained_flanks)| {
                if constrained_flanks {
                    let prefix = prefix % length as usize;
                    (length, prefix, suffix % (length as usize - prefix))
                } else {
                    (length, prefix, suffix)
                }
            },
        );

        (corpus, flanks, 0usize..=4, any::<bool>(), any::<bool>())
            .prop_flat_map(
                |(corpus, (length, prefix, suffix), max_mismatches, seed_wobble, no_max_prune)| {
                    let config = SeedConfig {
                        seed_length: Some(length),
                        max_mismatches,
                        min_prefix_matches: prefix,
                        min_suffix_matches: suffix,
                        seed_wobble,
                        no_max_prune,
                        ..Default::default()
                    };
                    let query_len = corpus[0].len();
                    let interval = (0..query_len).prop_flat_map(move |start| {
                        (Just(start), (start + 1)..=query_len, any::<bool>())
                    });
                    (
                        Just(corpus),
                        Just(config),
                        prop_oneof![Just(None::<(usize, usize, bool)>), interval.prop_map(Some)],
                    )
                },
            )
            .prop_map(|(corpus, config, interval)| {
                let Some((start, end, negative)) = interval else {
                    return (corpus.clone(), corpus, None, config);
                };
                let query_len = corpus[0].len() as i64;
                let config = SeedConfig {
                    seed_start: Some(if negative {
                        start as i64 - query_len
                    } else {
                        start as i64 + 1
                    }),
                    seed_end: Some(if negative {
                        end as i64 - 1 - query_len
                    } else {
                        end as i64
                    }),
                    seed_length: config
                        .seed_length
                        .map(|length| length.min((end - start) as i64)),
                    ..config
                };
                (corpus[..1].to_vec(), corpus, Some(start..end), config)
            })
    }

    proptest! {
        #[test]
        fn randomized_seeding_matches_brute_force_oracle(
            (queries, targets, bounds, config) in seeding_case(),
        ) {
            assert_seeds_match_enumeration(&queries, &targets, bounds, &config);
        }
    }
}
