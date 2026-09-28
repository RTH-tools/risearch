//! Benchmark seed detection: target and query scaling, mismatch seeding
//!
//! Run with: cargo bench --bench seed

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::config::SeedConfig;
use risearch::registry::QueryRegistry;
use risearch::seed::{SeedHit, SeedingEngine};
use risearch::seq::Sequence;
use risearch::types::Base;
use risearch::TargetRegistry;

struct SimpleLcg {
    state: u64,
}

impl SimpleLcg {
    const A: u64 = 6364136223846793005;
    const C: u64 = 1442695040888963407;

    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(Self::A).wrapping_add(Self::C);
        self.state
    }

    fn next_base(&mut self) -> Base {
        match (self.next() >> 32) % 4 {
            0 => Base::A,
            1 => Base::G,
            2 => Base::C,
            _ => Base::U,
        }
    }
}

struct ProductionSeedDataset {
    queries: QueryRegistry,
    store: TargetRegistry,
}

fn generate_sequence(len: usize, seed: u64) -> Sequence {
    let mut rng = SimpleLcg::new(seed);
    let bases: Vec<Base> = (0..len).map(|_| rng.next_base()).collect();
    Sequence::from(bases)
}

fn build_production_dataset(
    query_count: usize,
    query_len: usize,
    target_len: usize,
    seed_config: &SeedConfig,
) -> ProductionSeedDataset {
    let queries = (0..query_count)
        .map(|i| {
            (
                format!("q{i}"),
                generate_sequence(query_len, 1_000 + i as u64),
            )
        })
        .collect();
    let targets = [generate_sequence(target_len, 9_999)];

    let queries = QueryRegistry::build(queries, seed_config).expect("query registry");
    let named_targets = targets
        .iter()
        .enumerate()
        .map(|(idx, seq)| (format!("t{idx}"), seq.clone()))
        .collect();
    let store = TargetRegistry::build(named_targets, None).expect("build target index");

    ProductionSeedDataset { queries, store }
}

fn seed_count(groups: &[(usize, Vec<SeedHit>)]) -> usize {
    groups.iter().map(|(_, seeds)| seeds.len()).sum()
}

fn exact_seed_config(seed_length: i64) -> SeedConfig {
    SeedConfig {
        seed_length: Some(seed_length),
        min_prefix_matches: 1,
        seed_wobble: true,
        ..Default::default()
    }
}

fn mismatch_seed_config(
    seed_length: i64,
    max_mm: usize,
    prefix: usize,
    suffix: usize,
    wobble: bool,
) -> SeedConfig {
    SeedConfig {
        seed_length: Some(seed_length),
        seed_wobble: wobble,
        max_mismatches: max_mm,
        min_prefix_matches: prefix,
        min_suffix_matches: suffix,
        ..Default::default()
    }
}

fn bench_seed_exact(c: &mut Criterion) {
    let mut group = c.benchmark_group("exact_target_scaling");
    let seed_config = exact_seed_config(7);

    // Isolate exact seeding cost as target corpus size grows.
    for target_len in [1_000, 10_000, 100_000] {
        let dataset = build_production_dataset(1, 22, target_len, &seed_config);
        group.bench_with_input(
            BenchmarkId::from_parameter(target_len),
            &target_len,
            |b, _| {
                b.iter(|| {
                    let seeds =
                        SeedingEngine::new(black_box(&dataset.queries), black_box(&dataset.store))
                            .run(black_box(&seed_config))
                            .unwrap();
                    black_box(seed_count(&seeds));
                });
            },
        );
    }

    group.finish();
}

fn bench_seed_mismatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("mismatch");

    // Measure how allowing more mismatches changes seed enumeration cost.
    for max_mm in [1, 2, 3, 4, 5] {
        let seed_config = mismatch_seed_config(7, max_mm, 2, 2, false);
        let dataset = build_production_dataset(1, 22, 10_000, &seed_config);

        group.bench_with_input(BenchmarkId::from_parameter(max_mm), &max_mm, |b, _| {
            b.iter(|| {
                let seeds =
                    SeedingEngine::new(black_box(&dataset.queries), black_box(&dataset.store))
                        .run(black_box(&seed_config))
                        .unwrap();
                black_box(seed_count(&seeds));
            });
        });
    }

    group.finish();
}

fn bench_seed_prod_shaped_mismatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("prod_shaped_mismatch");
    group.sample_size(10);

    // Approximate the common production-shaped case: multiple short queries
    // against one large target corpus with wobble-enabled mismatch seeding.
    for max_mm in [1, 2, 3, 4, 5] {
        let seed_config = mismatch_seed_config(7, max_mm, 2, 2, true);
        let dataset = build_production_dataset(10, 22, 100_000, &seed_config);

        group.bench_with_input(BenchmarkId::from_parameter(max_mm), &max_mm, |b, _| {
            b.iter(|| {
                let seeds =
                    SeedingEngine::new(black_box(&dataset.queries), black_box(&dataset.store))
                        .run(black_box(&seed_config))
                        .unwrap();
                black_box(seed_count(&seeds));
            });
        });
    }

    group.finish();
}

fn bench_seed_query_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("query_scaling");
    group.sample_size(10);

    let seed_config = mismatch_seed_config(7, 1, 2, 2, true);

    // Hold target size fixed and measure how global query traversal scales
    // as more queries are packed into the combined query corpus.
    for query_count in [1, 10, 50, 100] {
        let dataset = build_production_dataset(query_count, 22, 100_000, &seed_config);
        group.bench_with_input(
            BenchmarkId::from_parameter(query_count),
            &query_count,
            |b, _| {
                b.iter(|| {
                    let seeds =
                        SeedingEngine::new(black_box(&dataset.queries), black_box(&dataset.store))
                            .run(black_box(&seed_config))
                            .unwrap();
                    black_box(seed_count(&seeds));
                });
            },
        );
    }

    group.finish();
}

fn bench_seed_long_mismatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("long_seed_mismatch");
    group.sample_size(10);

    // Stress test: 21nt seed on 22nt query with 2 mismatches.
    // This mirrors the scenario where C pulls ahead by 3x.
    let seed_config = mismatch_seed_config(21, 2, 2, 2, true);
    let dataset = build_production_dataset(1, 22, 100_000, &seed_config);

    group.bench_function("seed_21_mm_2_2", |b| {
        b.iter(|| {
            let seeds = SeedingEngine::new(black_box(&dataset.queries), black_box(&dataset.store))
                .run(black_box(&seed_config))
                .unwrap();
            black_box(seed_count(&seeds));
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_seed_exact,
    bench_seed_mismatch,
    bench_seed_prod_shaped_mismatch,
    bench_seed_query_scaling,
    bench_seed_long_mismatch
);
criterion_main!(benches);
