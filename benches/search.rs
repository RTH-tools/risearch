//! Benchmark the production-shaped end-to-end search pipeline
//!
//! Run with: cargo bench --bench search

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::fs;
use std::path::Path;

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rayon::ThreadPoolBuilder;
use risearch::config::{ExtendConfig, FilterConfig, ScoreConfig, SearchConfig, SeedConfig};
use risearch::fastx::read_sequences;
use risearch::registry::QueryRegistry;
use risearch::search::run_search;
use risearch::seed::SeedingEngine;
use risearch::seq::Sequence;
use risearch::types::{Base, DsmId, Strand};
use risearch::Energy;
use risearch::TargetRegistry;
use risearch::VecSink;
use tempfile::TempDir;

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

struct ProductionSearchDataset {
    _tmpdir: TempDir,
    queries: Vec<(String, Sequence)>,
    store: TargetRegistry,
}

fn generate_sequence(len: usize, seed: u64) -> Sequence {
    let mut rng = SimpleLcg::new(seed);
    let bases: Vec<Base> = (0..len).map(|_| rng.next_base()).collect();
    Sequence::from(bases)
}

fn write_fasta(path: &Path, prefix: &str, seqs: &[Sequence]) {
    let mut buf = Vec::new();
    for (idx, seq) in seqs.iter().enumerate() {
        buf.extend_from_slice(format!(">{prefix}{idx}\n").as_bytes());
        for &base in seq.iter() {
            buf.push(base.to_u8_upper());
        }
        buf.push(b'\n');
    }
    fs::write(path, buf).expect("write FASTA");
}

fn build_production_dataset(
    query_count: usize,
    query_len: usize,
    target_len: usize,
) -> ProductionSearchDataset {
    let tmpdir = TempDir::new().expect("tempdir");
    let queries_path = tmpdir.path().join("queries.fa");

    let queries: Vec<_> = (0..query_count)
        .map(|i| generate_sequence(query_len, 1_000 + i as u64))
        .collect();
    let mut targets: Vec<Vec<Base>> = [9_999, 19_999]
        .map(|seed| {
            generate_sequence(target_len, seed)
                .iter()
                .copied()
                .collect()
        })
        .into();
    for (i, query) in queries.iter().enumerate() {
        let start = 1_000 + i * 1_000;
        let reverse_complement = query[..]
            .iter()
            .rev()
            .map(|base| base.complement())
            .collect::<Vec<_>>();
        targets[0][start..start + query.len()].copy_from_slice(&reverse_complement);
        targets[1][start..start + query.len()].copy_from_slice(query);
    }
    let targets = targets.into_iter().map(Sequence::from).collect::<Vec<_>>();

    write_fasta(&queries_path, "q", &queries);

    let queries = read_sequences(&queries_path).expect("read queries");
    let named_targets = targets
        .iter()
        .enumerate()
        .map(|(idx, seq)| (format!("t{idx}"), seq.clone()))
        .collect();
    let store = TargetRegistry::build(named_targets, None).expect("build target index");

    ProductionSearchDataset {
        _tmpdir: tmpdir,
        queries,
        store,
    }
}

fn make_search_config(seed_config: &SeedConfig) -> SearchConfig {
    SearchConfig {
        seed: seed_config.clone(),
        score: ScoreConfig {
            dsm_id: DsmId::from("t04"),
            penalty: Energy::from_kcal(0.0),
            temperature: 37,
        },
        extend: ExtendConfig {
            max_extension: 20,
            build_alignment: false,
        },
        filter: FilterConfig {
            delta_g: Energy::MIN,
            seed_energy: Energy::from_kcal(0.0),
            no_dedup: false,
        },
    }
}

fn seed_counts(
    queries: &QueryRegistry,
    store: &TargetRegistry,
    config: &SeedConfig,
) -> (usize, usize) {
    let engine = SeedingEngine::new(queries, store);
    let (mut forward, mut reverse) = (0, 0);
    for query_idx in 0..queries.len() {
        engine
            .seed_query(query_idx, config, |hit| match hit.strand() {
                Strand::Forward => forward += 1,
                Strand::Reverse => reverse += 1,
            })
            .expect("stream seeds");
    }
    (forward, reverse)
}

fn bench_search_prod_shaped_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("prod_shaped");
    group.sample_size(20);

    let seed_config = SeedConfig {
        seed_start: None,
        seed_end: None,
        seed_length: Some(7),
        seed_wobble: true,
        no_max_prune: false,
        max_mismatches: 1,
        min_prefix_matches: 2,
        min_suffix_matches: 2,
    };
    let dataset = build_production_dataset(10, 22, 50_000);
    let prepared =
        QueryRegistry::build(dataset.queries.clone(), &seed_config).expect("prepare queries");
    let mut score_args = make_search_config(&seed_config);
    score_args.seed.no_max_prune = true;
    score_args.filter.no_dedup = true;
    let mut hit_args = make_search_config(&seed_config);
    hit_args.filter.delta_g = Energy::from_kcal(100.0);
    hit_args.extend.build_alignment = true;
    let pool = ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("single-thread benchmark pool");

    let (forward, reverse) = pool.install(|| seed_counts(&prepared, &dataset.store, &seed_config));
    let seed_total = forward + reverse;
    assert_eq!(
        (seed_total, forward, reverse),
        (118_701, 59_318, 59_383),
        "seed workload changed"
    );
    let setup_sink = VecSink::default();
    let retained = pool
        .install(|| run_search(&dataset.queries, &dataset.store, &hit_args, &setup_sink))
        .expect("run setup search");
    assert_eq!(retained, 110_823, "retained-hit workload changed");

    let case = "10q_x_2x50k";
    group.throughput(Throughput::Elements(seed_total as u64));

    group.bench_with_input(BenchmarkId::new("seed_stream", case), &case, |b, _| {
        b.iter(|| {
            let counts = pool.install(|| {
                seed_counts(
                    black_box(&prepared),
                    black_box(&dataset.store),
                    black_box(&seed_config),
                )
            });
            black_box(counts);
        });
    });

    // Reject only after scoring, with maximality pruning disabled, so every
    // collected seed exercises ungapped scoring and both DP directions.
    group.bench_with_input(BenchmarkId::new("score_all_seeds", case), &case, |b, _| {
        b.iter(|| {
            let sink = VecSink::default();
            let emitted = pool
                .install(|| {
                    run_search(
                        black_box(&dataset.queries),
                        black_box(&dataset.store),
                        black_box(&score_args),
                        &sink,
                    )
                })
                .expect("score every seed");
            black_box(emitted);
        });
    });

    group.throughput(Throughput::Elements(retained as u64));

    group.bench_with_input(BenchmarkId::new("trace_hits", case), &case, |b, _| {
        b.iter(|| {
            let sink = VecSink::default();
            let emitted = pool
                .install(|| {
                    run_search(
                        black_box(&dataset.queries),
                        black_box(&dataset.store),
                        black_box(&hit_args),
                        &sink,
                    )
                })
                .expect("run traceback search");
            black_box((emitted, sink.into_hits()));
        });
    });

    group.finish();
}

criterion_group!(benches, bench_search_prod_shaped_pipeline);
criterion_main!(benches);
