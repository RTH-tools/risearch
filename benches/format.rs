//! Benchmark hit formatting across the output formats
//!
//! Run with: cargo bench --bench format

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::alignment::AlignColumn;
use risearch::config::OutputFormat;
use risearch::output::format::format_hit_into;
use risearch::search::SearchHit;
use risearch::types::{Base, Energy, Strand};

fn setup_hit() -> (SearchHit, Vec<Base>, Vec<Base>) {
    let query = vec![Base::U, Base::G, Base::C, Base::U, Base::G, Base::C];
    let target = vec![Base::G, Base::U, Base::C, Base::G, Base::U, Base::C];

    let alignment = query
        .iter()
        .copied()
        .zip(target.iter().copied())
        .map(|(q, t)| AlignColumn::new(q, t))
        .collect::<Vec<_>>()
        .into_boxed_slice();

    let hit = SearchHit {
        query_idx: 0,
        target_idx: 0,
        q_start: 0,
        q_end: 5,
        t_start: 0,
        t_end: 5,
        strand: Strand::Forward,
        energy: Energy::from_kcal(-10.5),
        alignment: Some(alignment),
    };

    (hit, query, target)
}

fn bench_format_hit(c: &mut Criterion) {
    let (hit, _q_seq, target) = setup_hit();

    let mut group = c.benchmark_group("hit_formatting");

    for format in [
        OutputFormat::Minimal,
        OutputFormat::Detailed,
        OutputFormat::Cigar,
        OutputFormat::BindingSite,
    ] {
        let mut out = Vec::with_capacity(1024);

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{:?}", format)),
            &format,
            |b, &fmt| {
                b.iter(|| {
                    out.clear();
                    format_hit_into(
                        &mut out,
                        &hit,
                        "test_query",
                        "test_target",
                        &target,
                        0..target.len(),
                        fmt,
                    );
                    black_box(&out);
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_format_hit);
criterion_main!(benches);
