//! Benchmark Gotoh extension over the DP grid: left, right, repeated extensions
//!
//! Run with: cargo bench --bench dp

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::dp::gotoh::Gotoh;
use risearch::dp::DpGrid;
use risearch::dsm::{DsmRegistry, ScoringModel};
use risearch::types::{Base, DsmId, Energy};

// ============================================================================
// SIMPLE LINEAR CONGRUENTIAL GENERATOR (LCG)
// ============================================================================
// A deterministic PRNG without external dependencies. Provides good coverage
// of RNA bases for realistic DP workloads.

struct SimpleLcg {
    state: u64,
}

impl SimpleLcg {
    const A: u64 = 6364136223846793005;
    const C: u64 = 1442695040888963407;

    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Generate next pseudo-random u64
    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(Self::A).wrapping_add(Self::C);
        self.state
    }

    /// Generate a random RNA base (A, G, C, U, or N).
    fn next_base(&mut self) -> Base {
        match (self.next() >> 32) % 5 {
            0 => Base::A,
            1 => Base::G,
            2 => Base::C,
            3 => Base::U,
            _ => Base::N,
        }
    }
}

// ============================================================================
// SEQUENCE GENERATION
// ============================================================================

/// Generate a pseudo-random RNA sequence as the dense symbol ranks the DP kernel
/// consumes.
///
/// Orienting a flank into this form — reversing for a left extension — belongs to
/// the search engine, so benches consume ranks directly rather than keeping a
/// second copy of that polarity rule. Direction reaches the kernel only through
/// the transposed scoring model.
fn generate_ranks(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = SimpleLcg::new(seed);
    (0..len).map(|_| rng.next_base().as_u8()).collect()
}

/// The `n` ranks ending at `anchor`, as a left extension's window would cover.
/// Clamped at the sequence start so the window shrinks rather than panicking.
fn before(src: &[u8], anchor: usize, n: usize) -> &[u8] {
    &src[(anchor + 1).saturating_sub(n)..=anchor]
}

/// The `n` ranks starting at `anchor`, as a right extension's window would cover.
/// Clamped at the sequence end so the window shrinks rather than panicking.
fn from(src: &[u8], anchor: usize, n: usize) -> &[u8] {
    let tail = &src[anchor..];
    &tail[..tail.len().min(n)]
}

/// Copy a window of ranks into a reused DP buffer and return its length.
///
/// A forward `copy_from_slice`, not a reproduction of the engine's fill: the
/// engine reads `Base` per byte and reverses for a left extension. Kept inside
/// the timed region so an O(n) load is still charged against the O(n²) kernel,
/// which is what these numbers are compared on.
fn load(dst: &mut [u8], src: &[u8]) -> usize {
    let len = src.len().min(dst.len());
    dst[..len].copy_from_slice(&src[..len]);
    len
}

// ============================================================================
// BENCHMARK GROUPS
// ============================================================================

fn bench_extend_left(c: &mut Criterion) {
    let mut group = c.benchmark_group("extend_left");
    let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
    let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
    let left_model = model.transpose();
    let gotoh = Gotoh::new(&left_model);

    for len in [10, 20, 30, 50].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(len), len, |b, &len| {
            let query = generate_ranks(100, 12345);
            let target = generate_ranks(100, 54321);

            let mut grid = DpGrid::new(200);
            let mut q_buf = [0u8; 256];
            let mut t_buf = [0u8; 256];

            b.iter(|| {
                let q_len = load(&mut q_buf, black_box(&query[..len]));
                let t_len = load(&mut t_buf, black_box(&target[..len]));
                let _ = gotoh.extend(&q_buf[..q_len], &t_buf[..t_len], &mut grid);
            });
        });
    }

    group.finish();
}

fn bench_extend_right(c: &mut Criterion) {
    let mut group = c.benchmark_group("extend_right");
    let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
    let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
    let gotoh = Gotoh::new(&model);

    for len in [10, 20, 30, 50].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(len), len, |b, &len| {
            let query = generate_ranks(100, 12345);
            let target = generate_ranks(100, 54321);
            let q_end = 49;
            let t_end = 49;

            let mut grid = DpGrid::new(200);
            let mut q_buf = [0u8; 256];
            let mut t_buf = [0u8; 256];

            b.iter(|| {
                let q_len = load(&mut q_buf, black_box(from(&query, q_end, len)));
                let t_len = load(&mut t_buf, black_box(from(&target, t_end, len)));
                let _ = gotoh.extend(&q_buf[..q_len], &t_buf[..t_len], &mut grid);
            });
        });
    }

    group.finish();
}

/// Many extensions benchmark: simulates a search workload with 100 extensions.
/// Each extension pair (left + right) is a typical search seed extension pattern.
fn bench_many_extensions(c: &mut Criterion) {
    let mut group = c.benchmark_group("many_extensions");
    group.sample_size(10);
    let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
    let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
    let left_model = model.transpose();
    let gotoh_right = Gotoh::new(&model);
    let gotoh_left = Gotoh::new(&left_model);

    group.bench_function("100_extensions_len30", |b| {
        let queries: Vec<Vec<u8>> = (0..10)
            .map(|i| generate_ranks(100, 1000 + i as u64))
            .collect();
        let targets: Vec<Vec<u8>> = (0..10)
            .map(|i| generate_ranks(100, 5000 + i as u64))
            .collect();

        let mut grid = DpGrid::new(200);
        let mut q_buf = [0u8; 256];
        let mut t_buf = [0u8; 256];

        b.iter(|| {
            let mut total_score = 0i32;

            // Simulate 100 search hit extensions (10 query × 10 target combinations)
            for (q_idx, query) in queries.iter().enumerate() {
                for (t_idx, target) in targets.iter().enumerate() {
                    // Anchor positions based on indices for variety
                    let q_start = (50 + q_idx) % query.len();
                    let t_start = (50 + t_idx) % target.len();
                    let q_end = if q_start + 30 < query.len() {
                        q_start + 30
                    } else {
                        query.len() - 1
                    };
                    let t_end = if t_start + 30 < target.len() {
                        t_start + 30
                    } else {
                        target.len() - 1
                    };

                    // Left extension: the window is the 30 symbols before the anchor.
                    let q_len = load(&mut q_buf, black_box(before(query, q_start, 30)));
                    let t_len = load(&mut t_buf, black_box(before(target, t_start, 30)));
                    let left_result =
                        gotoh_left.extend(&q_buf[..q_len], &t_buf[..t_len], &mut grid);
                    total_score = total_score.wrapping_add(left_result.score);

                    // Right extension: whatever remains after the anchor, up to 30.
                    let q_len = load(&mut q_buf, black_box(from(query, q_end, 30)));
                    let t_len = load(&mut t_buf, black_box(from(target, t_end, 30)));
                    let right_result =
                        gotoh_right.extend(&q_buf[..q_len], &t_buf[..t_len], &mut grid);
                    total_score = total_score.wrapping_add(right_result.score);
                }
            }

            black_box(total_score);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_extend_left,
    bench_extend_right,
    bench_many_extensions
);
criterion_main!(benches);
