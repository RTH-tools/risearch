//! Search module - finds miRNA-target interactions.
//!
//! Pipeline (parallel over queries; each query is processed end-to-end by one
//! worker):
//! 1) Build the query's own suffix array, traverse it against the shared target
//!    SA to enumerate seeds
//! 2) Extend each seed (DP) into a final hit
//! 3) Collapse hits that share a final bounding box to the lowest-energy one
//!    (unless `--no-dedup`), then hand the query's hits to a [`HitSink`]

mod extension;
#[cfg(test)]
mod metamorphic_tests;
#[cfg(test)]
mod test_support;

use crate::error::{Error, Result};
use log::info;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::ops::Range;
use std::sync::Mutex;

use self::extension::{ExtensionEngine, SeedExtension};
use crate::alignment::{fingerprint_symbols, AlignColumn};
use crate::config::{SearchConfig, MAX_EXTENSION};
use crate::dp::MAX_EXT;
use crate::dsm::ScoringModel;
use crate::index::store::TargetRegistry;
use crate::index::view::TargetView;
use crate::registry::QueryRegistry;
use crate::seed::{SeedHit, SeedingEngine};
use crate::seq::Sequence;
use crate::types::{Base, Energy, Strand};

/// One accepted interaction: a span of a query paired against a span of a target.
///
/// Coordinates are 0-based and inclusive at both ends; the output formats add 1.
/// `t_start`/`t_end` always address the original FASTA sequence. Internally the
/// search uses physical duplex views ordered 3'->5' alongside the 5'->3' query;
/// that span is converted to FASTA coordinates once when the hit is assembled.
///
/// Bases and names are not carried, only indices into the queries passed to
/// [`run_search`] and into the target registry. [`query`](Self::query) slices a
/// supplied query sequence; [`target`](Self::target) resolves the paired span
/// through a [`TargetView`].
///
/// `alignment` is populated only under
/// [`ExtendConfig::build_alignment`](crate::ExtendConfig).
#[derive(Debug, Clone)]
pub struct SearchHit {
    /// Position in the queries passed to [`run_search`].
    pub query_idx: u32,
    /// Index into the target registry; see [`target_index`](Self::target_index).
    pub target_idx: u32,
    /// First paired query position.
    pub q_start: usize,
    /// Last paired query position.
    pub q_end: usize,
    /// First paired target position.
    pub t_start: usize,
    /// Last paired target position.
    pub t_end: usize,
    /// Target strand the duplex lies on.
    pub strand: Strand,
    /// Free energy of the duplex.
    pub energy: Energy,
    // Keep the variable-length traceback out of the densely stored hit. Compact
    // registry indices offset the boxed slice's metadata, preserving a 64-byte
    // SearchHit on 64-bit targets even when alignment building is disabled.
    /// Resolved duplex columns; `None` unless traceback ran.
    pub alignment: Option<Box<[AlignColumn]>>,
}

#[cfg(target_pointer_width = "64")]
const _: () = assert!(std::mem::size_of::<SearchHit>() == 64);

/// Where a search's hits go.
///
/// The driver runs the seed -> extend -> dedup pipeline and hands each query's
/// surviving hits to one of these; the destination (one text file, one file per
/// query, Arrow columns) is entirely the implementor's business.
///
/// `consume` is called at most once per `query_idx`, from a rayon worker, in
/// arbitrary order — hence `&self` plus interior mutability. Lock once per call,
/// never per hit, and do any expensive conversion before taking the lock.
pub trait HitSink: Sync {
    /// Take one query's surviving hits.
    fn consume(&self, query_idx: usize, hits: Vec<SearchHit>) -> Result<()>;

    /// Finalize the destination. The driver calls this once every query has been
    /// consumed and the search succeeded — including when no hits arrived at all,
    /// which is what lets a zero-hit run truncate a stale output file. A run with
    /// nothing to search (empty target or query set) returns without consuming or
    /// flushing anything.
    fn flush(&self) -> Result<()> {
        Ok(())
    }
}

/// Collects every hit into one `Vec`. Peak RAM scales with the hit count, so
/// prefer a sink that writes through for large searches.
#[derive(Default)]
pub struct VecSink(Mutex<Vec<SearchHit>>);

impl VecSink {
    /// Unwrap the collected hits.
    pub fn into_hits(self) -> Vec<SearchHit> {
        self.0.into_inner().unwrap()
    }
}

impl HitSink for VecSink {
    fn consume(&self, _query_idx: usize, hits: Vec<SearchHit>) -> Result<()> {
        self.0.lock().unwrap().extend(hits);
        Ok(())
    }
}

/// Collapse hits sharing a bounding box (target, strand and extended span) to
/// the single best one under [`SearchHit::cmp`]. Caller guarantees all hits
/// belong to one query, so this is globally exact for that query.
fn dedup_hits(hits: Vec<SearchHit>) -> Vec<SearchHit> {
    use std::collections::hash_map::Entry;
    let mut best = HashMap::with_capacity(hits.len());
    for hit in hits {
        let key = (
            hit.target_idx,
            hit.strand,
            hit.q_start,
            hit.q_end,
            hit.t_start,
            hit.t_end,
        );
        match best.entry(key) {
            Entry::Occupied(mut e) => {
                if hit.cmp(e.get()).is_lt() {
                    e.insert(hit);
                }
            }
            Entry::Vacant(e) => {
                e.insert(hit);
            }
        }
    }
    best.into_values().collect()
}

/// Run search, handing every query's hits to `sink` as they are produced.
/// Returns the number of hits handed over.
///
/// Parallelizes over queries: each rayon worker owns one query end-to-end —
/// seed it against the shared target SA, extend, dedup — then hands the hits to
/// the sink and frees its per-query state. Peak RAM is therefore bounded by the
/// concurrent workers' per-query state plus whatever the sink retains, not by the
/// total hit count.
///
/// Queries arrive at the sink in completion order, and hits within a query in
/// `HashMap` order, so nothing about the output order is stable across runs. The
/// hit *set* is. Callers needing an order must sort.
///
/// `queries` are `(name, sequence)` records, e.g. from
/// [`read_sequences`](crate::fastx::read_sequences), which guarantees what is
/// required here: at least one record, and every sequence non-empty. They are
/// prepared against `opts.seed` here, and duplicate names are rejected. The `query_idx` a sink
/// receives, like [`SearchHit::query_idx`], is a position in `queries`.
///
/// The config is validated before any query runs, so a sink that defers touching
/// its destination until first use will not have disturbed it if this returns an
/// error.
pub fn run_search(
    queries: &[(String, Sequence)],
    targets: &TargetRegistry,
    opts: &SearchConfig,
    sink: &dyn HitSink,
) -> Result<usize> {
    opts.validate()?;

    let queries = QueryRegistry::build(queries.to_vec(), &opts.seed)?;
    // Deliberately not flushed: an empty index leaves the destination
    // untouched, where a zero-hit run truncates it.
    if targets.is_empty() {
        return Ok(0);
    }

    check_unlimited_fits(&queries, opts)?;
    log_search_banner(&queries, targets, opts);

    let ctx = SearchContext {
        queries: &queries,
        targets,
        opts,
        model: ScoringModel::load(
            &opts.score.dsm_id,
            opts.score.temperature,
            opts.score.penalty,
        )?,
    };
    let engine = SeedingEngine::new(ctx.queries, ctx.targets);
    let counts: Vec<usize> = (0..ctx.queries.len())
        .into_par_iter()
        .map_init(
            || SearchWorker::new(ctx.opts, &ctx.model),
            |worker, qi| -> Result<usize> {
                let hits = worker.search_query(&ctx, &engine, qi)?;
                let emitted = hits.len();
                sink.consume(qi, hits)?;
                Ok(emitted)
            },
        )
        .collect::<Result<Vec<usize>>>()?;
    sink.flush()?;

    let total = counts.iter().sum();
    info!("Search complete: {} hits", total);
    Ok(total)
}

struct SearchContext<'a> {
    queries: &'a QueryRegistry,
    targets: &'a TargetRegistry,
    opts: &'a SearchConfig,
    model: ScoringModel,
}

/// Unlimited extension (`-l -1`) promises to span the whole query, but the DP
/// buffers cap each side at MAX_EXT. Refuse rather than silently clamp: a query
/// longer than the cap cannot be served as requested. (Mirrors clap rejecting an
/// explicit `-l > MAX_EXTENSION`.)
fn check_unlimited_fits(queries: &QueryRegistry, opts: &SearchConfig) -> Result<()> {
    if !opts.extend.is_unlimited() {
        return Ok(());
    }
    for q in queries.iter() {
        let n = q.sequence().len();
        if n > MAX_EXT {
            return Err(Error::Config(format!(
                "query '{}' is {n} nt; `-l -1` cannot extend across it ({MAX_EXT} nt cap). \
                 Pass an explicit `-l <={MAX_EXTENSION}` to accept the cap, or shorten the query.",
                q.name()
            )));
        }
    }
    Ok(())
}

fn log_search_banner(queries: &QueryRegistry, targets: &TargetRegistry, opts: &SearchConfig) {
    let max_ext = if opts.extend.is_unlimited() {
        format!("unlimited(<={MAX_EXT})")
    } else {
        opts.extend.max_extension.to_string()
    };
    info!(
        "Starting search: {} queries x {} targets, seed_length={:?}, max_ext={}, delta_g={}",
        queries.len(),
        targets.len(),
        opts.seed.seed_length,
        max_ext,
        opts.filter.delta_g
    );
}

struct SearchWorker {
    extension: ExtensionEngine,
}

impl SearchWorker {
    fn new(opts: &SearchConfig, model: &ScoringModel) -> Self {
        Self {
            extension: ExtensionEngine::new(
                opts.extend.max_window(),
                opts.extend.build_alignment,
                model,
            ),
        }
    }

    /// Seed one query, extend and score every seed, and collapse the survivors.
    fn search_query(
        &mut self,
        ctx: &SearchContext<'_>,
        engine: &SeedingEngine<'_>,
        query_idx: usize,
    ) -> Result<Vec<SearchHit>> {
        let opts = ctx.opts;
        let query = ctx.queries[query_idx].sequence();

        // Loop-invariant: rebuilding it per seed re-enters the rkyv root twice.
        let tview = ctx.targets.view();
        let mut hits = Vec::new();
        engine.seed_query(query_idx, &opts.seed, |seed| {
            let target = tview.target(seed.target_idx(), seed.strand());

            let ext = self.extension.extend_seed(query, target, &seed);

            if ext.energy <= opts.filter.delta_g {
                let alignment = self.extension.materialize_alignment(query, target, &seed);
                hits.push(SearchHit::new(query_idx, &seed, ext, alignment, tview));
            }
        })?;
        // Dedup is exact per query: every box-mate of `query_idx` is in `hits`.
        Ok(if opts.filter.no_dedup {
            hits
        } else {
            dedup_hits(hits)
        })
    }
}

/// Convert a non-empty internal half-open range to public inclusive bounds.
fn inclusive_bounds(mut range: Range<usize>) -> (usize, usize) {
    let start = range.start;
    let end = range
        .next_back()
        .expect("extended duplex ranges must be non-empty");
    (start, end)
}

impl SearchHit {
    /// Compare two candidates for the same bounding box.
    ///
    /// `Less` means `self` is retained, `Greater` means `other` is retained, and
    /// `Equal` means they are indistinguishable under the deduplication policy.
    /// Energy orders first; an exact tie is broken by the alignment fingerprint.
    fn cmp(&self, other: &Self) -> Ordering {
        self.energy.cmp(&other.energy).then_with(|| {
            match (self.alignment.as_deref(), other.alignment.as_deref()) {
                (Some(a), Some(b)) => fingerprint_symbols(a).cmp(fingerprint_symbols(b)),
                (None, None) => Ordering::Equal,
                (None, Some(_)) => Ordering::Less,
                (Some(_), None) => Ordering::Greater,
            }
        })
    }

    /// Slice the paired span out of a query sequence.
    pub fn query<'a>(&self, q_seq: &'a [Base]) -> &'a [Base] {
        &q_seq[self.q_start..=self.q_end]
    }

    /// Resolve the paired span against the strand-selected target, in duplex-column order.
    pub fn target<'a>(&self, targets: TargetView<'a>) -> &'a [Base] {
        let target = targets.target(self.target_index(), self.strand);
        &target[self.duplex_target_range(targets)]
    }

    /// Map the public inclusive FASTA coordinates into a half-open range over
    /// the strand-selected duplex target.
    pub(crate) fn duplex_target_range(&self, targets: TargetView<'_>) -> Range<usize> {
        let fasta_range = self.t_start
            ..self
                .t_end
                .checked_add(1)
                .expect("inclusive target end must fit a half-open range");
        targets.map_target_range(self.target_index(), self.strand, fasta_range)
    }

    /// [`target_idx`](Self::target_idx) as an index into the target registry.
    #[inline]
    pub fn target_index(&self) -> usize {
        usize::try_from(self.target_idx).expect("u32 target index must fit usize")
    }

    /// Convert the extension's duplex-frame target span to FASTA coordinates and
    /// record the hit.
    fn new(
        query_idx: usize,
        seed: &SeedHit,
        ext: SeedExtension,
        alignment: Option<Box<[AlignColumn]>>,
        targets: TargetView<'_>,
    ) -> Self {
        let SeedExtension {
            q_range,
            t_range,
            energy,
        } = ext;
        let (q_start, q_end) = inclusive_bounds(q_range);
        let fasta_t_range = targets.map_target_range(seed.target_idx(), seed.strand(), t_range);
        let (t_start, t_end) = inclusive_bounds(fasta_t_range);

        Self {
            query_idx: u32::try_from(query_idx).expect("query index exceeds u32"),
            target_idx: u32::try_from(seed.target_idx()).expect("target index exceeds u32"),
            q_start,
            q_end,
            t_start,
            t_end,
            strand: seed.strand(),
            energy,
            alignment,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use super::test_support::{config, records, search, QUERIES, TARGETS};
    use super::*;
    use crate::config::{OutputCompression, OutputConfig, OutputFormat};
    use crate::index::store::TargetRegistry;
    use crate::output::TextSink;
    use crate::types::DsmId;
    use crate::Sequence;
    use rstest::rstest;

    fn detailed_output() -> OutputConfig {
        OutputConfig {
            format: OutputFormat::Detailed,
            compress: OutputCompression::None,
            multifile: false,
        }
    }

    fn minimal_output() -> OutputConfig {
        OutputConfig {
            format: OutputFormat::Minimal,
            ..detailed_output()
        }
    }

    #[test]
    fn text_output_emits_one_line_per_hit() {
        let output = minimal_output();
        // Mirror what the CLI derives for Minimal, so this also covers dedup's
        // empty-fingerprint (first-wins) tie-break.
        let config = config(|c| c.extend.build_alignment = false);

        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        let hits = search_to_path(&QUERIES, &TARGETS, &config, &output, out.path());
        let file_hit_count = fs_err::read_to_string(out.path())
            .unwrap()
            .lines()
            .filter(|line| !line.trim().is_empty())
            .count();

        assert_eq!(
            hits, file_hit_count,
            "text output must emit exactly one line per retained hit"
        );
    }

    #[test]
    fn no_hits_below_an_unreachable_cutoff() {
        let config = config(|c| c.filter.delta_g = Energy::from_kcal(-100.0));
        let hits = search(&QUERIES, &TARGETS, &config).unwrap();

        assert!(hits.is_empty(), "no duplex reaches -100 kcal/mol");
    }

    #[test]
    fn run_search_validates_configs_constructed_without_clap() {
        let bad_penalty = config(|c| c.score.penalty = Energy::from_kcal(-0.1));
        let err = search(&QUERIES, &TARGETS, &bad_penalty).unwrap_err();
        assert!(err.to_string().contains("penalty"));

        let bad_window = config(|c| c.extend.max_extension = -2);
        let err = search(&QUERIES, &TARGETS, &bad_window).unwrap_err();
        assert!(err.to_string().contains("max extension"));
    }

    /// Zero-hit behaviour differs by topology and neither branch is covered by a
    /// CLI test: multifile creates no file at all, single-file truncates the
    /// output to nothing.
    #[test]
    fn zero_hit_run_writes_empty_file_but_no_multifile_entry() {
        let config = config(|c| c.filter.delta_g = Energy::from_kcal(-100.0));
        let mut output = detailed_output();

        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        fs_err::write(out.path(), b"stale").unwrap();
        search_to_path(&QUERIES, &TARGETS, &config, &output, out.path());
        assert_eq!(fs_err::metadata(out.path()).unwrap().len(), 0);

        output.multifile = true;
        let tmp = tempfile::tempdir().unwrap();
        // Not pre-created, so this also pins that the directory itself is made.
        let dir = tmp.path().join("multi");
        search_to_path(&QUERIES, &TARGETS, &config, &output, &dir);
        assert_eq!(fs_err::read_dir(&dir).unwrap().count(), 0);
    }

    #[test]
    fn compressed_output_decodes_to_the_plain_rows_including_empty_results() {
        let mut config = collision_config();
        for delta_g in [-5.0, -1000.0] {
            config.filter.delta_g = Energy::from_kcal(delta_g);
            let mut plain = search_to_lines(&QUERIES, &TARGETS, &config);
            plain.sort();
            assert_eq!(!plain.is_empty(), delta_g > -1000.0);
            for compress in [OutputCompression::Gzip(6), OutputCompression::Zstd(3)] {
                let output = OutputConfig {
                    compress,
                    ..minimal_output()
                };
                let out = tempfile::NamedTempFile::new().unwrap();
                search_to_path(&QUERIES, &TARGETS, &config, &output, out.path());
                let bytes = fs_err::read(out.path()).unwrap();
                let decoded = match compress {
                    OutputCompression::Gzip(_) => {
                        assert!(bytes.starts_with(&[0x1f, 0x8b]));
                        let mut text = String::new();
                        flate2::read::MultiGzDecoder::new(bytes.as_slice())
                            .read_to_string(&mut text)
                            .unwrap();
                        text
                    }
                    OutputCompression::Zstd(_) => {
                        assert!(bytes.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]));
                        String::from_utf8(zstd::stream::decode_all(bytes.as_slice()).unwrap())
                            .unwrap()
                    }
                    OutputCompression::None => unreachable!(),
                };
                let mut rows: Vec<String> = decoded
                    .lines()
                    .filter(|l| !l.trim().is_empty())
                    .map(str::to_string)
                    .collect();
                rows.sort();
                assert_eq!(rows, plain, "{compress:?} delta_g={delta_g}");
            }
        }
    }

    /// The sink is built before the search validates the config, so it must not
    /// touch its destination until the driver hands it something.
    #[test]
    fn rejected_config_leaves_no_multifile_directory() {
        let config = config(|c| c.extend.max_extension = -1);
        let output = OutputConfig {
            multifile: true,
            ..detailed_output()
        };
        let queries = records(&[("longq", &"A".repeat(MAX_EXT + 1))]);

        let tmp = tempfile::tempdir().unwrap();
        let dir = tmp.path().join("multi");
        let sink = TextSink::new(&queries, &TARGETS, &output, &dir).unwrap();
        assert!(run_search(&queries, &TARGETS, &config, &sink).is_err());
        assert!(!dir.exists());
    }

    /// Multifile addresses its writers by index (`paths[query_idx]`) while queries
    /// that produce nothing get no file at all, so a query dropping out in the
    /// middle must not shift the others' rows into the wrong file.
    #[test]
    fn multifile_writes_each_querys_rows_to_its_own_file() {
        let queries = records(&[("q1", "AAAA"), ("q2", "CCCC"), ("q3", "AAAAA")]);
        let targets = TargetRegistry::build(records(&[("t", "UUUUUUUU")]), None).unwrap();
        let config = config(|c| {
            c.seed.seed_length = Some(4);
            c.extend.max_extension = 0;
            c.extend.build_alignment = false;
            c.filter.delta_g = Energy::from_kcal(100.0);
        });
        let output = OutputConfig {
            multifile: true,
            ..minimal_output()
        };

        let run_with_workers = |workers: usize| {
            let tmp = tempfile::tempdir().unwrap();
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| search_to_path(&queries, &targets, &config, &output, tmp.path()));
            let mut files: Vec<(String, String)> = fs_err::read_dir(tmp.path())
                .unwrap()
                .map(|entry| {
                    let path = entry.unwrap().path();
                    let stem = path.file_stem().unwrap().to_string_lossy().into_owned();
                    let mut lines: Vec<String> = fs_err::read_to_string(&path)
                        .unwrap()
                        .lines()
                        .map(str::to_owned)
                        .collect();
                    lines.sort();
                    (stem, lines.join("\n"))
                })
                .collect();
            files.sort();
            files
        };
        let files = run_with_workers(1);
        assert_eq!(
            files,
            run_with_workers(4),
            "multifile output changed across worker counts"
        );

        let stems: Vec<&str> = files.iter().map(|(stem, _)| stem.as_str()).collect();
        assert_eq!(stems, ["q1", "q3"], "only queries with hits get a file");
        for (stem, body) in &files {
            assert!(!body.trim().is_empty(), "{stem} file is empty");
            for line in body.lines() {
                assert_eq!(
                    line.split('\t').next().unwrap(),
                    stem,
                    "row landed in the wrong query's file"
                );
            }
        }
    }

    /// The [`HitSink`] contract the in-tree sinks depend on: every query
    /// is consumed exactly once, and every hit in a batch carries the
    /// `query_idx` handed alongside it (multifile picks its file from that index,
    /// single-file looks up the query's name and sequence).
    #[test]
    fn streaming_consumes_each_query_once_with_matching_hit_indices() {
        #[derive(Default)]
        struct Spy(Mutex<Vec<usize>>);

        impl HitSink for Spy {
            fn consume(&self, query_idx: usize, hits: Vec<SearchHit>) -> Result<()> {
                assert!(hits
                    .iter()
                    .all(|hit| usize::try_from(hit.query_idx).unwrap() == query_idx));
                self.0.lock().unwrap().push(query_idx);
                Ok(())
            }
        }

        let config = collision_config();
        let spy = Spy::default();
        run_search(&QUERIES, &TARGETS, &config, &spy).unwrap();

        let mut seen = spy.0.into_inner().unwrap();
        seen.sort_unstable();
        assert_eq!(seen, (0..QUERIES.len()).collect::<Vec<_>>());
    }

    /// Search into a text file at `path` (a directory under `--multifile`),
    /// returning the hit count. The sink dies here, so the caller reads a file
    /// whose compression trailer is already written.
    fn search_to_path(
        queries: &[(String, Sequence)],
        targets: &TargetRegistry,
        config: &SearchConfig,
        output: &OutputConfig,
        path: &std::path::Path,
    ) -> usize {
        let sink = TextSink::new(queries, targets, output, path).unwrap();
        run_search(queries, targets, config, &sink).unwrap()
    }

    fn search_to_lines(
        queries: &[(String, Sequence)],
        targets: &TargetRegistry,
        config: &SearchConfig,
    ) -> Vec<String> {
        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        search_to_path(queries, targets, config, &minimal_output(), out.path());
        fs_err::read_to_string(out.path())
            .unwrap()
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(str::to_string)
            .collect()
    }

    /// The fields dedup collapses on, within one query.
    fn hit_box(hit: &SearchHit) -> (u32, u32, Strand, usize, usize, usize, usize) {
        (
            hit.query_idx,
            hit.target_idx,
            hit.strand,
            hit.q_start,
            hit.q_end,
            hit.t_start,
            hit.t_end,
        )
    }

    fn collision_config() -> SearchConfig {
        // Parameters that produce overlapping-seed box collisions on the
        // corpus (~33% redundant rows): short seed, mismatches, wide window,
        // permissive energy.
        config(|c| {
            c.seed.seed_length = Some(6);
            c.seed.seed_wobble = true;
            c.seed.max_mismatches = 2;
            c.seed.min_suffix_matches = 0;
            c.extend.max_extension = 30;
            c.extend.build_alignment = false;
            c.filter.delta_g = Energy::from_kcal(-5.0);
        })
    }

    /// Default dedup keeps exactly one hit per bounding box, and that hit carries
    /// the minimum energy among the box's collapsed candidates. `--no-dedup`
    /// reproduces the full per-seed hit set.
    #[test]
    fn dedup_keeps_one_min_energy_hit_per_box() {
        use std::collections::HashMap;
        let config = collision_config();
        let mut nodedup = config.clone();
        nodedup.filter.no_dedup = true;

        let raw = search(&QUERIES, &TARGETS, &nodedup).unwrap();
        let dedup = search(&QUERIES, &TARGETS, &config).unwrap();
        assert!(!dedup.is_empty(), "corpus must yield hits");
        assert!(
            dedup.len() < raw.len(),
            "dedup must remove overlapping hits"
        );

        let mut min_energy = HashMap::new();
        for hit in &raw {
            min_energy
                .entry(hit_box(hit))
                .and_modify(|energy| *energy = hit.energy.min(*energy))
                .or_insert(hit.energy);
        }
        let kept: HashMap<_, _> = dedup.iter().map(|hit| (hit_box(hit), hit.energy)).collect();

        assert_eq!(kept.len(), dedup.len(), "duplicate box survived dedup");
        assert_eq!(kept, min_energy, "each box must keep its minimum energy");
    }

    /// The deduped set must be stable across runs. `dedup_hits` drains a HashMap,
    /// so a tie-break that depended on hash/iteration order would surface here as
    /// run-to-run drift (two HashMaps use different random seeds). This covers
    /// behavior raw-hit comparisons cannot exercise because they bypass deduplication.
    #[test]
    fn dedup_result_is_deterministic() {
        let config = collision_config();
        let boxes = || {
            let hits = search(&QUERIES, &TARGETS, &config).unwrap();
            let mut boxes: Vec<_> = hits.iter().map(|hit| (hit_box(hit), hit.energy)).collect();
            boxes.sort();
            boxes
        };

        let first = boxes();
        assert!(!first.is_empty(), "corpus must yield hits");
        assert_eq!(
            first,
            boxes(),
            "deduped set must not depend on hash/iteration order"
        );
    }

    // Unlimited seed extension (`-l -1`).
    //
    // risearch2 (C) cannot serve as an oracle: it clamps `-l` with
    // `MAX(0, atoi(optarg))`, so `-l -1` silently becomes `0` (seed-only). The
    // unlimited window is a risearch3-only feature, so these tests use risearch3
    // itself as the oracle via metamorphic relations.
    //
    // Note on what is *not* tested: `-l -1` is **not** equivalent to a large fixed
    // window such as `-l 255`. Unlimited sizes each extension window to the query
    // bases available on that side (capped at the `dp::MAX_EXT` = 256 buffer
    // ceiling), whereas a fixed `-l k` permits up to `k` of extension including
    // large target-side bulges. The two therefore diverge on real data, so there
    // is no fixed-window oracle to compare against. Instead we assert the
    // feature's actual guarantees: reach and rejection.
    //
    // Extension only does work when the optimal duplex extends past the exact
    // seed, so the reach test breaks exact complementarity with a single mismatch
    // near the 3' end: the seeder stops there, and only DP extension can bridge it.

    fn complement(base: char) -> char {
        char::from(
            Base::try_from(base)
                .unwrap()
                .complement()
                .to_byte()
                .to_ascii_uppercase(),
        )
    }

    /// Reverse complement, so `revcomp(q)` forms a full antiparallel duplex with `q`.
    fn revcomp(seq: &str) -> String {
        seq.chars().rev().map(complement).collect()
    }

    /// Deterministic pseudo-random RNA (LCG), for inputs too long to write by hand.
    fn pseudo_random_rna(len: usize, seed: u64) -> String {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                match (state >> 33) % 4 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'U',
                }
            })
            .collect()
    }

    /// Reach: extension crosses a mismatch the exact-match seeder cannot, reaching
    /// the query 3' end — something seed-only (`-l 0`) never does.
    #[test]
    fn unlimited_extends_past_the_exact_seed_to_the_query_end() {
        // Perfect 30-nt duplex, then one base near the 3' end mutated to break exact
        // matching there. The seeder stops before the mismatch; only DP extension
        // can bridge it and pair the final bases.
        const N: usize = 30;
        const BREAK: usize = 27; // mismatch position; indices 28,29 stay complementary

        let mut query: Vec<char> = pseudo_random_rna(N, 0x07).chars().collect();
        let target = revcomp(&query.iter().collect::<String>());
        // In an antiparallel duplex query index i pairs target index N-1-i. Setting
        // the query base equal to its opposing target base guarantees a mismatch
        // (identical bases never pair), without touching neighboring positions.
        let opposing = target.as_bytes()[N - 1 - BREAK] as char;
        query[BREAK] = opposing;
        let query: String = query.into_iter().collect();

        // The default penalty of zero keeps added base pairs favorable, so
        // extension actually runs.
        let with_window = |window| {
            let config = config(|c| {
                c.seed.seed_length = Some(7);
                c.extend.max_extension = window;
                c.filter.delta_g = Energy::from_kcal(-8.0);
            });
            search_pair(&query, &target, &config).unwrap()
        };
        let seed_only = with_window(0);
        let unlimited = with_window(-1);

        assert!(!seed_only.is_empty(), "seed-only search must produce a hit");
        assert!(
            seed_only.iter().all(|h| h.q_end < N - 1),
            "seed-only (`-l 0`) must stop at the mismatch, not reach the 3' end; got {:#?}",
            seed_only
        );
        assert!(
            unlimited
                .iter()
                .any(|h| h.q_start == 0 && h.q_end == N - 1 && h.energy.to_kcal() < 0.0),
            "unlimited extension must cross the mismatch and reach the query 3' end; got {:#?}",
            unlimited
        );
    }

    /// Rejection: a query longer than the 256-nt buffer ceiling (`dp::MAX_EXT`) is
    /// refused up front rather than silently clamped — `-l -1` cannot honor
    /// "span the whole query" past the cap.
    #[test]
    fn unlimited_rejects_a_query_longer_than_the_ceiling() {
        const N: usize = MAX_EXT + 1;

        let query = pseudo_random_rna(N, 0x5151_2323);
        let target = revcomp(&query);

        let err = search_pair(&query, &target, &config(|c| c.extend.max_extension = -1))
            .expect_err("`-l -1` must reject a query longer than the 256-nt cap");

        assert!(
            err.to_string().contains("cannot extend across it"),
            "unexpected error message: {err}"
        );
    }

    /// The ceiling is inclusive: a query of exactly `MAX_EXT` is still servable,
    /// so the rejection must not fire one nucleotide early.
    #[test]
    fn unlimited_accepts_a_query_of_exactly_the_ceiling() {
        let query = pseudo_random_rna(MAX_EXT, 0x5151_2323);
        let target = revcomp(&query);

        search_pair(&query, &target, &config(|c| c.extend.max_extension = -1))
            .expect("a query of exactly MAX_EXT must be accepted");
    }

    // ── Coordinate conversion: oracle + planted duplexes ─────────────────

    /// Search one raw query against one raw target.
    fn search_pair(
        query_text: &str,
        target_text: &str,
        config: &SearchConfig,
    ) -> Result<Vec<SearchHit>> {
        let targets = TargetRegistry::build(records(&[("t", target_text)]), None).unwrap();
        search(&records(&[("q", query_text)]), &targets, config)
    }

    fn plant_site(bg: char, total_len: usize, offset: usize, site: &str) -> String {
        let mut target: Vec<char> = std::iter::repeat_n(bg, total_len).collect();
        for (i, c) in site.chars().enumerate() {
            target[offset + i] = c;
        }
        target.into_iter().collect()
    }

    /// A site planted at a known offset comes back at exactly that FASTA span,
    /// on the strand it was planted for.
    #[rstest]
    fn planted_duplex_lands_at_its_fasta_position(
        #[values(Strand::Forward, Strand::Reverse)] strand: Strand,
        #[values((50, 0), (50, 5), (50, 20), (50, 43), (30, 3), (80, 60), (15, 8))] placement: (
            usize,
            usize,
        ),
    ) {
        let query = "UGCAUGU";
        let (target_len, offset) = placement;
        // Forward FASTA is reverse(physical); reverse FASTA is complement(physical).
        let site = match strand {
            Strand::Forward => revcomp(query),
            Strand::Reverse => query.to_string(),
        };
        let target = plant_site('A', target_len, offset, &site);
        let config = config(|c| {
            c.seed.seed_length = Some(query.len() as i64);
            c.seed.no_max_prune = true;
            c.extend.max_extension = 0;
            c.filter.delta_g = Energy::from_kcal(100_000.0);
        });

        let hits = search_pair(query, &target, &config).unwrap();

        let span = (strand, offset, offset + query.len() - 1);
        assert!(
            hits.iter().any(|h| (h.strand, h.t_start, h.t_end) == span),
            "expected {span:?}, got {hits:?}"
        );
    }

    #[cfg(not(miri))]
    mod generated {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn random_configs_do_not_panic(
                query in prop::collection::vec(0usize..4, 4..17),
                target in prop::collection::vec(0usize..5, 4..25),
                offset in 0usize..32,
                requested_length in 2usize..=12,
                mismatch in 0usize..=2,
                prefix in 0usize..=3,
                suffix in 0usize..=3,
                wobble in any::<bool>(),
                prune in any::<bool>(),
                interval in any::<bool>(),
                reverse in any::<bool>(),
                model in 0usize..4,
                temp in 0usize..9,
                penalty in prop_oneof![Just(0i32), Just(1), Just(500000), 0i32..50000],
                window in prop_oneof![Just(-1i32), Just(0), Just(1), 2i32..20],
                threshold in prop_oneof![Just(i32::MAX), Just(i32::MIN), -100000i32..100000],
            ) {
                let alphabet = b"ACGUN";
                let q: Vec<_> = query.iter().map(|&b| alphabet[b]).collect();
                let length = requested_length.min(q.len()).min(target.len());
                let start = offset % (q.len() - length + 1);
                let mut config = SearchConfig::default();
                config.seed.seed_length = Some(length as i64);
                config.seed.max_mismatches = mismatch;
                config.seed.min_prefix_matches = prefix;
                config.seed.min_suffix_matches = suffix;
                config.seed.seed_wobble = wobble;
                config.seed.no_max_prune = !prune;
                if interval {
                    config.seed.seed_start = Some(start as i64 + 1);
                    config.seed.seed_end = Some((start + length) as i64);
                }
                let names = ["t04", "slh04", "s95-rna-dna", "s95-dna-rna"];
                config.score.dsm_id = DsmId::from(names[model]);
                config.score.temperature = [0, 12, 25, 31, 37, 40, 42, 46, 50][temp];
                config.score.penalty = Energy(penalty);
                config.extend.max_extension = window;
                config.filter.delta_g = Energy(threshold);
                let mut planted: Vec<_> = target.iter().map(|&b| alphabet[b]).collect();
                let site = offset % (planted.len() - length + 1);
                let complement = |b| Base::try_from(char::from(b)).unwrap().complement().to_byte().to_ascii_uppercase();
                let mut duplex: Vec<_> = q[start..start + length].iter().copied().map(complement).collect();
                // Forward FASTA is reverse(physical); reverse FASTA is complement(physical).
                if reverse { duplex.iter_mut().for_each(|b| *b = complement(*b)); }
                else { duplex.reverse(); }
                planted[site..site+length].copy_from_slice(&duplex);

                // No-panic: three non-planted variants must survive every config.
                for (qcase, tcase) in [
                    (q.clone(), target.iter().map(|&b| alphabet[b]).collect::<Vec<_>>()),
                    (vec![b'G'; q.len()], vec![b'C'; target.len()]),
                    (vec![b'N'; q.len()], vec![b'N'; target.len()]),
                ] {
                    search_pair(
                        std::str::from_utf8(&qcase).unwrap(),
                        std::str::from_utf8(&tcase).unwrap(),
                        &config,
                    )
                    .unwrap();
                }

                // Planted case: always runs (no-panic), and under permissive
                // configs asserts the planted site produces a covering hit.
                let planted_hits = search_pair(
                    std::str::from_utf8(&q).unwrap(),
                    std::str::from_utf8(&planted).unwrap(),
                    &config,
                )
                .unwrap();

                let permissive = threshold == i32::MAX
                    && mismatch == 0
                    && prefix <= 1
                    && suffix == 0
                    && !interval
                    && !prune;

                if permissive {
                    let expected_strand = if reverse { Strand::Reverse } else { Strand::Forward };
                    let t_end_inclusive = site + length - 1;
                    let covers_site = planted_hits.iter().any(|h| {
                        h.strand == expected_strand
                            && h.t_start <= site
                            && h.t_end >= t_end_inclusive
                    });
                    prop_assert!(
                        covers_site,
                        "planted site at FASTA {}..={} ({:?}) not covered by {} hits: {:?}",
                        site, t_end_inclusive, expected_strand,
                        planted_hits.len(),
                        planted_hits.iter()
                            .map(|h| (h.t_start, h.t_end, char::from(h.strand)))
                            .collect::<Vec<_>>(),
                    );
                }
            }
        }
    }
}
