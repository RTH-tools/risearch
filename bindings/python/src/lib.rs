mod arrow;
mod error;

use std::num::NonZeroUsize;
use std::path::PathBuf;

use arrow_array::ffi_stream::FFI_ArrowArrayStream;
use arrow_array::{RecordBatch, RecordBatchIterator};
use arrow_schema::SchemaRef;
use pyo3::prelude::*;
use pyo3::types::{PyCapsule, PyDict};
use risearch::fastx::read_sequences;
use risearch::{
    run_search, DsmId, Energy, ExtendConfig, FilterConfig, ScoreConfig, SearchConfig, SeedConfig,
    TargetRegistry,
};

use crate::arrow::{search_result_schema, ArrowSink};
use crate::error::Result;

// =============================================================================
// PySearchResult — Arrow C Stream Interface producer
// =============================================================================

/// The result of a `search()` call.
///
/// Consumed by `pl.DataFrame(result)` via the Arrow PyCapsule Interface.
/// One-shot: the stream is consumed on the first call to `__arrow_c_stream__`.
#[pyclass(name = "SearchResult", module = "risearch._native")]
struct PySearchResult {
    batch: Option<RecordBatch>,
    schema: SchemaRef,
}

impl PySearchResult {
    fn new(batch: RecordBatch, schema: SchemaRef) -> Self {
        Self {
            batch: Some(batch),
            schema,
        }
    }
}

#[pymethods]
impl PySearchResult {
    /// Arrow PyCapsule Interface producer (`__arrow_c_stream__` protocol).
    ///
    /// Called automatically by `pl.DataFrame(result)` — do not call directly.
    // `-> "object"`: the inferred `types.CapsuleType` needs Python 3.13.
    #[pyo3(signature = (requested_schema = None) -> "object")]
    fn __arrow_c_stream__<'py>(
        &mut self,
        py: Python<'py>,
        requested_schema: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let _ = requested_schema;
        let batch = self
            .batch
            .take()
            .ok_or_else(|| pyo3::exceptions::PyIOError::new_err("Arrow stream already consumed"))?;

        let reader = RecordBatchIterator::new(std::iter::once(Ok(batch)), self.schema.clone());
        let stream = FFI_ArrowArrayStream::new(Box::new(reader));

        PyCapsule::new_with_value(py, stream, c"arrow_array_stream")
    }
}

// =============================================================================
// PyTargetRegistry
// =============================================================================

/// An mmap-backed target index for RNA-RNA interaction search.
///
/// Build once with `build_index()`, then reuse across many `search()` calls.
#[pyclass(name = "TargetRegistry", module = "risearch")]
pub struct PyTargetRegistry(TargetRegistry);

#[pymethods]
impl PyTargetRegistry {
    /// Load a pre-built index from disk.
    ///
    /// Parameters
    /// ----------
    /// path : str | os.PathLike
    ///     Path to the `.idx` file produced by `build_index()`.
    #[staticmethod]
    fn open(py: Python<'_>, path: PathBuf) -> Result<Self> {
        // Archive validation is rayon-parallel: a worker that logs deadlocks if
        // this thread holds the GIL.
        let store = py.detach(|| TargetRegistry::open(&path))?;
        Ok(PyTargetRegistry(store))
    }

    fn __repr__(&self) -> String {
        format!("TargetRegistry(targets={})", self.0.len())
    }
}

// =============================================================================
// Module-level functions
// =============================================================================

/// Run `f` on `pool`, or on rayon's global pool when no width was pinned.
fn in_pool<R: Send>(pool: Option<&rayon::ThreadPool>, f: impl FnOnce() -> R + Send) -> R {
    match pool {
        Some(p) => p.install(f),
        None => f(),
    }
}

/// Rejects 0 and sklearn's `-1` rather than aliasing them: `None` already
/// means every core, and rayon would silently read 0 as "auto".
fn thread_count(threads: Option<i64>) -> PyResult<Option<NonZeroUsize>> {
    match threads {
        None => Ok(None),
        Some(n) if n >= 1 => Ok(NonZeroUsize::new(n as usize)),
        Some(n) => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "threads must be >= 1, or None for every core; got {n}"
        ))),
    }
}

/// Build a binary index from a FASTA file of target sequences.
///
/// Parameters
/// ----------
/// fasta : str | os.PathLike
///     Input FASTA file containing target sequences.
/// output : str | os.PathLike
///     Destination path for the binary index (`.idx`).
/// threads : int, optional
///     OpenMP threads for suffix-array construction; `None` uses the OpenMP
///     default (OMP_NUM_THREADS, else every core).
///     Only takes effect when risearch is built with the `openmp` feature,
///     which the published wheels are not; without it construction is
///     single-threaded and this is ignored.
///
/// Notes
/// -----
/// The GIL is released during index construction so other Python threads
/// may run concurrently.
#[pyfunction]
#[pyo3(signature = (fasta, output, *, threads = None))]
fn build_index(
    py: Python<'_>,
    fasta: PathBuf,
    output: PathBuf,
    threads: Option<i64>,
) -> Result<()> {
    let threads = thread_count(threads)?;
    py.detach(|| {
        let targets = read_sequences(&fasta)?;
        TargetRegistry::build(targets, threads)?.save(&output)
    })?;
    Ok(())
}

/// Search for RNA-RNA interactions between query sequences and an indexed target set.
///
/// Parameters
/// ----------
/// query : list[str | os.PathLike]
///     One or more FASTA files containing query sequences.
/// target : TargetRegistry
///     Pre-loaded target index (from `TargetRegistry.open()`).
/// alignment : bool
///     Populate the `alignment` column. Set False to skip DP traceback when
///     only coordinates and energies are needed; the column becomes all-null.
/// threads : int, optional
///     Worker threads for the search. `None` leaves the choice to rayon,
///     which honours RAYON_NUM_THREADS.
///
/// Returns
/// -------
/// SearchResult
///     Arrow C Stream producer; pass to `pl.DataFrame()` to obtain a DataFrame.
///
/// Notes
/// -----
/// The GIL is released during the search so other Python threads may run.
#[pyfunction]
// No defaults here: the Python wrapper owns them.
#[pyo3(signature = (
    query,
    target,
    *,
    seed_length,
    seed_start,
    seed_end,
    energy_threshold,
    mismatches,
    mismatch_prefix,
    mismatch_suffix,
    seed_wobble,
    matrix,
    penalty,
    temperature,
    max_extension,
    seed_energy,
    no_max_prune,
    no_dedup,
    alignment,
    threads,
))]
// Binding surface: each argument maps to a documented Python keyword parameter,
// so the flat signature is the public API and grouping would break it.
#[allow(clippy::too_many_arguments)]
fn search(
    py: Python<'_>,
    query: Vec<PathBuf>,
    target: &PyTargetRegistry,
    seed_length: Option<i64>,
    seed_start: Option<i64>,
    seed_end: Option<i64>,
    energy_threshold: f64,
    mismatches: usize,
    mismatch_prefix: usize,
    mismatch_suffix: usize,
    seed_wobble: bool,
    matrix: PathBuf,
    penalty: f64,
    temperature: Option<i32>,
    max_extension: i32,
    seed_energy: f64,
    no_max_prune: bool,
    no_dedup: bool,
    alignment: bool,
    threads: Option<i64>,
) -> Result<PySearchResult> {
    let config = SearchConfig {
        seed: SeedConfig {
            seed_start,
            seed_end,
            seed_length,
            seed_wobble,
            no_max_prune,
            max_mismatches: mismatches,
            min_prefix_matches: mismatch_prefix,
            min_suffix_matches: mismatch_suffix,
        },
        score: ScoreConfig::new(
            DsmId::from(matrix.to_string_lossy().as_ref()),
            Energy::try_from(penalty).map_err(pyo3::exceptions::PyValueError::new_err)?,
            temperature,
        ),
        extend: ExtendConfig {
            max_extension,
            build_alignment: alignment,
        },
        filter: FilterConfig {
            delta_g: Energy::try_from(energy_threshold)
                .map_err(pyo3::exceptions::PyValueError::new_err)?,
            seed_energy: Energy::try_from(seed_energy)
                .map_err(pyo3::exceptions::PyValueError::new_err)?,
            no_dedup,
        },
    };
    config.validate()?;

    let threads = thread_count(threads)?;
    // Scoped, not build_global: build_global errors on a second call. Detached
    // because building it spawns OS threads.
    let pool = py
        .detach(|| {
            threads
                .map(|n| rayon::ThreadPoolBuilder::new().num_threads(n.get()).build())
                .transpose()
        })
        .map_err(|err| {
            pyo3::exceptions::PyRuntimeError::new_err(format!(
                "failed to build the thread pool: {err}"
            ))
        })?;

    let queries = py.detach(|| -> Result<Vec<_>> {
        let mut records = Vec::new();
        for path in &query {
            records.extend(read_sequences(path)?);
        }
        Ok(records)
    })?;

    let sink = ArrowSink::new(&queries, &target.0);
    // Detached: a rayon worker that logs deadlocks against a held GIL.
    py.detach(|| {
        in_pool(pool.as_ref(), || {
            run_search(&queries, &target.0, &config, &sink)
        })
    })?;
    let schema = search_result_schema().clone();
    Ok(PySearchResult::new(sink.into_batch(&schema), schema))
}

// =============================================================================
// Module registration
// =============================================================================

/// The Rust config defaults, keyed by their Python keyword name.
///
/// `alignment` is absent on purpose: `OutputFormat::needs_alignment` owns
/// `build_alignment`, so `ExtendConfig::default()` is not its canonical value.
/// `temperature` is `None`: the wrapper passes it through unresolved so an
/// explicit value can be detected, and `ScoreConfig::new` applies the default.
#[pyfunction]
fn _default_options(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    let seed = SeedConfig::default();
    let score = ScoreConfig::default();
    let extend = ExtendConfig::default();
    let filter = FilterConfig::default();

    let d = PyDict::new(py);
    d.set_item("seed_length", seed.seed_length)?;
    d.set_item("seed_start", seed.seed_start)?;
    d.set_item("seed_end", seed.seed_end)?;
    d.set_item("energy_threshold", f64::from(filter.delta_g))?;
    d.set_item("mismatches", seed.max_mismatches)?;
    d.set_item("mismatch_prefix", seed.min_prefix_matches)?;
    d.set_item("mismatch_suffix", seed.min_suffix_matches)?;
    d.set_item("seed_wobble", seed.seed_wobble)?;
    d.set_item("matrix", score.dsm_id.0)?;
    d.set_item("penalty", f64::from(score.penalty))?;
    d.set_item("temperature", None::<i32>)?;
    d.set_item("max_extension", extend.max_extension)?;
    d.set_item("seed_energy", f64::from(filter.seed_energy))?;
    d.set_item("no_max_prune", seed.no_max_prune)?;
    d.set_item("no_dedup", filter.no_dedup)?;
    Ok(d)
}

// Inline, not a function: stub generation only introspects inline modules.
#[pymodule]
mod _native {
    use pyo3::prelude::*;

    #[pymodule_export]
    use super::{_default_options, build_index, search, PySearchResult, PyTargetRegistry};

    #[pymodule_init]
    fn init(_m: &Bound<'_, PyModule>) -> PyResult<()> {
        // try_init, not init: init panics when another extension already installed
        // a logger, and this body re-runs on module reload.
        let _ = pyo3_log::try_init();
        Ok(())
    }
}
