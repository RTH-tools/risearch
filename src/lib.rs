//! RNA-RNA interaction search: suffix-array seeding with energy-based extension.
//!
//! The library behind the `risearch` binary. A run has two stages — build a
//! target index once, then search queries against it.
//!
//! # Indexing
//!
//! ```no_run
//! use std::path::Path;
//!
//! use risearch::fastx::read_sequences;
//! use risearch::TargetRegistry;
//!
//! let targets = read_sequences(Path::new("targets.fa"))?;
//! TargetRegistry::build(targets, None)?.save(Path::new("targets.idx"))?;
//! # Ok::<(), risearch::Error>(())
//! ```
//!
//! # Searching
//!
//! ```no_run
//! use std::path::Path;
//!
//! use risearch::fastx::read_sequences;
//! use risearch::{run_search, SearchConfig, TargetRegistry, VecSink};
//!
//! let opts = SearchConfig::default();
//! let queries = read_sequences(Path::new("queries.fa"))?;
//! let targets = TargetRegistry::open(Path::new("targets.idx"))?;
//!
//! let sink = VecSink::default();
//! run_search(&queries, &targets, &opts, &sink)?;
//!
//! for hit in sink.into_hits() {
//!     println!("{}\t{}", targets.get_name(hit.target_index()), hit.energy);
//! }
//! # Ok::<(), risearch::Error>(())
//! ```
//!
//! [`VecSink`] keeps every hit in memory. Implement [`HitSink`] to stream them
//! instead; that is what [`TextSink`] and the Python bindings' Arrow sink do.
//! Output order is not stable across runs — the hit set is. Sort if you need one.
//!
//! [`SearchConfig`] and [`OutputConfig`] are built directly.
//! [`SearchConfig::validate`] checks every part that has invalid states, and
//! [`SearchConfig::default`] supplies the same defaults the CLI applies.
//!
//! # Threads
//!
//! Index construction is parallel only through OpenMP (the `openmp` feature);
//! its thread count is the `threads` argument of [`TargetRegistry::build`].
//! Everything else parallel — query preparation, index validation, search — runs
//! on the caller's rayon pool: the global pool, or the one entered with
//! [`rayon::ThreadPool::install`]. The library never creates or configures a
//! pool.

#[doc(hidden)]
pub mod adapter;
pub mod alignment;
pub mod config;
#[doc(hidden)]
pub mod dp;
#[doc(hidden)]
pub mod dsm;
pub mod error;
pub mod fastx;
pub(crate) mod index;
pub mod output;
#[doc(hidden)]
pub mod registry;
pub mod search;
#[doc(hidden)]
pub mod seed;
pub mod seq;
pub mod types; // Core domain types

// Re-exports for convenience
pub use alignment::{AlignColumn, PairClass};
pub use error::{Error, Result};
pub use index::store::TargetRegistry;
pub use index::TargetView;
pub use seq::Sequence;
pub use types::{Base, Energy, Strand};

// Search API re-exports for library usage
pub use config::{
    ExtendConfig, FilterConfig, OutputCompression, OutputConfig, OutputFormat, ScoreConfig,
    SearchConfig, SeedConfig,
};
pub use output::TextSink;
pub use search::{run_search, HitSink, SearchHit, VecSink};
pub use types::DsmId;

#[cfg(kani)]
mod proofs;
