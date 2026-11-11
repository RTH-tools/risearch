//! Shared inputs, configuration and driver for the search tests.

use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use super::{run_search, SearchHit, VecSink};
use crate::error::Result;
use crate::fastx::read_sequences;
use crate::{SearchConfig, Sequence, TargetRegistry};

/// A `tests/data` file, resolved against the crate root the way
/// `tests/cli_library_agreement.rs` resolves it.
fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// Query records of the shared corpus.
pub(super) static QUERIES: LazyLock<Vec<(String, Sequence)>> =
    LazyLock::new(|| read_sequences(data("query.fa")).unwrap());

/// Target records of the shared corpus.
pub(super) static TARGET_RECORDS: LazyLock<Vec<(String, Sequence)>> =
    LazyLock::new(|| read_sequences(data("target.fa")).unwrap());

/// Index built over [`TARGET_RECORDS`].
pub(super) static TARGETS: LazyLock<TargetRegistry> =
    LazyLock::new(|| TargetRegistry::build(TARGET_RECORDS.clone(), None).unwrap());

/// Records from `(name, sequence text)` pairs.
pub(super) fn records(pairs: &[(&str, &str)]) -> Vec<(String, Sequence)> {
    pairs
        .iter()
        .map(|&(name, text)| {
            (
                name.to_string(),
                Sequence::try_from(text.as_bytes()).unwrap(),
            )
        })
        .collect()
}

/// The default configuration with `edit` applied.
pub(super) fn config(edit: impl FnOnce(&mut SearchConfig)) -> SearchConfig {
    let mut config = SearchConfig::default();
    edit(&mut config);
    config
}

/// Run the production search and collect its hits.
pub(super) fn search(
    queries: &[(String, Sequence)],
    targets: &TargetRegistry,
    config: &SearchConfig,
) -> Result<Vec<SearchHit>> {
    let sink = VecSink::default();
    run_search(queries, targets, config, &sink)?;
    Ok(sink.into_hits())
}
