//! Renders search hits as text and hands the bytes to the output writer.

use std::path::Path;

use crate::error::Result;

use crate::config::{OutputConfig, OutputFormat};
use crate::index::store::TargetRegistry;
use crate::output::format::format_hit_into;
use crate::output::writer::{build_multifile_paths, OutputWriter};
use crate::search::{HitSink, SearchHit};
use crate::seq::Sequence;

/// Turns each query's hits into text rows.
///
/// Knows nothing about files: it renders a query's whole block and hands it to the
/// writer keyed by query index, leaving destinations, paths and open/flush timing
/// to `OutputWriter`. Because the writer opens lazily, constructing a `TextSink`
/// touches nothing a rejected config would need left intact.
///
/// gzip and zstd write their trailer when the writer drops, so drop the sink only
/// after the search returns.
pub struct TextSink<'a> {
    queries: &'a [(String, Sequence)],
    store: &'a TargetRegistry,
    format: OutputFormat,
    out: OutputWriter,
}

impl<'a> TextSink<'a> {
    /// `output_path` is a file, or a directory when `output.multifile` is set.
    ///
    /// `queries` and `store` must be the ones later passed to
    /// [`run_search`](crate::run_search): `consume` indexes them by the driver's
    /// `query_idx`, so a different pair silently attributes rows to the wrong
    /// sequences.
    pub fn new(
        queries: &'a [(String, Sequence)],
        store: &'a TargetRegistry,
        output: &OutputConfig,
        output_path: &Path,
    ) -> Result<Self> {
        let out = if output.multifile {
            // Named up front so a `_1` collision suffix never depends on the order
            // queries finish in.
            let paths = build_multifile_paths(
                queries.iter().map(|(name, _)| name.as_str()),
                output_path,
                output.compress.extension(),
            );
            OutputWriter::per_key(output_path, paths, output.compress)
        } else {
            OutputWriter::single(output_path, output.compress)
        };

        Ok(Self {
            queries,
            store,
            format: output.format,
            out,
        })
    }
}

impl HitSink for TextSink<'_> {
    fn consume(&self, query_idx: usize, hits: Vec<SearchHit>) -> Result<()> {
        if hits.is_empty() {
            return Ok(());
        }
        let q_name = &self.queries[query_idx].0;

        // Loop-invariant: rebuilding it per hit re-enters the rkyv root.
        let tview = self.store.view();

        // `format_hit_into` reserves per hit, so the buffer sizes itself.
        let mut block = Vec::new();
        for hit in &hits {
            let target_idx = hit.target_index();
            let target = tview.target(target_idx, hit.strand);
            let target_range = hit.duplex_target_range(tview);
            let t_name = self.store.get_name(target_idx);
            format_hit_into(
                &mut block,
                hit,
                q_name,
                t_name,
                target,
                target_range,
                self.format,
            );
        }

        // Rendering above is off-lock; the writer sees one bulk write per query.
        self.out.write_block(query_idx, &block)
    }

    fn flush(&self) -> Result<()> {
        self.out.finish()
    }
}
