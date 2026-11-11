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
    targets: &'a TargetRegistry,
    format: OutputFormat,
    out: OutputWriter,
}

impl<'a> TextSink<'a> {
    /// `output_path` is a file, or a directory when `output.multifile` is set.
    ///
    /// `queries` and `targets` must be the ones later passed to
    /// [`run_search`](crate::run_search): `consume` indexes them by the driver's
    /// `query_idx`, so a different pair silently attributes rows to the wrong
    /// sequences.
    pub fn new(
        queries: &'a [(String, Sequence)],
        targets: &'a TargetRegistry,
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
            targets,
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
        let tview = self.targets.view();

        // `format_hit_into` reserves per hit, so the buffer sizes itself.
        let mut block = Vec::new();
        for hit in &hits {
            let target_idx = hit.target_index();
            let target = tview.target(target_idx, hit.strand);
            let target_range = hit.duplex_target_range(tview);
            let t_name = self.targets.get_name(target_idx);
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

#[cfg(test)]
mod tests {
    use std::io::Read;

    use rstest::rstest;

    use super::*;
    use crate::config::OutputCompression;
    use crate::fastx::read_sequences_from;
    use crate::types::Strand;
    use crate::Energy;

    fn queries() -> Vec<(String, Sequence)> {
        read_sequences_from(b">q1\nAAAA\n>q2\nCCCC\n>q3\nAAAAA\n".as_slice(), "inline").unwrap()
    }

    fn targets() -> TargetRegistry {
        let records = read_sequences_from(b">t\nUUUUUUUU\n".as_slice(), "inline").unwrap();
        TargetRegistry::build(records, None).unwrap()
    }

    /// A four-base hit of query `query_idx` on the only target, from `t_start`.
    fn hit(query_idx: u32, t_start: usize) -> SearchHit {
        SearchHit {
            query_idx,
            target_idx: 0,
            q_start: 0,
            q_end: 3,
            t_start,
            t_end: t_start + 3,
            strand: Strand::Forward,
            energy: Energy::from_kcal(-1.5),
            alignment: None,
        }
    }

    fn minimal_output() -> OutputConfig {
        OutputConfig {
            format: OutputFormat::Minimal,
            compress: OutputCompression::None,
            multifile: false,
        }
    }

    fn multifile_output() -> OutputConfig {
        OutputConfig {
            multifile: true,
            ..minimal_output()
        }
    }

    /// Hand `batches` to a fresh sink at `path` and flush it. The sink dies here,
    /// so the caller reads a file whose compression trailer is already written.
    fn write(output: &OutputConfig, path: &Path, batches: &[(usize, Vec<SearchHit>)]) {
        let (queries, targets) = (queries(), targets());
        let sink = TextSink::new(&queries, &targets, output, path).unwrap();
        for (query_idx, hits) in batches {
            sink.consume(*query_idx, hits.clone()).unwrap();
        }
        sink.flush().unwrap();
    }

    fn lines(path: &Path) -> Vec<String> {
        fs_err::read_to_string(path)
            .unwrap()
            .lines()
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn each_hit_becomes_one_row_under_its_query_and_target_names() {
        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        write(
            &minimal_output(),
            out.path(),
            &[(0, vec![hit(0, 0), hit(0, 4)]), (2, vec![hit(2, 2)])],
        );

        assert_eq!(
            lines(out.path()),
            [
                "q1\t1\t4\tt\t1\t4\t+\t-1.50",
                "q1\t1\t4\tt\t5\t8\t+\t-1.50",
                "q3\t1\t4\tt\t3\t6\t+\t-1.50",
            ]
        );
    }

    /// Zero-hit behaviour differs by topology: single-file truncates the output
    /// to nothing, multifile creates no file at all.
    #[test]
    fn flush_without_hits_empties_a_single_file_and_leaves_no_multifile_entry() {
        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        fs_err::write(out.path(), b"stale").unwrap();
        write(&minimal_output(), out.path(), &[(0, Vec::new())]);
        assert_eq!(fs_err::metadata(out.path()).unwrap().len(), 0);

        let tmp = tempfile::tempdir().unwrap();
        // Not pre-created, so this also pins that the directory itself is made.
        let dir = tmp.path().join("multi");
        write(&multifile_output(), &dir, &[(0, Vec::new())]);
        assert_eq!(fs_err::read_dir(&dir).unwrap().count(), 0);
    }

    #[rstest]
    #[case::with_hits(vec![(0, vec![hit(0, 0), hit(0, 4)]), (2, vec![hit(2, 2)])])]
    #[case::without_hits(Vec::new())]
    fn compressed_output_decodes_to_the_plain_text(
        #[case] batches: Vec<(usize, Vec<SearchHit>)>,
        #[values(OutputCompression::Gzip(6), OutputCompression::Zstd(3))]
        compress: OutputCompression,
    ) {
        let plain = tempfile::NamedTempFile::new().unwrap();
        write(&minimal_output(), plain.path(), &batches);

        let packed = tempfile::NamedTempFile::new().unwrap();
        let output = OutputConfig {
            compress,
            ..minimal_output()
        };
        write(&output, packed.path(), &batches);
        let bytes = fs_err::read(packed.path()).unwrap();
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
                String::from_utf8(zstd::stream::decode_all(bytes.as_slice()).unwrap()).unwrap()
            }
            OutputCompression::None => unreachable!(),
        };

        assert_eq!(decoded, fs_err::read_to_string(plain.path()).unwrap());
    }

    /// The driver builds its sink before it validates the config, so building one
    /// must not touch the destination.
    #[test]
    fn building_a_sink_touches_nothing() {
        let (queries, targets) = (queries(), targets());
        let tmp = tempfile::tempdir().unwrap();

        let file = tmp.path().join("hits.tsv");
        fs_err::write(&file, b"stale").unwrap();
        TextSink::new(&queries, &targets, &minimal_output(), &file).unwrap();
        assert_eq!(fs_err::read(&file).unwrap(), b"stale");

        let dir = tmp.path().join("multi");
        TextSink::new(&queries, &targets, &multifile_output(), &dir).unwrap();
        assert!(!dir.exists());
    }

    /// Multifile addresses its writers by query index, queries without hits get no
    /// file, and queries finish in any order, so a query dropping out in the middle
    /// must not shift the others' rows into the wrong file.
    #[test]
    fn multifile_writes_each_querys_rows_to_its_own_file() {
        let files = |batches: &[(usize, Vec<SearchHit>)]| {
            let tmp = tempfile::tempdir().unwrap();
            write(&multifile_output(), tmp.path(), batches);
            let mut files: Vec<(String, Vec<String>)> = fs_err::read_dir(tmp.path())
                .unwrap()
                .map(|entry| {
                    let path = entry.unwrap().path();
                    let stem = path.file_stem().unwrap().to_string_lossy().into_owned();
                    (stem, lines(&path))
                })
                .collect();
            files.sort();
            files
        };
        let q1 = (0, vec![hit(0, 0), hit(0, 4)]);
        let q2 = (1, Vec::new());
        let q3 = (2, vec![hit(2, 2)]);

        let in_order = files(&[q1.clone(), q2.clone(), q3.clone()]);
        assert_eq!(
            in_order,
            files(&[q3, q2, q1]),
            "files depend on the order queries finish in"
        );
        assert_eq!(
            in_order,
            [
                (
                    "q1".to_string(),
                    vec![
                        "q1\t1\t4\tt\t1\t4\t+\t-1.50".to_string(),
                        "q1\t1\t4\tt\t5\t8\t+\t-1.50".to_string(),
                    ]
                ),
                (
                    "q3".to_string(),
                    vec!["q3\t1\t4\tt\t3\t6\t+\t-1.50".to_string()]
                ),
            ]
        );
    }
}
