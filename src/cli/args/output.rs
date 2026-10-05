use risearch::config::{OutputCompression, OutputConfig, OutputFormat, STDIO};

use anyhow::{bail, ensure, Context, Error, Result};
use clap::ValueEnum;
use flate2::Compression;
use std::path::{Path, PathBuf};

/// CLI-facing codec selector — what the `--compress` flag parses into.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
#[clap(rename_all = "lowercase")]
pub(crate) enum OutputCodec {
    None,
    #[value(alias = "gz")]
    Gzip,
    #[value(alias = "zst")]
    Zstd,
}

impl From<&std::path::Path> for OutputCodec {
    /// Infer codec from file extension. Unrecognised or absent → `None`.
    fn from(path: &std::path::Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some(ext) => match ext.to_ascii_lowercase().as_str() {
                "gz" | "gzip" => Self::Gzip,
                "zst" | "zstd" => Self::Zstd,
                _ => Self::None,
            },
            None => Self::None,
        }
    }
}

/// Reject an output path whose parent directory is missing or is not a directory.
///
/// A CLI-boundary check: it exists so a bad `-o` fails before the work that
/// would fill it, not to guard the write itself.
pub(crate) fn validate_output_parent(path: &Path) -> Result<()> {
    let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) else {
        return Ok(());
    };

    let md = fs_err::metadata(parent)
        .with_context(|| format!("output directory '{}' does not exist", parent.display()))?;
    if !md.is_dir() {
        bail!("output path '{}' is not a directory", parent.display());
    }
    Ok(())
}

/// Boundary CLI arguments for output destination, formatting, and compression.
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct OutputArgs {
    /// Output file, '-' for stdout; a .gz or .zst file name compresses it
    ///
    /// With --multifile, the directory that receives one file per query.
    #[arg(short = 'o', long = "output", value_name = "FILE", default_value = STDIO)]
    pub(crate) path: PathBuf,

    /// Layout of each hit
    #[arg(
        short = 'f',
        long = "format",
        value_name = "FORMAT",
        value_enum,
        default_value_t = OutputFormat::default()
    )]
    pub(crate) report_format: OutputFormat,

    /// Compression codec, overriding the one chosen from the -o file name
    ///
    /// gz and zst are accepted as short forms.
    #[arg(long = "compress", value_name = "CODEC", value_enum)]
    pub(crate) output_compress: Option<OutputCodec>,

    #[arg(
        long = "compress-level",
        value_name = "LEVEL",
        allow_hyphen_values = true,
        hide_short_help = true,
        help = format!(
            "Output compression level: gzip {}–{}; zstd up to {}, negative levels trade ratio for speed [default: {} for gzip, {} for zstd]",
            Compression::none().level(),
            Compression::best().level(),
            zstd::compression_level_range().end(),
            Compression::default().level(),
            zstd::DEFAULT_COMPRESSION_LEVEL
        )
    )]
    pub(crate) output_level: Option<i32>,

    /// Write one file per query into the directory given by -o
    ///
    /// Each file is named after its query and ends in .tsv, .tsv.gz or .tsv.zst. Characters that are
    /// unsafe in file names become '_', repeated names get _1, _2 and so on, and the directory is
    /// created if missing.
    #[arg(long = "multifile", action = clap::ArgAction::SetTrue)]
    pub(crate) output_multifile: bool,
}

impl TryFrom<OutputArgs> for OutputConfig {
    type Error = Error;

    fn try_from(value: OutputArgs) -> Result<Self, Self::Error> {
        if value.output_multifile && value.path == Path::new(STDIO) {
            bail!("--multifile requires -o/--output to be a directory path; '-' (stdout) is not allowed.");
        }

        if value.path != Path::new(STDIO) {
            validate_output_parent(&value.path)?;
        }

        let codec = value
            .output_compress
            .unwrap_or_else(|| OutputCodec::from(value.path.as_path()));

        let compress = match (codec, value.output_level) {
            (OutputCodec::None, None) => OutputCompression::None,
            (OutputCodec::None, Some(_)) => bail!("--compress-level requires compressed output"),
            (OutputCodec::Gzip, level) => {
                let (none, best) = (
                    Compression::none().level() as i32,
                    Compression::best().level() as i32,
                );
                let level = level.unwrap_or(Compression::default().level() as i32);
                ensure!(
                    (none..=best).contains(&level),
                    "gzip level must be {none}-{best} (got {level})"
                );
                OutputCompression::Gzip(level as u8)
            }
            (OutputCodec::Zstd, level) => {
                let range = zstd::compression_level_range();
                let level = level.unwrap_or(zstd::DEFAULT_COMPRESSION_LEVEL);
                ensure!(
                    range.contains(&level),
                    "zstd level must be {}..{} (got {level})",
                    range.start(),
                    range.end()
                );
                OutputCompression::Zstd(level)
            }
        };

        Ok(OutputConfig {
            format: value.report_format,
            compress,
            multifile: value.output_multifile,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(
        path: &str,
        output_compress: Option<OutputCodec>,
        output_multifile: bool,
    ) -> OutputArgs {
        OutputArgs {
            path: path.into(),
            report_format: OutputFormat::default(),
            output_compress,
            output_level: None,
            output_multifile,
        }
    }

    #[test]
    fn codec_is_inferred_from_extension_unless_given_explicitly() {
        for (path, codec, expected) in [
            ("out.gz", None, OutputCompression::Gzip(6)),
            ("out.zst", None, OutputCompression::Zstd(3)),
            ("out.tsv", None, OutputCompression::None),
            (
                "out.tsv",
                Some(OutputCodec::Gzip),
                OutputCompression::Gzip(6),
            ),
            (
                "out.tsv",
                Some(OutputCodec::Zstd),
                OutputCompression::Zstd(3),
            ),
        ] {
            let config = OutputConfig::try_from(args(path, codec, false)).unwrap();
            assert_eq!(config.compress, expected, "{path} {codec:?}");
        }
    }

    #[test]
    fn multifile_rejects_stdout_output() {
        let err = OutputConfig::try_from(args("-", None, true)).unwrap_err();
        assert!(err
            .to_string()
            .contains("--multifile requires -o/--output to be a directory path"));
    }

    #[test]
    fn compress_level_is_range_checked_per_codec() {
        for (codec, level, expected) in [
            (OutputCodec::Gzip, 0, Some(OutputCompression::Gzip(0))),
            (OutputCodec::Gzip, 9, Some(OutputCompression::Gzip(9))),
            (OutputCodec::Gzip, 10, None),
            (OutputCodec::Gzip, -1, None),
            (
                OutputCodec::Zstd,
                -131072,
                Some(OutputCompression::Zstd(-131072)),
            ),
            (OutputCodec::Zstd, 22, Some(OutputCompression::Zstd(22))),
            (OutputCodec::Zstd, -131073, None),
            (OutputCodec::Zstd, 23, None),
            (OutputCodec::None, 1, None),
        ] {
            let mut cfg = args("out.tsv", Some(codec), false);
            cfg.output_level = Some(level);
            match expected {
                Some(compress) => assert_eq!(
                    OutputConfig::try_from(cfg).unwrap().compress,
                    compress,
                    "{codec:?} {level}"
                ),
                None => assert!(OutputConfig::try_from(cfg).is_err(), "{codec:?} {level}"),
            }
        }
    }

    #[test]
    fn output_parent_directory_must_exist() {
        let err = OutputConfig::try_from(args("no_such_dir/out.tsv", None, false)).unwrap_err();
        assert!(err.to_string().contains("does not exist"), "{err}");
    }
}
