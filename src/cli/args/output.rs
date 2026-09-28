use risearch::config::{OutputCompression, OutputConfig, OutputFormat};

use anyhow::{bail, Context, Error, Result};
use clap::ValueEnum;
use log::warn;
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
    /// Output file for search results (use '-' for stdout)
    #[arg(short = 'o', long = "output", value_name = "FILE", default_value = "-")]
    pub(crate) path: PathBuf,

    /// Output format
    #[arg(
        short = 'f',
        long = "format",
        value_name = "FORMAT",
        default_missing_value = "minimal",
        value_enum
    )]
    pub(crate) report_format: Option<OutputFormat>,

    /// DEPRECATED: Legacy argument for output format (1=detailed, 2=cigar, 3=binding_site, 4=minimal)
    #[arg(
        short = 'p',
        long = "report-alignment",
        value_name = "MODE",
        num_args = 0..=1,
        default_missing_value = "1",
        value_parser = parse_legacy_format,
        help_heading = "Deprecated"
    )]
    pub(crate) report_legacy: Option<OutputFormat>,

    /// Output compression codec (overrides file extension inference; gzip/gz, zstd/zst accepted)
    #[arg(long = "compress", value_enum)]
    pub(crate) output_compress: Option<OutputCodec>,

    /// Output compression level (codec-specific: gzip 0–9, zstd -7..22)
    #[arg(long = "compress-level", value_name = "LEVEL")]
    pub(crate) output_level: Option<i32>,

    /// Write one output file per query into the directory given by -o
    #[arg(long = "multifile", action = clap::ArgAction::SetTrue)]
    pub(crate) output_multifile: bool,
}

fn parse_legacy_format(s: &str) -> Result<OutputFormat, String> {
    match s.parse::<u8>().map_err(|e| e.to_string())? {
        1 => Ok(OutputFormat::Detailed),
        2 => Ok(OutputFormat::Cigar),
        3 => Ok(OutputFormat::BindingSite),
        4 => Ok(OutputFormat::Minimal),
        n => Err(format!("unknown format mode {n}, expected 1–4")),
    }
}

impl TryFrom<OutputArgs> for OutputConfig {
    type Error = Error;

    fn try_from(value: OutputArgs) -> Result<Self, Self::Error> {
        if value.output_multifile && value.path.as_os_str() == "-" {
            bail!("--multifile requires -o/--output to be a directory path; '-' (stdout) is not allowed.");
        }

        if value.path.as_os_str() != "-" {
            validate_output_parent(&value.path)?;
        }

        let codec = value
            .output_compress
            .unwrap_or_else(|| OutputCodec::from(value.path.as_path()));

        if let Some(level) = value.output_level {
            match codec {
                OutputCodec::None => bail!("--compress-level requires compressed output"),
                OutputCodec::Gzip if !(0..=9).contains(&level) => {
                    bail!("gzip level must be 0-9 (got {})", level)
                }
                OutputCodec::Zstd if !(-7..=22).contains(&level) => {
                    bail!("zstd level must be -7..22 (got {})", level)
                }
                _ => {}
            }
        }

        let compress = match (codec, value.output_level) {
            (OutputCodec::None, _) => OutputCompression::None,
            (OutputCodec::Gzip, lvl) => OutputCompression::Gzip(lvl.unwrap_or(6) as u8),
            (OutputCodec::Zstd, lvl) => OutputCompression::Zstd(lvl.unwrap_or(3)),
        };

        let format = match (value.report_format, value.report_legacy) {
            (Some(fmt), Some(_)) => {
                warn!("Both legacy -p/--report-alignment and --format were provided; --format takes precedence.");
                fmt
            }
            (Some(fmt), None) => fmt,
            (None, Some(fmt)) => {
                warn!(
                    "Legacy -p/--report-alignment is deprecated; use --format {}.",
                    fmt.to_possible_value().unwrap().get_name()
                );
                fmt
            }
            (None, None) => OutputFormat::default(),
        };

        Ok(OutputConfig {
            format,
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
            report_format: None,
            report_legacy: None,
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
    fn legacy_report_modes_map_to_formats_and_yield_to_format_flag() {
        assert_eq!(parse_legacy_format("1").unwrap(), OutputFormat::Detailed);
        assert_eq!(parse_legacy_format("2").unwrap(), OutputFormat::Cigar);
        assert_eq!(parse_legacy_format("3").unwrap(), OutputFormat::BindingSite);
        assert_eq!(parse_legacy_format("4").unwrap(), OutputFormat::Minimal);
        assert!(parse_legacy_format("5").is_err());

        let mut legacy = args("out.tsv", None, false);
        legacy.report_legacy = Some(OutputFormat::Cigar);
        assert_eq!(
            OutputConfig::try_from(legacy.clone()).unwrap().format,
            OutputFormat::Cigar
        );
        legacy.report_format = Some(OutputFormat::Minimal);
        assert_eq!(
            OutputConfig::try_from(legacy).unwrap().format,
            OutputFormat::Minimal
        );
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
            (OutputCodec::Zstd, -7, Some(OutputCompression::Zstd(-7))),
            (OutputCodec::Zstd, 22, Some(OutputCompression::Zstd(22))),
            (OutputCodec::Zstd, -8, None),
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
