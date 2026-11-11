use clap::ValueEnum;
use log::warn;
use risearch::config::OutputFormat;
use risearch::dsm::DsmRegistry;
use risearch::types::DsmId;
use std::path::PathBuf;
use std::str::FromStr;

use super::SearchArgs;

/// Legacy CLI parser boundary for `-m/--mismatch`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct LegacyMismatchSpec {
    pub(crate) max_mismatches: usize,
    pub(crate) min_prefix_matches: usize,
    pub(crate) min_suffix_matches: usize,
}

impl FromStr for LegacyMismatchSpec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s
            .trim()
            .split(':')
            .map(|part| {
                part.parse::<usize>()
                    .map_err(|e| format!("invalid mismatch spec '{s}': {e}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let (max_mismatches, min_prefix_matches, min_suffix_matches) = match parts[..] {
            [c] => (c, c, c),
            [c, p] => (c, p, p),
            [c, ps, pe] => (c, ps, pe),
            _ => {
                return Err(format!(
                    "invalid mismatch spec '{s}': expected 'c', 'c:p', or 'c:ps:pe' format"
                ))
            }
        };
        Ok(Self {
            max_mismatches,
            min_prefix_matches,
            min_suffix_matches,
        })
    }
}

/// Legacy CLI parser boundary for `-s/--seed`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct LegacySeedSpec {
    pub(crate) seed_start: Option<i64>,
    pub(crate) seed_end: Option<i64>,
    pub(crate) seed_length: Option<i64>,
}

impl FromStr for LegacySeedSpec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parse =
            |val: &str, field: &str| val.parse::<i64>().map_err(|e| format!("bad {field}: {e}"));

        let s = s.trim();
        let Some((start, rest)) = s.split_once(':') else {
            return Ok(Self {
                seed_start: None,
                seed_end: None,
                seed_length: Some(parse(s, "length")?),
            });
        };
        let (end, length) = match rest.split_once('/') {
            Some((end, length)) => (end, Some(parse(length, "length")?)),
            None => (rest, None),
        };
        Ok(Self {
            seed_start: Some(parse(start, "start")?),
            seed_end: Some(parse(end, "end")?),
            seed_length: length,
        })
    }
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

fn parse_matrix(s: &str) -> Result<DsmId, String> {
    DsmRegistry::parse_id(s).map_err(|e| e.to_string())
}

#[derive(clap::Args, Debug, Clone, Default)]
pub(crate) struct LegacyArgs {
    /// DEPRECATED: legacy alias for -t/--target
    #[arg(
        short = 'i',
        long = "index",
        value_name = "TARGET",
        hide = true,
        conflicts_with = "target"
    )]
    legacy_target: Option<PathBuf>,

    /// DEPRECATED (will be removed in a future release): legacy seed spec
    /// Formats: "l", "m:n", "m:n/l"
    #[arg(
        short = 's',
        long = "seed",
        value_name = "start:end/length",
        conflicts_with_all = ["seed_start", "seed_end", "seed_length"],
        hide = true
    )]
    seed_legacy: Option<LegacySeedSpec>,

    /// DEPRECATED: no-op; wobble is off unless --seed-wobble is given
    #[arg(long = "noGUseed", hide = true, action = clap::ArgAction::SetTrue)]
    no_guseed_legacy: bool,

    /// DEPRECATED (will be removed in a future release): legacy mismatch shorthand
    /// Set max mismatches (c) and min consecutive matches at seed start/end (p)
    /// These seeds will not overlap with perfect complementary seeds.
    /// Prefer --mismatch-max/--mismatch-prefix/--mismatch-suffix.
    #[arg(
        short = 'm',
        long = "mismatch",
        value_name = "c[:ps[:pe]]",
        conflicts_with_all = ["mismatch_max", "mismatch_prefix", "mismatch_suffix"],
        hide = true
    )]
    mismatch_legacy: Option<LegacyMismatchSpec>,

    /// DEPRECATED: Legacy argument for output format (1=detailed, 2=cigar, 3=bindingsite, 4=minimal)
    #[arg(
        short = 'p',
        long = "report-alignment",
        value_name = "MODE",
        num_args = 0..=1,
        default_missing_value = "1",
        value_parser = parse_legacy_format,
        conflicts_with = "report_format",
        hide = true
    )]
    report_legacy: Option<OutputFormat>,

    /// DEPRECATED (will be removed in a future release): use -P/--params or --params-file
    #[arg(
        short = 'z',
        long = "matrix",
        value_name = "MATRIX",
        value_parser = parse_matrix,
        conflicts_with_all = ["dsm_id", "params_file"],
        hide = true
    )]
    matrix_legacy: Option<DsmId>,
}

impl SearchArgs {
    pub(crate) fn resolve_legacy(mut self) -> Self {
        let legacy = std::mem::take(&mut self.legacy);

        if let Some(fmt) = legacy.report_legacy {
            warn!(
                    "Legacy -p/--report-alignment is deprecated and will be removed in a future release; use --format {}.",
                    fmt.to_possible_value().unwrap().get_name()
                );
            self.output.report_format = fmt;
        }

        if legacy.no_guseed_legacy {
            warn!("Legacy --noGUseed is deprecated and will be removed in a future release; it has no effect, since G-U wobble is off unless --seed-wobble is given, so drop it.");
        }

        if let Some(spec) = legacy.seed_legacy {
            let suggestion = match (spec.seed_start, spec.seed_end, spec.seed_length) {
                (None, None, Some(len)) => format!("--seed-length {}", len),
                (Some(start), Some(end), None) => {
                    format!("--seed-start {} --seed-end {}", start, end)
                }
                (Some(start), Some(end), Some(length)) => format!(
                    "--seed-start {} --seed-end {} --seed-length {}",
                    start, end, length
                ),
                _ => "equivalent named flags".into(),
            };
            warn!(
                "Legacy -s/--seed is deprecated and will be removed in a future release; use {}.",
                suggestion
            );
            self.seed.seed_start = spec.seed_start;
            self.seed.seed_end = spec.seed_end;
            self.seed.seed_length = spec.seed_length;
        }

        if let Some(m) = legacy.mismatch_legacy {
            warn!(
                "Legacy -m/--mismatch is deprecated and will be removed in a future release; use --mismatch-max {} --mismatch-prefix {} --mismatch-suffix {}.",
                m.max_mismatches, m.min_prefix_matches, m.min_suffix_matches
            );
            self.seed.mismatch_max = m.max_mismatches;
            self.seed.mismatch_prefix = m.min_prefix_matches;
            self.seed.mismatch_suffix = m.min_suffix_matches;
        }

        if let Some(id) = legacy.matrix_legacy {
            let flag = if DsmRegistry::all_names().contains(&id.0.as_str()) {
                "-P/--params"
            } else {
                "--params-file"
            };
            warn!("Legacy -z/--matrix is deprecated and will be removed in a future release; use {flag} {id}.");
            self.score.dsm_id = id;
        }

        if let Some(target) = legacy.legacy_target {
            warn!("Legacy -i/--index is deprecated and will be removed in a future release; use -t/--target.");
            self.input.target = Some(target);
        }

        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};
    use clap::error::ErrorKind;
    use clap::Parser;
    use risearch::config::{OutputConfig, ScoreConfig, SeedConfig};
    use rstest::rstest;
    use std::path::Path;

    fn parse(args: &[&str]) -> Result<SearchArgs, clap::Error> {
        let base = ["risearch", "search", "-q", "q.fa"];
        let Commands::Search(search) = Cli::try_parse_from([&base[..], args].concat())?.command
        else {
            unreachable!()
        };
        Ok(search)
    }

    #[test]
    fn legacy_i_flag_is_an_alias_for_target() {
        let args = parse(&["-i", "t.idx"]).unwrap();
        assert_eq!(args.input.target, None);
        assert_eq!(
            args.legacy.legacy_target.as_deref(),
            Some(Path::new("t.idx"))
        );
        assert_eq!(
            args.resolve_legacy().input.target.as_deref(),
            Some(Path::new("t.idx"))
        );
    }

    #[test]
    fn parses_seed_cli_specs() {
        let s = LegacySeedSpec::from_str("10").expect("parse");
        assert_eq!(s.seed_start, None);
        assert_eq!(s.seed_end, None);
        assert_eq!(s.seed_length, Some(10));

        let s = LegacySeedSpec::from_str("10:20").expect("parse");
        assert_eq!(s.seed_start, Some(10));
        assert_eq!(s.seed_end, Some(20));
        assert_eq!(s.seed_length, None);

        let s = LegacySeedSpec::from_str("10:20/5").expect("parse");
        assert_eq!(s.seed_start, Some(10));
        assert_eq!(s.seed_end, Some(20));
        assert_eq!(s.seed_length, Some(5));

        assert!(LegacySeedSpec::from_str("abc").is_err());
    }

    #[test]
    fn parses_mismatch_cli_specs() {
        let m = LegacyMismatchSpec::from_str("1").expect("parse");
        assert_eq!(
            (m.max_mismatches, m.min_prefix_matches, m.min_suffix_matches),
            (1, 1, 1)
        );

        let m = LegacyMismatchSpec::from_str("1:3").expect("parse");
        assert_eq!(
            (m.max_mismatches, m.min_prefix_matches, m.min_suffix_matches),
            (1, 3, 3)
        );

        let m = LegacyMismatchSpec::from_str("1:3:5").expect("parse");
        assert_eq!(
            (m.max_mismatches, m.min_prefix_matches, m.min_suffix_matches),
            (1, 3, 5)
        );

        assert!(LegacyMismatchSpec::from_str("1:2:3:4").is_err());
    }

    #[test]
    fn legacy_no_guseed_is_a_noop() {
        let args = parse(&["-t", "t.idx", "--noGUseed"])
            .unwrap()
            .resolve_legacy();
        let config = SeedConfig::from(args.seed);
        assert!(!config.seed_wobble);
    }

    #[test]
    fn legacy_report_modes_map_to_formats_and_conflict_with_format_flag() {
        assert_eq!(parse_legacy_format("1").unwrap(), OutputFormat::Detailed);
        assert_eq!(parse_legacy_format("2").unwrap(), OutputFormat::Cigar);
        assert_eq!(parse_legacy_format("3").unwrap(), OutputFormat::BindingSite);
        assert_eq!(parse_legacy_format("4").unwrap(), OutputFormat::Minimal);
        assert!(parse_legacy_format("5").is_err());

        let format = |args: &[&str]| {
            let args = parse(&[&["-t", "t.idx"][..], args].concat()).unwrap();
            OutputConfig::try_from(args.resolve_legacy().output)
                .unwrap()
                .format
        };
        assert_eq!(format(&["-p", "2"]), OutputFormat::Cigar);
        assert_eq!(
            parse(&["-t", "t.idx", "-p", "2", "-f", "minimal"])
                .unwrap_err()
                .kind(),
            ErrorKind::ArgumentConflict
        );
    }

    fn dsm_id(args: &[&str]) -> Result<DsmId, clap::Error> {
        let args = parse(&[&["-t", "t.idx"][..], args].concat())?;
        Ok(ScoreConfig::from(args.resolve_legacy().score).dsm_id)
    }

    #[rstest]
    #[case::legacy_name(&["-z", "slh04"], "slh04")]
    #[case::legacy_file(&["-z", "Cargo.toml"], "Cargo.toml")]
    fn selected_set_reaches_the_config(#[case] args: &[&str], #[case] expected: &str) {
        assert_eq!(dsm_id(args).unwrap(), DsmId::from(expected));
    }

    #[rstest]
    #[case::legacy_and_params(&["-z", "t04", "--params", "slh04"])]
    #[case::legacy_and_file(&["-z", "t04", "--params-file", "Cargo.toml"])]
    fn invalid_selections_are_rejected(#[case] args: &[&str]) {
        assert_eq!(
            dsm_id(args).unwrap_err().kind(),
            ErrorKind::ArgumentConflict
        );
    }
}
