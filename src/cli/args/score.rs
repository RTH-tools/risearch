use clap::builder::{PossibleValue, PossibleValuesParser, TypedValueParser};
use risearch::config::{
    ScoreConfig, MAX_PENALTY_KCAL, MAX_TEMPERATURE_C, MIN_PENALTY_KCAL, MIN_TEMPERATURE_C,
};
use risearch::dsm::DsmRegistry;
use risearch::types::{DsmId, Energy};
use std::path::Path;

fn parse_penalty(s: &str) -> Result<Energy, String> {
    let v: f64 = s.parse().map_err(|e| format!("{e}"))?;
    if !(MIN_PENALTY_KCAL..=MAX_PENALTY_KCAL).contains(&v) {
        return Err(format!(
            "penalty must be between {} and {}, got {v}",
            MIN_PENALTY_KCAL, MAX_PENALTY_KCAL
        ));
    }
    Energy::try_from(v)
}

fn parse_temperature(s: &str) -> Result<i32, String> {
    let v: i32 = s.parse().map_err(|e| format!("{e}"))?;
    if !(MIN_TEMPERATURE_C..=MAX_TEMPERATURE_C).contains(&v) {
        return Err(format!(
            "temperature must be between {} and {}, got {v}",
            MIN_TEMPERATURE_C, MAX_TEMPERATURE_C
        ));
    }
    Ok(v)
}

fn parse_params_file(s: &str) -> Result<DsmId, String> {
    if !Path::new(s).is_file() {
        return Err("not an existing file".into());
    }
    std::path::absolute(s)
        .map_err(|e| e.to_string())?
        .into_os_string()
        .into_string()
        .map(DsmId)
        .map_err(|p| format!("{p:?} is not valid UTF-8"))
}

/// Arguments for global scoring model
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct ScoreArgs {
    /// Nearest-neighbor energy parameter set (use --params-file for a TSV table)
    #[arg(
        short = 'P',
        long = "params",
        value_name = "NAME",
        default_value_t = ScoreConfig::default().dsm_id,
        value_parser = PossibleValuesParser::new(
            DsmRegistry::all_models()
                .iter()
                .map(|&(name, help)| PossibleValue::new(name).help(help))
        ).map(DsmId),
        conflicts_with = "params_file"
    )]
    pub(crate) dsm_id: DsmId,

    /// Read energy parameters from a TSV table, used as-is at any temperature; cannot be combined
    /// with --params
    ///
    /// Tab-separated, with the header q1, q2, t1, t2, delta_g_kcal_per_mol. Each row is the
    /// stacking free energy in kcal/mol of query bases q1 q2 (5' to 3') over target bases t1 t2
    /// (3' to 5'). The README describes the format in full.
    #[arg(long = "params-file", value_name = "FILE", value_parser = parse_params_file)]
    pub(crate) params_file: Option<DsmId>,

    #[arg(
        short = 'd',
        long = "penalty",
        value_name = "PENALTY",
        help = format!(
            "Penalty per duplex nucleotide, on both strands and including the seed, added to the reported binding energy (in kcal/mol, {MIN_PENALTY_KCAL}–{MAX_PENALTY_KCAL})"
        ),
        default_value_t = ScoreConfig::default().penalty,
        value_parser = parse_penalty
    )]
    pub(crate) penalty: Energy,

    #[arg(
        short = 'T',
        long = "temperature",
        value_name = "TEMP",
        help = format!(
            "Temperature for energy calculations (degrees Celsius, {MIN_TEMPERATURE_C}–{MAX_TEMPERATURE_C}; ignored by custom TSV tables) [default: {}]",
            ScoreConfig::default().temperature
        ),
        value_parser = parse_temperature
    )]
    pub(crate) temperature: Option<i32>,
}

impl From<ScoreArgs> for ScoreConfig {
    fn from(value: ScoreArgs) -> Self {
        ScoreConfig::new(
            value.params_file.unwrap_or(value.dsm_id),
            value.penalty,
            value.temperature,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};
    use clap::error::ErrorKind;
    use clap::Parser;
    use rstest::rstest;

    #[test]
    fn penalty_parses_inside_the_range_and_is_rejected_outside() {
        assert_eq!(parse_penalty("5").unwrap(), Energy::try_from(5.0).unwrap());
        for bad in ["-1", "51", "abc"] {
            assert!(parse_penalty(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn score_args_carry_every_field_into_the_config() {
        let config = ScoreConfig::from(ScoreArgs {
            dsm_id: DsmId::from("slh04"),
            params_file: None,
            penalty: Energy::try_from(7.0).unwrap(),
            temperature: Some(42),
        });

        assert_eq!(config.dsm_id, DsmId::from("slh04"));
        assert_eq!(config.penalty, Energy::try_from(7.0).unwrap());
        assert_eq!(config.temperature, 42);
    }

    #[test]
    fn omitted_temperature_falls_back_to_the_config_default() {
        let config = ScoreConfig::from(ScoreArgs {
            dsm_id: DsmId::from("t04"),
            params_file: None,
            penalty: Energy::try_from(7.0).unwrap(),
            temperature: None,
        });

        // 37 is also stated in bindings/python/risearch/__init__.py and its README.
        assert_eq!(config.temperature, 37);
    }

    fn dsm_id(args: &[&str]) -> Result<DsmId, clap::Error> {
        let base = ["risearch", "search", "-q", "q.fa", "-t", "t.idx"];
        let Commands::Search(search) = Cli::try_parse_from([&base[..], args].concat())?.command
        else {
            unreachable!()
        };
        Ok(ScoreConfig::from(search.score).dsm_id)
    }

    #[rstest]
    #[case::default(&[], "t04")]
    #[case::params(&["--params", "slh04"], "slh04")]
    fn selected_set_reaches_the_config(#[case] args: &[&str], #[case] expected: &str) {
        assert_eq!(dsm_id(args).unwrap(), DsmId::from(expected));
    }

    #[test]
    fn params_file_reaches_the_config_as_an_absolute_path() {
        let expected = std::env::current_dir().unwrap().join("Cargo.toml");
        assert_eq!(
            dsm_id(&["--params-file", "Cargo.toml"]).unwrap(),
            DsmId::from(expected.to_str().unwrap())
        );
    }

    #[rstest]
    #[case::missing_file(&["--params-file", "missing.tsv"], ErrorKind::ValueValidation)]
    #[case::params_and_file(&["--params", "t04", "--params-file", "Cargo.toml"], ErrorKind::ArgumentConflict)]
    fn invalid_selections_are_rejected(#[case] args: &[&str], #[case] kind: ErrorKind) {
        assert_eq!(dsm_id(args).unwrap_err().kind(), kind);
    }
}
