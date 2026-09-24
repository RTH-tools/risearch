use risearch::config::{
    ScoreConfig, MAX_PENALTY_KCAL, MAX_TEMPERATURE_C, MIN_PENALTY_KCAL, MIN_TEMPERATURE_C,
};
use risearch::dsm::DsmRegistry;
use risearch::types::{DsmId, Energy};

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

fn parse_matrix(s: &str) -> Result<DsmId, String> {
    DsmRegistry::parse_id(s).map_err(|e| e.to_string())
}

/// Arguments for global scoring model
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct ScoreArgs {
    #[arg(
        short = 'z',
        long = "matrix",
        value_name = "MATRIX",
        help = format!(
            "Dinucleotide stacking model: bundled id ({}) or path to a TSV table",
            DsmRegistry::all_names().join(", ")
        ),
        default_value_t = ScoreConfig::default().dsm_id,
        value_parser = parse_matrix
    )]
    pub(crate) dsm_id: DsmId,

    /// Per-nucleotide extension penalty, included in the reported binding energy (in kcal/mol, 0–50)
    #[arg(
        short = 'd',
        long = "penalty",
        value_name = "PENALTY",
        default_value_t = ScoreConfig::default().penalty,
        value_parser = parse_penalty
    )]
    pub(crate) penalty: Energy,

    #[arg(
        short = 'T',
        long = "temperature",
        value_name = "TEMP",
        help = format!(
            "Temperature for energy calculations (degrees Celsius, 0–100; default {})",
            ScoreConfig::default().temperature
        ),
        value_parser = parse_temperature
    )]
    pub(crate) temperature: Option<i32>,
}

impl From<ScoreArgs> for ScoreConfig {
    fn from(value: ScoreArgs) -> Self {
        ScoreConfig::new(value.dsm_id, value.penalty, value.temperature)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            penalty: Energy::try_from(7.0).unwrap(),
            temperature: None,
        });

        // 37 is also stated in bindings/python/risearch/__init__.py and its README.
        assert_eq!(config.temperature, 37);
    }
}
