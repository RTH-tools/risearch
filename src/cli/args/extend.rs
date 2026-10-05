use risearch::config::{ExtendConfig, MAX_EXTENSION, UNLIMITED_EXTENSION};
use risearch::dp::MAX_EXT;

/// Arguments for seed extension strategy
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct ExtendArgs {
    #[arg(
        short = 'l',
        long = "extension",
        value_name = "LENGTH",
        help = format!(
            "Bases to extend on each seed side; 0 keeps the seed alone, {UNLIMITED_EXTENSION} spans \
             the whole query"
        ),
        long_help = format!(
            "Bases to extend on each seed side; 0 keeps the seed alone, {UNLIMITED_EXTENSION} spans \
             the whole query\n\n\
             -l N adds at most N bases on each side; \
             {UNLIMITED_EXTENSION} rejects queries longer than {MAX_EXT} nt, so pass \
             -l {MAX_EXTENSION} or less for those."
        ),
        default_value_t = ExtendConfig::default().max_extension,
        allow_hyphen_values = true,
        value_parser = clap::value_parser!(i32)
            .range(UNLIMITED_EXTENSION as i64..=MAX_EXTENSION as i64)
    )]
    pub(crate) max_extension: i32,
}

impl From<ExtendArgs> for ExtendConfig {
    fn from(value: ExtendArgs) -> Self {
        ExtendConfig {
            max_extension: value.max_extension,
            ..Default::default()
        }
    }
}
