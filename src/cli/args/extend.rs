use risearch::config::{ExtendConfig, MAX_EXTENSION, UNLIMITED_EXTENSION};

/// Arguments for seed extension strategy
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct ExtendArgs {
    #[arg(
        short = 'l',
        long = "extension",
        value_name = "LENGTH",
        help = format!(
            "Extension window per seed side in nt, counting the seed's edge base; 0 or 1 keeps the \
             seed alone, {UNLIMITED_EXTENSION} spans the whole query"
        ),
        long_help = format!(
            "Extension window per seed side in nt, counting the seed's edge base; 0 or 1 keeps the \
             seed alone, {UNLIMITED_EXTENSION} spans the whole query\n\n\
             -l N adds at most N-1 bases on each side, as in RIsearch2. {UNLIMITED_EXTENSION} is \
             capped at {MAX_EXTENSION} nt per side; a longer query is rejected, so pass \
             -l {MAX_EXTENSION} or less for it."
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
