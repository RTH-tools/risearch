use risearch::config::{ExtendConfig, MAX_EXTENSION, UNLIMITED_EXTENSION};

/// Arguments for seed extension strategy
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct ExtendArgs {
    /// Max DP extension length up- and downstream of the seed.
    /// -1 extends across the whole query, up to a 256 nt cap per side;
    /// a query longer than that is rejected (pass an explicit -l <=256).
    #[arg(
        short = 'l',
        long = "extension",
        value_name = "LENGTH",
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
