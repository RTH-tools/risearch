use risearch::config::FilterConfig;
use risearch::types::Energy;

/// Arguments for filtering and pruning policies
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct FilterArgs {
    /// Report hits with a binding energy at or below this, in kcal/mol
    #[arg(
        short = 'e',
        long = "energy",
        value_name = "dG",
        default_value_t = FilterConfig::default().delta_g,
        allow_hyphen_values = true
    )]
    pub(crate) total_energy: Energy,

    /// Accepted for compatibility; the search does not use it yet
    #[arg(
        long = "seed-energy",
        value_name = "THRESHOLD",
        default_value_t = FilterConfig::default().seed_energy,
        allow_hyphen_values = true
    )]
    pub(crate) seed_energy: Energy,

    /// Report one hit per seed instead of one per duplex
    ///
    /// By default, hits whose extended duplexes cover the same query and target span on the same
    /// strand collapse to the one with the lowest energy.
    #[arg(long = "no-dedup", action = clap::ArgAction::SetTrue)]
    pub(crate) no_dedup: bool,
}

impl From<FilterArgs> for FilterConfig {
    fn from(value: FilterArgs) -> Self {
        FilterConfig {
            delta_g: value.total_energy,
            seed_energy: value.seed_energy,
            no_dedup: value.no_dedup,
        }
    }
}
