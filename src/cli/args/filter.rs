use risearch::config::FilterConfig;
use risearch::types::Energy;

/// Arguments for filtering and pruning policies
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct FilterArgs {
    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    #[arg(
        short = 'e',
        long = "energy",
        value_name = "dG",
        default_value_t = FilterConfig::default().delta_g,
        allow_hyphen_values = true
    )]
    pub(crate) total_energy: Energy,

    /// Energy per length threshold that filters seeds (in kcal/mol)
    #[arg(
        long = "seed-energy",
        value_name = "THRESHOLD",
        default_value_t = FilterConfig::default().seed_energy
    )]
    pub(crate) seed_energy: Energy,

    /// Report every maximal seed as its own hit instead of collapsing hits
    /// that share a final bounding box to the lowest-energy alignment
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
