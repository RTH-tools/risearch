use anyhow::Error;
use risearch::config::{ExtendConfig, OutputConfig, SearchConfig};

use super::legacy::LegacyArgs;
use super::{ExtendArgs, FilterArgs, InputArgs, OutputArgs, ScoreArgs, SeedArgs};

/// Arguments that apply to the `search` subcommand
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct SearchArgs {
    #[command(flatten)]
    pub(crate) input: InputArgs,

    #[command(flatten, next_help_heading = "Seed options")]
    pub(crate) seed: SeedArgs,

    #[command(flatten, next_help_heading = "Energy model")]
    pub(crate) score: ScoreArgs,

    #[command(flatten, next_help_heading = "Extension options")]
    pub(crate) extend: ExtendArgs,

    #[command(flatten, next_help_heading = "Filter options")]
    pub(crate) filter: FilterArgs,

    #[command(flatten, next_help_heading = "Output options")]
    pub(crate) output: OutputArgs,

    #[command(flatten)]
    pub(crate) legacy: LegacyArgs,
}

impl SearchArgs {
    /// Resolve CLI-only output policy alongside the frontend-independent search
    /// configuration. Output format controls whether traceback is worth doing,
    /// but the output settings themselves do not belong in [`SearchConfig`].
    pub(crate) fn try_into_configs(self) -> Result<(SearchConfig, OutputConfig), Error> {
        let output: OutputConfig = self.output.try_into()?;
        let mut extend: ExtendConfig = self.extend.into();
        extend.build_alignment = output.format.needs_alignment();

        let search = SearchConfig {
            seed: self.seed.try_into()?,
            score: self.score.into(),
            extend,
            filter: self.filter.into(),
        };
        search.validate().map_err(Error::msg)?;

        Ok((search, output))
    }
}
