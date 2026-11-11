use anyhow::Error;
use risearch::config::{SeedConfig, DEFAULT_SEED_LEN};

/// Arguments for seed generation
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct SeedArgs {
    /// First query position a seed may cover (1-based; negative counts from the 3' end, -1 is the
    /// last base)
    ///
    /// Requires --seed-end with the same sign; both ends are inclusive. Seeds must lie entirely
    /// inside the interval and, without --seed-length, span all of it.
    #[arg(
        long = "seed-start",
        value_name = "START",
        allow_hyphen_values = true,
        requires = "seed_end"
    )]
    pub(crate) seed_start: Option<i64>,

    /// Last query position a seed may cover (inclusive; same sign as --seed-start)
    #[arg(
        long = "seed-end",
        value_name = "END",
        allow_hyphen_values = true,
        requires = "seed_start"
    )]
    pub(crate) seed_end: Option<i64>,

    #[arg(
        long = "seed-length",
        value_name = "LENGTH",
        help = format!(
            "Minimum seed length [default: {DEFAULT_SEED_LEN}, or the interval's full width with --seed-start/--seed-end]"
        ),
        long_help = format!(
            "Minimum seed length\n\n\
             A seed is a stretch of complementary pairs shared by query and target; only seeds at least \
             LENGTH long are extended. A query shorter than LENGTH uses its full length. With an \
             interval, 0 also means its full width, and LENGTH may not exceed it.\n\n\
             [default: {DEFAULT_SEED_LEN}, or the interval's full width with --seed-start/--seed-end]"
        )
    )]
    pub(crate) seed_length: Option<i64>,

    /// Allow G-U wobble pairs in seeds (extension is not affected)
    #[arg(long = "seed-wobble", action = clap::ArgAction::SetTrue)]
    pub(crate) seed_wobble: bool,

    /// Also extend non-maximal seeds, which produces redundant hits
    ///
    /// By default a seed that one more pairing column would lengthen is dropped, because it is a
    /// shorter copy of a longer seed.
    #[arg(
        long = "no-max-prune",
        action = clap::ArgAction::SetTrue,
        hide_short_help = true
    )]
    pub(crate) no_max_prune: bool,

    /// Mismatches allowed in a seed
    ///
    /// A seed with mismatches is dropped when it contains an uninterrupted run of the minimum seed
    /// length, since perfect seeds already cover that region.
    #[arg(
        long = "mismatch-max",
        value_name = "COUNT",
        default_value_t = SeedConfig::default().max_mismatches
    )]
    pub(crate) mismatch_max: usize,

    /// Matches required at the seed start (5') before any mismatch
    #[arg(
        long = "mismatch-prefix",
        value_name = "LENGTH",
        requires = "mismatch_max",
        default_value_t = SeedConfig::default().min_prefix_matches,
        hide_short_help = true
    )]
    pub(crate) mismatch_prefix: usize,

    /// Matches required at the seed end (3') after any mismatch
    #[arg(
        long = "mismatch-suffix",
        value_name = "LENGTH",
        requires = "mismatch_max",
        default_value_t = SeedConfig::default().min_suffix_matches,
        hide_short_help = true
    )]
    pub(crate) mismatch_suffix: usize,
}

impl TryFrom<SeedArgs> for SeedConfig {
    type Error = Error;

    fn try_from(value: SeedArgs) -> Result<Self, Self::Error> {
        let config = SeedConfig {
            seed_start: value.seed_start,
            seed_end: value.seed_end,
            seed_length: value.seed_length,
            seed_wobble: value.seed_wobble,
            no_max_prune: value.no_max_prune,
            max_mismatches: value.mismatch_max,
            min_prefix_matches: value.mismatch_prefix,
            min_suffix_matches: value.mismatch_suffix,
        };
        config.validate().map_err(Error::msg)?;
        Ok(config)
    }
}

#[cfg(test)]
mod tests {
    use super::SeedArgs;
    use risearch::config::SeedConfig;

    #[test]
    fn seed_args_conversion_is_fallible_for_invalid_partial_seed_bounds() {
        let args = SeedArgs {
            seed_start: Some(3),
            seed_end: None,
            seed_length: None,
            seed_wobble: false,
            no_max_prune: false,
            mismatch_max: 0,
            mismatch_prefix: 1,
            mismatch_suffix: 1,
        };

        assert!(SeedConfig::try_from(args).is_err());
    }
}
