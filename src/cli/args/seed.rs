use anyhow::Error;
use log::warn;
use risearch::config::SeedConfig;
use std::str::FromStr;

/// Resolved (seed_start, seed_end, seed_length) bounds from CLI parsing.
type SeedBounds = (Option<i64>, Option<i64>, Option<i64>);

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
        let s = s.trim();
        if s.is_empty() {
            return Err("empty mismatch spec".into());
        }

        let parts: Vec<&str> = s.split(':').collect();
        match parts.len() {
            1 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                Ok(Self {
                    max_mismatches: max,
                    min_prefix_matches: max,
                    min_suffix_matches: max,
                })
            }
            2 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                let min = parts[1]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min consecutive: {}", e))?;
                Ok(Self {
                    max_mismatches: max,
                    min_prefix_matches: min,
                    min_suffix_matches: min,
                })
            }
            3 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                let min_prefix = parts[1]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min prefix matches: {}", e))?;
                let min_suffix = parts[2]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min suffix matches: {}", e))?;
                Ok(Self {
                    max_mismatches: max,
                    min_prefix_matches: min_prefix,
                    min_suffix_matches: min_suffix,
                })
            }
            _ => Err(format!(
                "invalid mismatch spec '{}': expected 'c', 'c:p', or 'c:ps:pe' format",
                s
            )),
        }
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
        let s = s.trim();
        if s.is_empty() {
            return Err("empty seed spec".into());
        }

        let parse_i64 = |val: &str, field: &str| -> Result<i64, String> {
            val.parse::<i64>()
                .map_err(|e| format!("bad {}: {}", field, e))
        };

        match s.find(':') {
            None => {
                let len = parse_i64(s, "length")?;
                Ok(Self {
                    seed_start: None,
                    seed_end: None,
                    seed_length: Some(len),
                })
            }
            Some(colon_idx) => {
                let (start_str, rest) = s.split_at(colon_idx);
                let rest = &rest[1..];

                if rest.is_empty() {
                    return Err("missing end in interval".into());
                }

                let start = parse_i64(start_str, "start")?;

                match rest.find('/') {
                    None => {
                        let end = parse_i64(rest, "end")?;
                        Ok(Self {
                            seed_start: Some(start),
                            seed_end: Some(end),
                            seed_length: None,
                        })
                    }
                    Some(slash_idx) => {
                        let (end_str, len_part) = rest.split_at(slash_idx);
                        let len_str = &len_part[1..];

                        if end_str.is_empty() {
                            return Err("missing end in interval".into());
                        }
                        if len_str.is_empty() {
                            return Err("missing length after '/'".into());
                        }

                        let end = parse_i64(end_str, "end")?;
                        let length = parse_i64(len_str, "length")?;
                        Ok(Self {
                            seed_start: Some(start),
                            seed_end: Some(end),
                            seed_length: Some(length),
                        })
                    }
                }
            }
        }
    }
}

/// Arguments for seed generation
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct SeedArgs {
    /// DEPRECATED (will be removed in a future release): legacy seed spec
    /// Formats: "l", "m:n", "m:n/l"
    #[arg(
        short = 's',
        long = "seed",
        value_name = "start:end/length",
        conflicts_with_all = ["seed_start", "seed_end", "seed_length"],
        help_heading = "Deprecated"
    )]
    pub(crate) seed_legacy: Option<LegacySeedSpec>,

    /// Seed interval start (1-based, can be negative)
    /// TODO: Consider explicit one-sided bounds (e.g. --seed-to-end/--seed-from-start)
    /// instead of inferring missing start/end. Keep strict parsing for now.
    #[arg(
        long = "seed-start",
        value_name = "START",
        allow_hyphen_values = true,
        requires = "seed_end"
    )]
    pub(crate) seed_start: Option<i64>,

    /// Seed interval end (1-based, can be negative)
    #[arg(
        long = "seed-end",
        value_name = "END",
        allow_hyphen_values = true,
        requires = "seed_start"
    )]
    pub(crate) seed_end: Option<i64>,

    /// Seed length (use alone, or with seed-start/seed-end to constrain interval).
    /// If 0 or omitted when using an interval, the full interval width is used.
    /// If omitted completely (no interval provided), defaults to 6.
    #[arg(long = "seed-length", value_name = "LENGTH")]
    pub(crate) seed_length: Option<i64>,

    /// Allow G-U wobble pairs when locating and maximizing seeds
    #[arg(long = "seed-wobble", action = clap::ArgAction::SetTrue)]
    pub(crate) seed_wobble: bool,

    /// Disable maximality check (allows redundant seeds)
    #[arg(long = "no-max-prune", action = clap::ArgAction::SetTrue)]
    pub(crate) no_max_prune: bool,

    /// DEPRECATED: no-op; wobble is off unless --seed-wobble is given
    #[arg(long = "noGUseed", hide = true, action = clap::ArgAction::SetTrue)]
    pub(crate) no_guseed_legacy: bool,

    /// DEPRECATED (will be removed in a future release): legacy mismatch shorthand
    /// Set max mismatches (c) and min consecutive matches at seed start/end (p)
    /// These seeds will not overlap with perfect complementary seeds.
    /// Prefer --mismatch-max/--mismatch-prefix/--mismatch-suffix.
    #[arg(
        short = 'm',
        long = "mismatch",
        value_name = "c[:ps[:pe]]",
        conflicts_with_all = ["mismatch_max", "mismatch_prefix", "mismatch_suffix"],
        help_heading = "Deprecated"
    )]
    pub(crate) mismatch_legacy: Option<LegacyMismatchSpec>,

    /// Max number of mismatches allowed in the seed (preferred)
    #[arg(long = "mismatch-max", value_name = "C")]
    pub(crate) mismatch_max: Option<usize>,

    /// Min consecutive matches at seed start (prefix / 5')
    #[arg(long = "mismatch-prefix", value_name = "PS", requires = "mismatch_max")]
    pub(crate) mismatch_prefix: Option<usize>,

    /// Min consecutive matches at seed end (suffix / 3')
    #[arg(
        long = "mismatch-suffix",
        value_name = "PE",
        requires_all = ["mismatch_max", "mismatch_prefix"]
    )]
    pub(crate) mismatch_suffix: Option<usize>,
}

impl SeedArgs {
    fn resolve_mismatches(&self) -> Result<(usize, usize, usize), String> {
        match (self.mismatch_max, self.mismatch_prefix, self.mismatch_suffix) {
            (None, None, None) => {
                let defaults = SeedConfig::default();
                Ok((
                    defaults.max_mismatches,
                    defaults.min_prefix_matches,
                    defaults.min_suffix_matches,
                ))
            }
            (Some(max), None, None) => Ok((max, max, max)),
            (Some(max), Some(prefix), None) => Ok((max, prefix, prefix)),
            (Some(max), Some(prefix), Some(suffix)) => Ok((max, prefix, suffix)),
            _ => Err(
                "use mismatch-max alone, or mismatch-max + mismatch-prefix (optionally with mismatch-suffix)"
                    .into(),
            ),
        }
    }

    fn resolve_seed_bounds(&self) -> Result<SeedBounds, String> {
        match (self.seed_start, self.seed_end, self.seed_length) {
            (Some(s), Some(e), length) => Ok((Some(s), Some(e), length)),
            (None, None, length) => Ok((None, None, length)),
            _ => Err(
                "use seed_length alone, or seed_start + seed_end (optionally with seed_length)"
                    .into(),
            ),
        }
    }
}

impl TryFrom<SeedArgs> for SeedConfig {
    type Error = Error;

    fn try_from(value: SeedArgs) -> Result<Self, Self::Error> {
        if value.no_guseed_legacy {
            warn!("'--noGUseed' is deprecated and now a no-op; G-U wobble is off by default, enable it with --seed-wobble.");
        }

        // clap makes -s exclusive with the named bounds, so the legacy branch is
        // reached only when it is the sole source.
        let (seed_start, seed_end, seed_length) = if value.seed_start.is_some()
            || value.seed_end.is_some()
            || value.seed_length.is_some()
        {
            value.resolve_seed_bounds().map_err(Error::msg)?
        } else if let Some(spec) = value.seed_legacy.as_ref() {
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
            warn!("Legacy -s/--seed is deprecated; use {}.", suggestion);
            (spec.seed_start, spec.seed_end, spec.seed_length)
        } else {
            (None, None, None)
        };

        let (max_mismatches, min_prefix_matches, min_suffix_matches) = if let Some(m) =
            value.mismatch_legacy.as_ref()
        {
            warn!(
                    "Legacy -m/--mismatch is deprecated; use --mismatch-max {} --mismatch-prefix {} --mismatch-suffix {}.",
                    m.max_mismatches, m.min_prefix_matches, m.min_suffix_matches
                );
            (m.max_mismatches, m.min_prefix_matches, m.min_suffix_matches)
        } else {
            value.resolve_mismatches().map_err(Error::msg)?
        };

        let config = SeedConfig {
            seed_start,
            seed_end,
            seed_length,
            seed_wobble: value.seed_wobble,
            no_max_prune: value.no_max_prune,
            max_mismatches,
            min_prefix_matches,
            min_suffix_matches,
        };
        config.validate().map_err(Error::msg)?;
        Ok(config)
    }
}

#[cfg(test)]
mod tests {
    use super::{LegacyMismatchSpec, LegacySeedSpec, SeedArgs, SeedBounds};
    use risearch::config::SeedConfig;
    use std::str::FromStr;

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

    fn test_args_mismatch(
        max: Option<usize>,
        prefix: Option<usize>,
        suffix: Option<usize>,
    ) -> Result<(usize, usize, usize), String> {
        let args = SeedArgs {
            seed_legacy: None,
            seed_start: None,
            seed_end: None,
            seed_length: None,
            seed_wobble: false,
            no_guseed_legacy: false,
            no_max_prune: false,
            mismatch_legacy: None,
            mismatch_max: max,
            mismatch_prefix: prefix,
            mismatch_suffix: suffix,
        };
        args.resolve_mismatches()
    }

    fn test_args_seed(
        start: Option<i64>,
        end: Option<i64>,
        length: Option<i64>,
    ) -> Result<SeedBounds, String> {
        let args = SeedArgs {
            seed_legacy: None,
            seed_start: start,
            seed_end: end,
            seed_length: length,
            seed_wobble: false,
            no_guseed_legacy: false,
            no_max_prune: false,
            mismatch_legacy: None,
            mismatch_max: None,
            mismatch_prefix: None,
            mismatch_suffix: None,
        };
        args.resolve_seed_bounds()
    }

    #[test]
    fn builds_mismatch_from_named_args() {
        assert_eq!(
            test_args_mismatch(None, None, None).expect("parse"),
            (0, 1, 0)
        );
        assert_eq!(
            test_args_mismatch(Some(1), None, None).expect("parse"),
            (1, 1, 1)
        );
        assert_eq!(
            test_args_mismatch(Some(1), Some(3), None).expect("parse"),
            (1, 3, 3)
        );
        assert_eq!(
            test_args_mismatch(Some(1), Some(3), Some(5)).expect("parse"),
            (1, 3, 5)
        );
    }

    #[test]
    fn rejects_partial_named_mismatch_args() {
        assert!(test_args_mismatch(None, Some(3), None).is_err());
        assert!(test_args_mismatch(None, None, Some(5)).is_err());
        assert!(test_args_mismatch(Some(1), None, Some(5)).is_err());
    }

    #[test]
    fn builds_seed_from_named_args() {
        assert_eq!(
            test_args_seed(None, None, None).expect("parse"),
            (None, None, None)
        );
        assert_eq!(
            test_args_seed(None, None, Some(6)).expect("parse"),
            (None, None, Some(6))
        );
        assert_eq!(
            test_args_seed(Some(3), Some(12), None).expect("parse"),
            (Some(3), Some(12), None)
        );
        assert_eq!(
            test_args_seed(Some(3), Some(12), Some(7)).expect("parse"),
            (Some(3), Some(12), Some(7))
        );
    }

    #[test]
    fn rejects_partial_named_seed_args() {
        assert!(test_args_seed(Some(3), None, None).is_err());
        assert!(test_args_seed(None, Some(12), None).is_err());
    }

    #[test]
    fn seed_args_conversion_is_fallible_for_invalid_partial_seed_bounds() {
        let args = SeedArgs {
            seed_legacy: None,
            seed_start: Some(3),
            seed_end: None,
            seed_length: None,
            seed_wobble: false,
            no_guseed_legacy: false,
            no_max_prune: false,
            mismatch_legacy: None,
            mismatch_max: None,
            mismatch_prefix: None,
            mismatch_suffix: None,
        };

        assert!(SeedConfig::try_from(args).is_err());
    }

    #[test]
    fn legacy_no_guseed_is_a_noop() {
        let args = SeedArgs {
            seed_legacy: None,
            seed_start: None,
            seed_end: None,
            seed_length: None,
            seed_wobble: false,
            no_guseed_legacy: true,
            no_max_prune: false,
            mismatch_legacy: None,
            mismatch_max: None,
            mismatch_prefix: None,
            mismatch_suffix: None,
        };
        let config = SeedConfig::try_from(args).unwrap();
        assert!(!config.seed_wobble);
    }
}
