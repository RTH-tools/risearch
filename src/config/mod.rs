//! Frontend-independent search and output configuration.
//!
//! [`SearchConfig`] is what [`run_search`](crate::run_search) consumes;
//! [`OutputConfig`] covers only what the output writer needs. A caller can
//! build either by hand without going through the CLI:
//! [`SearchConfig::validate`] checks the parts that have invalid states
//! ([`SeedConfig`], [`ScoreConfig`], [`ExtendConfig`]), and the bounds it
//! enforces are the constants below.

use std::path::Path;

use crate::dp::MAX_EXT;
use crate::dsm::DsmRegistry;
use crate::error::{Error, Result};
use crate::types::{DsmId, Energy};

// `Default` implementations below own frontend default policy. Keep constants
// here only for effective fallbacks and shared validity bounds or sentinels.
/// Seed length applied when [`SeedConfig::seed_length`] is `None` and no
/// interval was given.
pub const DEFAULT_SEED_LEN: i64 = 6;

/// Lower bound on [`ScoreConfig::penalty`], in kcal/mol.
pub const MIN_PENALTY_KCAL: f64 = 0.0;
/// Upper bound on [`ScoreConfig::penalty`], in kcal/mol.
pub const MAX_PENALTY_KCAL: f64 = 50.0;
/// Lower bound on [`ScoreConfig::temperature`], in °C.
pub const MIN_TEMPERATURE_C: i32 = 0;
/// Upper bound on [`ScoreConfig::temperature`], in °C.
pub const MAX_TEMPERATURE_C: i32 = 50;

/// [`ExtendConfig::max_extension`] sentinel: extend across the whole query.
pub const UNLIMITED_EXTENSION: i32 = -1;
/// Largest [`ExtendConfig::max_extension`] the DP grid accepts per side.
pub const MAX_EXTENSION: i32 = MAX_EXT as i32;

/// Path that stands for stdin when reading and stdout when writing.
pub const STDIO: &str = "-";

// =============================================================================
// ENUMS (shared by config and CLI via clap derives)
// =============================================================================

/// Layout of each reported hit.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
#[cfg_attr(
    feature = "cli",
    derive(clap::ValueEnum),
    clap(rename_all = "lowercase")
)]
pub enum OutputFormat {
    /// Minimal columns, each hit preceded by a three-line alignment
    Detailed,
    /// Minimal columns plus a pairing string
    Cigar,
    /// Cigar columns plus the aligned target and its 5' and 3' flanks
    BindingSite,
    /// Query, start, end, target, start, end, strand and energy
    #[default]
    Minimal,
}

impl OutputFormat {
    /// Whether this format prints pairing data, and so needs the extension to
    /// run traceback. The single source of truth for
    /// [`ExtendConfig::build_alignment`].
    pub const fn needs_alignment(self) -> bool {
        !matches!(self, Self::Minimal)
    }
}

/// Config-facing compression value — codec and level bound together.
/// Illegal states (level without codec) are unrepresentable.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum OutputCompression {
    /// Uncompressed.
    #[default]
    None,
    /// gzip at the given level, 0–9.
    Gzip(u8),
    /// zstd at the given level, within [`zstd::compression_level_range`].
    Zstd(i32),
}

impl OutputCompression {
    /// Return a suitable file extension (including leading dot) for this compression
    /// configuration, assuming TSV format, e.g. `".tsv"`, `".tsv.gz"`, `".tsv.zst"`.
    pub fn extension(&self) -> &'static str {
        match self {
            Self::None => ".tsv",
            Self::Gzip(_) => ".tsv.gz",
            Self::Zstd(_) => ".tsv.zst",
        }
    }
}

/// Arguments for seed generation.
#[derive(Debug, Clone)]
pub struct SeedConfig {
    /// Seed interval start (1-based, can be negative). `None` = full query.
    pub seed_start: Option<i64>,
    /// Seed interval end (1-based, can be negative). `None` = full query.
    pub seed_end: Option<i64>,
    /// Seed length. In length-only mode this is the desired seed length.
    /// In interval mode, `None` means use the full interval width.
    pub seed_length: Option<i64>,

    /// Allow G-U wobble pairs when locating and maximizing seeds. Off by default.
    pub seed_wobble: bool,

    /// Emit non-maximal seeds too — those that could be grown by one more
    /// pairing column, and so are shorter copies of a longer match.
    pub no_max_prune: bool,

    /// Maximum number of mismatches allowed in seed. Mismatch-bearing seeds
    /// containing an uninterrupted pairing run of the minimum seed length are
    /// omitted: that region is represented by perfect seeds. This selection
    /// applies even with [`Self::no_max_prune`].
    pub max_mismatches: usize,
    /// Minimum consecutive matches at seed start (prefix, 5').
    pub min_prefix_matches: usize,
    /// Minimum consecutive matches at seed end (suffix, 3').
    pub min_suffix_matches: usize,
}

impl Default for SeedConfig {
    fn default() -> Self {
        Self {
            seed_start: None,
            seed_end: None,
            seed_length: None,
            seed_wobble: false,
            no_max_prune: false,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        }
    }
}

impl SeedConfig {
    /// Validate query-independent seed invariants.
    ///
    /// Bounds against a particular sequence length, and a requested length
    /// against that normalized interval, are checked by [`Self::resolve`].
    pub fn validate(&self) -> Result<()> {
        match (self.seed_start, self.seed_end) {
            (None, None) => {
                if self.seed_length.is_some_and(|length| length <= 0) {
                    return Err(Error::Config(
                        "Invalid seed length: must be positive without an interval".into(),
                    ));
                }
            }
            (Some(start), Some(end)) => {
                if start == 0 || end == 0 {
                    return Err(Error::Config(
                        "Invalid seed interval: coordinates are 1-based and cannot be zero".into(),
                    ));
                }
                if start.signum() != end.signum() {
                    return Err(Error::Config("Invalid seed interval: mixed sign".into()));
                }
                if end < start {
                    return Err(Error::Config(
                        "Invalid seed interval: end precedes start".into(),
                    ));
                }

                // In interval mode, zero or a negative length means the full
                // interval, preserving the legacy CLI's established behavior.
                if let Some(length) = self.seed_length.filter(|length| *length > 0) {
                    let interval_len = i128::from(end) - i128::from(start) + 1;
                    if i128::from(length) > interval_len {
                        return Err(Error::Config(
                            "Invalid seed length (exceeds interval)".into(),
                        ));
                    }
                }
            }
            _ => {
                return Err(Error::Config(
                    "use seed_start + seed_end together, or neither".into(),
                ))
            }
        }

        Ok(())
    }

    /// Resolve this config's seed bounds against a specific query length.
    pub fn resolve(&self, query_len: usize) -> Result<(usize, usize, usize)> {
        self.validate()?;

        let n = query_len as i64;
        if n <= 0 {
            return Err(Error::Config("query length must be positive".into()));
        }

        match (self.seed_start, self.seed_end) {
            (None, None) => {
                let effective_length = self.seed_length.unwrap_or(DEFAULT_SEED_LEN);
                if effective_length <= 0 {
                    return Err(Error::Config("Invalid seed length".into()));
                }
                let length = effective_length.min(n) as usize;
                Ok((1, query_len, length))
            }
            (Some(start), Some(end)) => {
                let (s_pos, e_pos) = match (start.signum(), end.signum()) {
                    (1, 1) if end >= start && end <= n => (start as usize, end as usize),
                    (-1, -1) if end >= start && start >= -n => {
                        let to_pos = |v: i64| -> usize {
                            let x = v + 1;
                            let r = (x + n) % n;
                            if r == 0 {
                                n as usize
                            } else {
                                r as usize
                            }
                        };
                        (to_pos(start), to_pos(end))
                    }
                    _ => {
                        return Err(Error::Config(
                            "Invalid seed interval: outside query bounds".into(),
                        ))
                    }
                };

                if s_pos == 0 || e_pos == 0 {
                    return Err(Error::Config("Invalid seed interval".into()));
                }

                let interval_len = e_pos.saturating_sub(s_pos) + 1;
                if interval_len == 0 {
                    return Err(Error::Config("Invalid seed interval: empty".into()));
                }

                let final_len = match self.seed_length {
                    None => interval_len,
                    // 0 or negative length explicitly requested -> use full interval
                    Some(l) if l <= 0 => interval_len,
                    Some(l) if (l as usize) > interval_len => {
                        return Err(Error::Config(
                            "Invalid seed length (exceeds interval)".into(),
                        ))
                    }
                    Some(l) => l as usize,
                };

                Ok((s_pos, e_pos, final_len))
            }
            _ => Err(Error::Config(
                "use seed_start + seed_end together, or neither".into(),
            )),
        }
    }
}

/// Global scoring model shared by seed scoring and DP extension.
#[derive(Debug, Clone)]
pub struct ScoreConfig {
    /// Which bundled dinucleotide stacking model to score with.
    pub dsm_id: DsmId,
    /// Per-nucleotide extension penalty, included in the reported binding energy.
    pub penalty: Energy,
    /// Temperature in °C the model is evaluated at.
    pub temperature: i32,
}

impl Default for ScoreConfig {
    fn default() -> Self {
        Self {
            dsm_id: DsmId::from("t04"),
            penalty: Energy::default(),
            temperature: 37,
        }
    }
}

impl ScoreConfig {
    /// Build a scoring config. `None` takes the default temperature; an explicit
    /// one is kept, with a warning when `dsm_id` is a user TSV table, which is
    /// used as-is.
    pub fn new(dsm_id: DsmId, penalty: Energy, temperature: Option<i32>) -> Self {
        let custom = !DsmRegistry::all_names().contains(&dsm_id.0.as_str())
            && Path::new(&dsm_id.0).is_file();
        if temperature.is_some() && custom {
            log::warn!(
                "temperature has no effect on the custom energy parameter file '{dsm_id}'; it is used as-is."
            );
        }
        ScoreConfig {
            dsm_id,
            penalty,
            temperature: temperature.unwrap_or(Self::default().temperature),
        }
    }

    /// Check the model id resolves and the penalty and temperature are in range.
    pub fn validate(&self) -> Result<()> {
        DsmRegistry::parse_id(self.dsm_id.0.as_str())?;

        let penalty = self.penalty.to_kcal();
        if !(MIN_PENALTY_KCAL..=MAX_PENALTY_KCAL).contains(&penalty) {
            return Err(Error::Config(format!(
                "penalty must be between {MIN_PENALTY_KCAL} and {MAX_PENALTY_KCAL}, got {penalty}"
            )));
        }
        if !(MIN_TEMPERATURE_C..=MAX_TEMPERATURE_C).contains(&self.temperature) {
            return Err(Error::Config(format!(
                "temperature must be between {MIN_TEMPERATURE_C} and {MAX_TEMPERATURE_C}, got {}",
                self.temperature
            )));
        }

        Ok(())
    }
}

/// Arguments for extension strategy.
#[derive(Debug, Clone)]
pub struct ExtendConfig {
    /// Max extension length on the seed (do DP for max this length up- and downstream of seed).
    /// Negative = unlimited: extend across the whole query. Queries longer than
    /// MAX_EXT per side are rejected rather than silently clamped.
    pub max_extension: i32,

    /// Record the pairing alignment during DP traceback, populating
    /// [`SearchHit::alignment`](crate::SearchHit). Off skips that work entirely.
    ///
    /// Dedup's tie-break compares fingerprints, so with this off a different
    /// member of an energy-tied bounding box may survive. Tied hits share a box
    /// and an energy, so every other field is identical either way.
    pub build_alignment: bool,
}

impl Default for ExtendConfig {
    fn default() -> Self {
        Self {
            max_extension: 20,
            build_alignment: true,
        }
    }
}

impl ExtendConfig {
    /// Whether the extension window is unlimited (extends across the whole query).
    pub fn is_unlimited(&self) -> bool {
        self.max_extension == UNLIMITED_EXTENSION
    }

    /// The fixed window size, or `None` for unlimited (follow the query).
    pub fn max_window(&self) -> Option<usize> {
        (!self.is_unlimited()).then_some(self.max_extension as usize)
    }

    /// Check `max_extension` is the unlimited sentinel or within the DP cap.
    pub fn validate(&self) -> Result<()> {
        if !(UNLIMITED_EXTENSION..=MAX_EXTENSION).contains(&self.max_extension) {
            return Err(Error::Config(format!(
                "max extension must be {UNLIMITED_EXTENSION} (unlimited) or between 0 and {MAX_EXTENSION}, got {}",
                self.max_extension
            )));
        }
        Ok(())
    }
}

/// Hit acceptance and pruning policies.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    pub delta_g: Energy,

    /// Energy per length threshold that filters seeds
    pub seed_energy: Energy,

    /// Report every maximal seed as its own hit. When unset (the default),
    /// hits whose extension resolves to the same final bounding box are
    /// collapsed to the single lowest-energy alignment.
    pub no_dedup: bool,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            delta_g: Energy::from_kcal(-20.0),
            seed_energy: Energy::default(),
            no_dedup: false,
        }
    }
}

/// Options that apply to the `search` subcommand
#[derive(Debug, Clone, Default)]
pub struct SearchConfig {
    /// Where seeds are taken from and how they may mismatch.
    pub seed: SeedConfig,
    /// Scoring model, penalty, and temperature.
    pub score: ScoreConfig,
    /// How far a seed may extend, and whether traceback runs.
    pub extend: ExtendConfig,
    /// Which finished hits are reported.
    pub filter: FilterConfig,
}

impl SearchConfig {
    /// Validate every part. `filter` has no invalid states of its own.
    pub fn validate(&self) -> Result<()> {
        self.seed.validate()?;
        self.score.validate()?;
        self.extend.validate()?;
        Ok(())
    }
}

/// Where and how hits are written. Separate from [`SearchConfig`] because none
/// of it changes which hits are found.
#[derive(Debug, Clone)]
pub struct OutputConfig {
    /// Output format
    pub format: OutputFormat,

    /// Compression codec and level, bound together (resolved at arg boundary)
    pub compress: OutputCompression,

    /// Write one output file per query (directory mode)
    pub multifile: bool,
}

#[cfg(test)]
mod spec_tests;
