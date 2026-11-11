//! Core domain types used across the codebase

/// Nucleotide/gap representation for DSM indexing and sequence operations.
///
/// The discriminants are the canonical internal rank. DSM tables, DP scoring,
/// and suffix-array construction/partitioning all depend on this exact order.
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum Base {
    /// Gap column. Rank 0 doubles as the DSM lookup and SA-partitioning sentinel.
    #[default]
    Gap = 0,
    /// Adenine.
    A = 1,
    /// Cytosine.
    C = 2,
    /// Guanine.
    G = 3,
    /// Ambiguous or unknown base.
    N = 4,
    /// Uracil. `T` parses to this variant as well.
    U = 5,
}

/// Base → uppercase ASCII byte
static BASE_TO_UPPER: [u8; 6] = *b"-ACGNU";

/// Base → lowercase ASCII byte (for to_byte())
static BASE_TO_BYTE: [u8; 6] = *b"-acgnu";

impl TryFrom<char> for Base {
    type Error = String;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c.to_ascii_uppercase() {
            'A' => Ok(Base::A),
            'C' => Ok(Base::C),
            'G' => Ok(Base::G),
            'U' | 'T' => Ok(Base::U),
            'N' => Ok(Base::N),
            '-' => Ok(Base::Gap),
            _ => Err(format!("invalid base '{c}'")),
        }
    }
}

impl Base {
    /// Convert to array index
    #[inline]
    pub const fn as_usize(self) -> usize {
        self as usize
    }

    /// Convert raw u8 byte representation back to Base.
    #[inline]
    pub fn from_u8(v: u8) -> Self {
        match v {
            0 => Base::Gap,
            1 => Base::A,
            2 => Base::C,
            3 => Base::G,
            4 => Base::N,
            5 => Base::U,
            _ => panic!("invalid base rank {v}"),
        }
    }

    /// Convert raw u8 byte representation.
    #[inline]
    pub const fn as_u8(self) -> u8 {
        self as u8
    }

    /// Get standard uppercase ASCII byte (A, G, C, U, N)
    #[inline]
    pub fn to_u8_upper(self) -> u8 {
        BASE_TO_UPPER[self.as_usize()]
    }

    /// Convert Base to lowercase ASCII byte (a, g, c, t, n, -)
    /// Used for display and output formatting.
    #[inline]
    pub fn to_byte(self) -> u8 {
        BASE_TO_BYTE[self.as_usize()]
    }

    /// Watson-Crick complement (A <-> U/T, G <-> C)
    #[inline(always)]
    pub const fn complement(self) -> Self {
        match self {
            Base::Gap => Base::Gap,
            Base::A => Base::U,
            Base::G => Base::C,
            Base::C => Base::G,
            Base::U => Base::A,
            Base::N => Base::N,
        }
    }

    /// True for the four matchable bases (A, C, G, U); false for N and Gap.
    #[inline(always)]
    pub(crate) const fn is_matchable(self) -> bool {
        matches!(self, Base::A | Base::C | Base::G | Base::U)
    }

    /// Classify pairing between two bases.
    ///
    /// Watson-Crick: A↔U, C↔G.  Wobble: G↔U.  Everything else: mismatch.
    #[inline(always)]
    pub(crate) const fn pair_type(self, other: Base) -> PairType {
        match (self, other) {
            (Base::A, Base::U) | (Base::U, Base::A) | (Base::C, Base::G) | (Base::G, Base::C) => {
                PairType::Canonical
            }
            (Base::G, Base::U) | (Base::U, Base::G) => PairType::Wobble,
            _ => PairType::Mismatch,
        }
    }
}

/// Chemistry of a base pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PairType {
    /// Watson-Crick: A↔U, C↔G
    Canonical,
    /// G↔U wobble
    Wobble,
    /// Non-pairing
    Mismatch,
}

impl PairType {
    /// Whether this pair type counts as a match under the given wobble mode.
    #[inline(always)]
    pub(crate) const fn is_match(self, wobble: bool) -> bool {
        match self {
            PairType::Canonical => true,
            PairType::Wobble => wobble,
            PairType::Mismatch => false,
        }
    }
}

/// Number of nucleotide types (Gap, A, C, G, N, U)
pub(crate) const BASE_COUNT: usize = 6;

/// Constant for Gap index used in array indexing and DSM lookups.
pub(crate) const GAP: u8 = Base::Gap.as_u8();

/// Strand direction for search

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    /// Reported as `+`. Selects `R(T)` in duplex-column order.
    Forward,
    /// Reported as `-`. Selects `C(T)` in duplex-column order.
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", char::from(*self))
    }
}

impl TryFrom<char> for Strand {
    type Error = String;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            '+' => Ok(Strand::Forward),
            '-' => Ok(Strand::Reverse),
            _ => Err(format!("invalid strand '{c}'")),
        }
    }
}

impl From<Strand> for char {
    fn from(s: Strand) -> char {
        match s {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }
}

/// Represents thermodynamic energy, internally stored as an integer (scaling kcal/mol by 10,000).
/// All conversions from floating-point values expect inputs in kcal/mol.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Energy(pub(crate) i32);
impl Energy {
    pub(crate) const SCALE: f64 = 10000.0;

    /// Most-negative threshold. `energy <= Energy::MIN` is ~never true, so used
    /// as a filter bound it records nothing — the integer-domain stand-in for
    /// a `-inf` threshold (which `from_kcal` rightly rejects as non-finite).
    pub const MIN: Energy = Energy(i32::MIN);

    /// Create Energy from a value in kcal/mol. Panics if the value is non-finite or overflows.
    /// Use this for trusted internal constants and tests.
    pub fn from_kcal(kcal: f64) -> Self {
        Self::try_from(kcal).expect("invalid energy value")
    }

    /// Convert the internal integer representation back to kcal/mol.
    pub fn to_kcal(self) -> f64 {
        f64::from(self.0) / Self::SCALE
    }
}

impl TryFrom<f64> for Energy {
    type Error = String;

    /// Create Energy from a value in kcal/mol.
    /// Returns an error if the value is non-finite or would overflow the internal representation.
    fn try_from(kcal: f64) -> Result<Self, Self::Error> {
        if !kcal.is_finite() {
            return Err(format!("non-finite energy: {kcal}"));
        }
        let score = kcal * Self::SCALE;
        if score < f64::from(i32::MIN) || score > f64::from(i32::MAX) {
            return Err(format!("energy overflow: {kcal} kcal/mol"));
        }
        Ok(Energy(score.round() as i32))
    }
}

impl From<Energy> for f64 {
    fn from(e: Energy) -> Self {
        e.to_kcal()
    }
}

impl std::str::FromStr for Energy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let kcal = s.parse::<f64>().map_err(|e| e.to_string())?;
        Self::try_from(kcal)
    }
}

impl std::ops::Add for Energy {
    type Output = Self;

    #[track_caller]
    fn add(self, rhs: Self) -> Self::Output {
        Self(
            self.0
                .checked_add(rhs.0)
                .expect("Energy addition overflowed i32"),
        )
    }
}

impl std::ops::Sub for Energy {
    type Output = Self;

    #[track_caller]
    fn sub(self, rhs: Self) -> Self::Output {
        Self(
            self.0
                .checked_sub(rhs.0)
                .expect("Energy subtraction overflowed i32"),
        )
    }
}

impl std::ops::Mul<usize> for Energy {
    type Output = Self;

    #[track_caller]
    fn mul(self, rhs: usize) -> Self::Output {
        let rhs_i32 = i32::try_from(rhs).expect("Energy multiplier overflows i32");
        Self(
            self.0
                .checked_mul(rhs_i32)
                .expect("Energy multiplication overflowed i32"),
        )
    }
}

impl std::fmt::Display for Energy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.2}", self.to_kcal())
    }
}

/// Dinucleotide stacking model (DSM) selector: a bundled identifier from
/// `DsmRegistry::all_names()` or a path to a TSV table.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct DsmId(pub String);

impl std::fmt::Display for DsmId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for DsmId {
    fn from(s: &str) -> Self {
        Self(s.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::{Base, DsmId, Energy, PairType, Strand};

    const BASES: [Base; 6] = [Base::Gap, Base::A, Base::C, Base::G, Base::N, Base::U];
    const CANONICAL: [(Base, Base); 4] = [
        (Base::A, Base::U),
        (Base::U, Base::A),
        (Base::C, Base::G),
        (Base::G, Base::C),
    ];
    const WOBBLE: [(Base, Base); 2] = [(Base::G, Base::U), (Base::U, Base::G)];

    #[test]
    fn base_discriminants_are_canonical_internal_rank() {
        assert_eq!(Base::Gap.as_u8(), 0);
        assert_eq!(Base::A.as_u8(), 1);
        assert_eq!(Base::C.as_u8(), 2);
        assert_eq!(Base::G.as_u8(), 3);
        assert_eq!(Base::N.as_u8(), 4);
        assert_eq!(Base::U.as_u8(), 5);
    }

    #[test]
    fn from_u8_round_trips_every_rank() {
        for base in BASES {
            assert_eq!(Base::from_u8(base.as_u8()), base);
        }
    }

    #[test]
    #[should_panic(expected = "invalid base rank 6")]
    fn from_u8_rejects_ranks_beyond_the_alphabet() {
        Base::from_u8(6);
    }

    #[test]
    fn upper_bytes_are_the_display_alphabet() {
        assert_eq!(Base::Gap.to_u8_upper(), b'-');
        assert_eq!(Base::A.to_u8_upper(), b'A');
        assert_eq!(Base::C.to_u8_upper(), b'C');
        assert_eq!(Base::G.to_u8_upper(), b'G');
        assert_eq!(Base::N.to_u8_upper(), b'N');
        assert_eq!(Base::U.to_u8_upper(), b'U');
    }

    #[test]
    fn lower_bytes_are_the_output_alphabet() {
        assert_eq!(BASES.map(Base::to_byte), *b"-acgnu");
    }

    #[test]
    fn chars_parse_case_insensitively_with_t_as_u() {
        for (c, base) in [
            ('A', Base::A),
            ('C', Base::C),
            ('G', Base::G),
            ('U', Base::U),
            ('T', Base::U),
            ('N', Base::N),
            ('-', Base::Gap),
        ] {
            assert_eq!(Base::try_from(c), Ok(base));
            assert_eq!(Base::try_from(c.to_ascii_lowercase()), Ok(base));
        }
        assert_eq!(Base::try_from('x'), Err("invalid base 'x'".to_string()));
    }

    #[test]
    fn complement_swaps_watson_crick_partners_and_fixes_the_rest() {
        assert_eq!(
            BASES.map(Base::complement),
            [Base::Gap, Base::U, Base::G, Base::C, Base::N, Base::A]
        );
    }

    #[test]
    fn only_the_four_nucleotides_are_matchable() {
        for base in [Base::A, Base::C, Base::G, Base::U] {
            assert!(base.is_matchable());
        }
        assert!(!Base::N.is_matchable());
        assert!(!Base::Gap.is_matchable());
    }

    #[test]
    fn every_pair_classifies_by_the_literal_pairing_table() {
        for q in BASES {
            for t in BASES {
                let expected = if CANONICAL.contains(&(q, t)) {
                    PairType::Canonical
                } else if WOBBLE.contains(&(q, t)) {
                    PairType::Wobble
                } else {
                    PairType::Mismatch
                };
                assert_eq!(q.pair_type(t), expected, "{q:?}-{t:?}");
            }
        }
    }

    #[test]
    fn wobble_counts_as_a_match_only_when_enabled() {
        assert!(PairType::Canonical.is_match(false));
        assert!(PairType::Canonical.is_match(true));
        assert!(!PairType::Wobble.is_match(false));
        assert!(PairType::Wobble.is_match(true));
        assert!(!PairType::Mismatch.is_match(false));
        assert!(!PairType::Mismatch.is_match(true));
    }

    #[test]
    fn strand_round_trips_through_its_symbol() {
        assert_eq!(Strand::try_from('+'), Ok(Strand::Forward));
        assert_eq!(Strand::try_from('-'), Ok(Strand::Reverse));
        assert_eq!(Strand::Forward.to_string(), "+");
        assert_eq!(Strand::Reverse.to_string(), "-");
    }

    #[test]
    fn strand_rejects_other_symbols_and_orders_forward_first() {
        assert_eq!(Strand::try_from('x'), Err("invalid strand 'x'".to_string()));
        assert!(Strand::Forward < Strand::Reverse);
    }

    #[test]
    fn energy_accepts_the_exact_i32_bounds_and_rejects_beyond() {
        assert!(Energy::try_from(f64::from(i32::MIN) / Energy::SCALE).is_ok());
        assert!(Energy::try_from(f64::from(i32::MAX) / Energy::SCALE).is_ok());
        assert!(Energy::try_from(-1e6).is_err());
        assert!(Energy::try_from(1e6).is_err());
    }

    #[test]
    fn energy_scales_kcal_by_ten_thousand_and_rounds_away_from_zero() {
        assert_eq!(Energy::try_from(-1.2345), Ok(Energy(-12345)));
        assert_eq!(Energy::try_from(0.00006), Ok(Energy(1)));
        assert_eq!(Energy::try_from(-0.00006), Ok(Energy(-1)));
        assert_eq!(Energy::try_from(0.00004), Ok(Energy(0)));
        assert_eq!(Energy::from_kcal(-1.2345).to_kcal(), -1.2345);
        assert_eq!(f64::from(Energy(-12345)), -1.2345);
    }

    #[test]
    fn energy_rejects_non_finite_input() {
        assert_eq!(
            Energy::try_from(f64::NAN),
            Err("non-finite energy: NaN".to_string())
        );
        assert_eq!(
            Energy::try_from(f64::INFINITY),
            Err("non-finite energy: inf".to_string())
        );
        assert_eq!(
            "inf".parse::<Energy>(),
            Err("non-finite energy: inf".to_string())
        );
    }

    #[test]
    fn energy_parses_kcal_text() {
        assert_eq!("-1.5".parse::<Energy>(), Ok(Energy(-15000)));
        assert_eq!(
            "abc".parse::<Energy>(),
            Err("invalid float literal".to_string())
        );
    }

    #[test]
    fn energy_arithmetic_is_checked_integer_arithmetic() {
        assert_eq!(Energy(3) + Energy(4), Energy(7));
        assert_eq!(Energy(3) - Energy(4), Energy(-1));
        assert_eq!(Energy(-3) * 4, Energy(-12));
        assert!(Energy(-1) < Energy(0));
        assert!(Energy::MIN < Energy(i32::MIN + 1));
        for overflow in [
            (|| Energy(i32::MAX) + Energy(1)) as fn() -> Energy,
            || Energy(i32::MIN) - Energy(1),
            || Energy(i32::MAX) * 2,
            || Energy(1) * usize::MAX,
        ] {
            assert!(std::panic::catch_unwind(overflow).is_err());
        }
    }

    #[test]
    fn energy_displays_two_decimals_of_kcal() {
        assert_eq!(Energy(-12300).to_string(), "-1.23");
        assert_eq!(Energy(123456).to_string(), "12.35");
        assert_eq!(Energy(0).to_string(), "0.00");
    }

    #[test]
    fn dsm_id_displays_its_name() {
        assert_eq!(DsmId::from("t04").to_string(), "t04");
        assert_eq!(DsmId::default(), DsmId(String::new()));
    }
}
