//! Sequence newtype wrapping Vec<Base>
//!
//! This module provides the core `Sequence` type that wraps `Vec<Base>`,
//! ensuring sequences are normalized once at the input boundary and then
//! work with Base enums throughout the codebase.

use crate::types::Base;
use std::ops::{Deref, Index, Range, RangeFrom, RangeFull, RangeTo};

/// A normalized RNA sequence stored as Vec<Base>.
///
/// Sequences are normalized once at input (FASTA parsing) and stored as Base enums.
/// This eliminates repeated ASCII→Base conversions during search and alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence(Vec<Base>);

impl Sequence {
    /// Length of the sequence
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Check if sequence is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Iterator over bases
    pub fn iter(&self) -> impl Iterator<Item = &Base> {
        self.0.iter()
    }

    /// Unchecked access to base at index (for DP hot path).
    ///
    /// # Safety
    /// Caller must ensure `index < self.len()`.
    /// Use only in performance-critical code where bounds are pre-validated.
    #[inline(always)]
    pub unsafe fn get_unchecked(&self, index: usize) -> Base {
        // SAFETY: Caller guarantees index is in bounds
        unsafe { *self.0.get_unchecked(index) }
    }
}

// Implement Deref to allow treating Sequence as &[Base]
impl Deref for Sequence {
    type Target = [Base];

    fn deref(&self) -> &[Base] {
        &self.0
    }
}

// Implement Index to allow seq[i] access
impl Index<usize> for Sequence {
    type Output = Base;

    fn index(&self, index: usize) -> &Base {
        &self.0[index]
    }
}

// Implement Index for ranges
impl Index<Range<usize>> for Sequence {
    type Output = [Base];

    fn index(&self, index: Range<usize>) -> &[Base] {
        &self.0[index]
    }
}

impl Index<RangeFrom<usize>> for Sequence {
    type Output = [Base];

    fn index(&self, index: RangeFrom<usize>) -> &[Base] {
        &self.0[index]
    }
}

impl Index<RangeTo<usize>> for Sequence {
    type Output = [Base];

    fn index(&self, index: RangeTo<usize>) -> &[Base] {
        &self.0[index]
    }
}

impl Index<RangeFull> for Sequence {
    type Output = [Base];

    fn index(&self, index: RangeFull) -> &[Base] {
        &self.0[index]
    }
}

/// Normalizes raw sequence text: case is folded, `T` reads as `U`, `-` and `.`
/// are dropped, any other letter becomes `N`, and any other byte is an error.
impl TryFrom<&[u8]> for Sequence {
    type Error = String;

    fn try_from(raw: &[u8]) -> Result<Self, Self::Error> {
        let mut bases = Vec::with_capacity(raw.len());
        for &b in raw {
            match Base::try_from(b as char) {
                Ok(Base::Gap) => {}
                Ok(base) => bases.push(base),
                Err(_) if b == b'.' => {}
                Err(_) if b.is_ascii_alphabetic() => bases.push(Base::N),
                Err(_) => return Err(format!("byte=0x{b:02X} ('{}')", b as char)),
            }
        }
        Ok(Self(bases))
    }
}

impl From<Vec<Base>> for Sequence {
    fn from(bases: Vec<Base>) -> Self {
        Self(bases)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[test]
    fn normalize_folds_case_and_maps_t_to_u() {
        let seq = Sequence::try_from(b"aCgUtN".as_slice()).unwrap();
        assert_eq!(
            &*seq,
            &[Base::A, Base::C, Base::G, Base::U, Base::U, Base::N]
        );
    }

    #[rstest]
    #[case::gaps_between_bases(b"AC-.GU", &[Base::A, Base::C, Base::G, Base::U])]
    #[case::only_gaps(b"-.", &[])]
    fn normalize_strips_gaps(#[case] raw: &[u8], #[case] expected: &[Base]) {
        let seq = Sequence::try_from(raw).unwrap();
        assert_eq!(&*seq, expected);
    }

    // Alphabetic ambiguity is deliberately accepted;
    // punctuation/non-ASCII bytes are rejected by a separate contract.
    #[rstest]
    #[case::single_code(b"ACRGU", &[Base::A, Base::C, Base::N, Base::G, Base::U])]
    #[case::every_code(b"RYSWKMBDHVXZJryswkmbdhvxzj", &[Base::N; 26])]
    fn normalize_converts_ambiguous_to_n(#[case] raw: &[u8], #[case] expected: &[Base]) {
        let seq = Sequence::try_from(raw).unwrap();
        assert_eq!(&*seq, expected);
    }

    #[rstest]
    #[case::punctuation(b'?')]
    #[case::digit(b'1')]
    #[case::nul(0)]
    #[case::high_byte(0xff)]
    fn rejects_non_base_bytes(#[case] byte: u8) {
        let message = Sequence::try_from([byte].as_slice()).unwrap_err();
        assert!(message.contains(&format!("0x{byte:02X}")));
    }

    #[test]
    fn deref_coercion_to_slice() {
        let seq = Sequence::try_from(b"ACGU".as_slice()).unwrap();
        let slice: &[Base] = &seq; // Deref coercion

        assert_eq!(slice.len(), 4);
        assert_eq!(slice[0], Base::A);
    }

    #[test]
    fn indexing_returns_the_underlying_bases() {
        let seq = Sequence::try_from(b"ACGU".as_slice()).unwrap();

        assert_eq!(seq.len(), 4);
        assert!(!seq.is_empty());
        assert!(Sequence::from(Vec::new()).is_empty());
        assert_eq!(seq[0], Base::A);
        assert_eq!(seq[1], Base::C);
        assert_eq!(seq[2], Base::G);
        assert_eq!(seq[3], Base::U);
        assert_eq!(&seq[1..3], &[Base::C, Base::G]);
        assert_eq!(&seq[2..], &[Base::G, Base::U]);
        assert_eq!(&seq[..2], &[Base::A, Base::C]);
        assert_eq!(&seq[..], &[Base::A, Base::C, Base::G, Base::U]);
    }

    #[test]
    fn iter_yields_all_bases() {
        let seq = Sequence::try_from(b"ACG".as_slice()).unwrap();
        let bases: Vec<Base> = seq.iter().copied().collect();

        assert_eq!(bases, vec![Base::A, Base::C, Base::G]);
    }
}
