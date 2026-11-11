//! Sequence newtype wrapping Vec<Base>
//!
//! This module provides the core `Sequence` type that wraps `Vec<Base>`,
//! ensuring sequences are normalized once at the input boundary and then
//! work with Base enums throughout the codebase.

use crate::types::Base;
use std::ops::Deref;

/// A normalized RNA sequence stored as Vec<Base>.
///
/// Sequences are normalized once at input (FASTA parsing) and stored as Base enums.
/// This eliminates repeated ASCII→Base conversions during search and alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence(Vec<Base>);

// Implement Deref to allow treating Sequence as &[Base]
impl Deref for Sequence {
    type Target = [Base];

    fn deref(&self) -> &[Base] {
        &self.0
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
}
