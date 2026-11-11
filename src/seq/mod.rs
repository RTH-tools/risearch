//! Stores sequences as normalized `Base` values.
//!
//! `Sequence` holds a `Vec<Base>`, where `Base` is A, C, G, U, N or Gap.
//! `Sequence::try_from(&[u8])` normalizes raw bytes once, when they are read.
//! Code elsewhere in the crate works on `Base` values, not on text.

pub mod sequence;

pub use sequence::Sequence;
