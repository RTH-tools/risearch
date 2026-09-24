//! Suffix-index construction and borrowed traversal views.
//!
//! A suffix index keeps an indexed sequence and its lexicographically ordered
//! suffix positions together. Consumers must not mix the positions from one
//! index with the sequence from another.

use std::num::NonZeroUsize;
use std::ops::Range;

use libsais::SuffixArrayConstruction;
#[cfg(debug_assertions)]
use rayon::prelude::*;

use crate::error::Error;
use crate::types::Base;

/// An owned sequence and the lexicographically ordered suffixes built over it.
///
/// `suffix_array` may contain either every suffix or an order-preserving subset
/// selected for a particular search. Every retained position refers to
/// `sequence`; the two fields therefore form one logical value and must not be
/// separated.
///
/// The sequence is stored as bytes so it can be handed to the archive writer
/// and reinterpreted from a mapping without copying. All bytes are valid
/// [`Base`] discriminants.
///
/// This is the in-memory form only. The persisted form lives in
/// [`crate::index::archive`] and always holds the complete suffix array; the
/// subset case below is never archived.
pub(crate) struct SuffixIndex {
    pub(super) sequence: Vec<u8>,
    pub(super) suffix_array: Vec<u64>,
}

/// A borrowed sequence paired with its lexicographically ordered suffixes.
///
/// This is the traversal-facing form of [`SuffixIndex`]. The suffix slice may be
/// an order-preserving subset of the complete suffix array. Every position must
/// still refer to the paired sequence.
#[derive(Clone, Copy)]
pub(crate) struct SuffixIndexView<'a> {
    sequence: &'a [Base],
    suffix_array: &'a [u64],
}

impl<'a> SuffixIndexView<'a> {
    /// Project a view from an invariant-bearing byte representation.
    ///
    /// This reinterprets the sequence bytes as [`Base`] values without copying.
    /// It is restricted to owners such as [`SuffixIndex`] and validated archive
    /// boundaries that already establish the sequence/suffix-array pairing.
    ///
    /// # Safety
    ///
    /// Every byte in `sequence` must be a valid [`Base`] discriminant.
    #[inline]
    pub(super) unsafe fn from_bytes_unchecked(sequence: &'a [u8], suffix_array: &'a [u64]) -> Self {
        let sequence =
            unsafe { std::slice::from_raw_parts(sequence.as_ptr().cast(), sequence.len()) };
        Self {
            sequence,
            suffix_array,
        }
    }

    /// Return the indexed sequence.
    #[inline]
    pub(crate) fn sequence(self) -> &'a [Base] {
        self.sequence
    }

    /// Return the number of retained suffixes.
    #[inline]
    pub(crate) fn len(self) -> usize {
        self.suffix_array.len()
    }

    /// Return whether the view contains no retained suffixes.
    #[inline]
    pub(crate) fn is_empty(self) -> bool {
        self.suffix_array.is_empty()
    }

    /// Return a range of positions in lexicographic suffix order.
    ///
    /// # Panics
    ///
    /// Panics if `interval` is not a valid range within `0..self.len()`.
    #[inline]
    pub(crate) fn suffix_positions(self, interval: Range<usize>) -> &'a [u64] {
        &self.suffix_array[interval]
    }

    /// Read `sequence[suffix_array[sa_idx] + depth]`, or `Base::Gap` past its end.
    ///
    /// A depth beyond the sequence reads as `Gap` so descent stops there; a
    /// suffix array that is not sorted cannot walk off the sequence.
    ///
    /// # Safety
    ///
    /// `sa_idx` must be in bounds.
    #[inline(always)]
    pub(crate) unsafe fn base_unchecked(self, sa_idx: usize, depth: usize) -> Base {
        let suffix_pos = unsafe { *self.suffix_array.get_unchecked(sa_idx) as usize };
        self.sequence
            .get(suffix_pos + depth)
            .copied()
            .unwrap_or(Base::Gap)
    }
}

impl SuffixIndex {
    /// Build the complete suffix array over exactly `sequence`.
    ///
    /// No sentinel or padding is added. Callers that require a terminal
    /// [`Base::Gap`] must include it in `sequence` before calling this method.
    ///
    /// When the `openmp` feature is enabled, `threads` selects the libsais
    /// worker count and `None` uses the OpenMP default. Without that feature,
    /// `threads` is ignored.
    ///
    /// # Errors
    ///
    /// Returns an error if libsais cannot construct the suffix array.
    pub(crate) fn build(sequence: Vec<Base>, threads: Option<NonZeroUsize>) -> Result<Self, Error> {
        let suffix_array = build_suffix_array(&sequence, threads)?;
        Ok(Self {
            sequence: sequence.into_iter().map(Base::as_u8).collect(),
            suffix_array,
        })
    }

    /// Build the seed-traversal index over `seed`, terminated by a [`Base::Gap`].
    ///
    /// Suffixes too short to reach `min_len` are dropped, which also drops the
    /// terminal `Gap` suffix; `min_len` must be at least 1.
    ///
    /// # Errors
    ///
    /// Returns an error if libsais cannot construct the suffix array.
    pub(crate) fn build_for_seed(seed: &[Base], min_len: usize) -> Result<Self, Error> {
        let mut sequence = Vec::with_capacity(seed.len() + 1);
        sequence.extend_from_slice(seed);
        sequence.push(Base::Gap);

        // A per-query seed is tens of bases; an OpenMP team would cost more than the sort.
        let mut index = Self::build(sequence, Some(NonZeroUsize::MIN))?;
        let max_valid_start = seed.len().saturating_sub(min_len) as u64;
        index
            .suffix_array
            .retain(|&position| position <= max_valid_start);
        Ok(index)
    }

    /// Borrow the sequence and its current suffix selection as one view.
    #[inline]
    pub(crate) fn view(&self) -> SuffixIndexView<'_> {
        // SAFETY: both constructors store only valid Base discriminants.
        unsafe { SuffixIndexView::from_bytes_unchecked(&self.sequence, &self.suffix_array) }
    }
}

fn build_suffix_array(bases: &[Base], threads: Option<NonZeroUsize>) -> Result<Vec<u64>, Error> {
    // SAFETY: Base is #[repr(u8)] and its enum discriminants perfectly
    // match the required suffix array lexicographical sort order.
    let sort_bytes: &[u8] =
        unsafe { std::slice::from_raw_parts(bases.as_ptr().cast::<u8>(), bases.len()) };

    // Build the suffix array. Use OpenMP if feature is enabled.
    #[cfg(feature = "openmp")]
    let sa_raw = {
        // Cap at u16::MAX, the libsais count width.
        let thread_count = match threads {
            Some(n) => libsais::ThreadCount::fixed(n.get().min(u16::MAX as usize) as u16),
            None => libsais::ThreadCount::openmp_default(),
        };
        SuffixArrayConstruction::for_text(sort_bytes)
            .in_owned_buffer()
            .multi_threaded(thread_count)
            .run()
            .map_err(|e| Error::Index(format!("{e:?}")))?
    };

    #[cfg(not(feature = "openmp"))]
    let sa_raw = {
        let _ = threads; // thread count only applies with the `openmp` feature
        SuffixArrayConstruction::for_text(sort_bytes)
            .in_owned_buffer()
            .single_threaded()
            .run()
            .map_err(|e| Error::Index(format!("{e:?}")))?
    };

    // Extract the SA vec
    let mut sa: Vec<i64> = sa_raw.into_vec();

    // Verify no negative indices were generated by libsais before casting
    #[cfg(debug_assertions)]
    {
        sa.par_iter().for_each(|&pos_i64| {
            debug_assert!(pos_i64 >= 0, "libsais returned negative suffix position");
        });
    }

    // SAFETY: i64 and u64 have identical memory layout (size and alignment).
    // We take ownership of the buffer, stop the old Vec from running its Drop logic,
    // and reconstruct a new Vec<u64> from the exact same memory region.
    // This eliminates an O(N) allocation, slashing peak memory usage by 50%.
    let positions: Vec<u64> = unsafe {
        let ptr = sa.as_mut_ptr() as *mut u64;
        let len = sa.len();
        let cap = sa.capacity();
        std::mem::forget(sa);
        Vec::from_raw_parts(ptr, len, cap)
    };

    Ok(positions)
}
