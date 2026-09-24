//! Runtime handle over a mapped target index file.

use std::num::NonZeroUsize;
use std::path::Path;

use fs_err::File;
use memmap2::Mmap;

use crate::error::{Error, Result};
use crate::index::archive::{self, ArchivedTargetStore, TargetRecord, TargetStore};
use crate::index::sa::{SuffixIndex, SuffixIndexView};
use crate::index::view::TargetView;
use crate::seq::Sequence;
use crate::types::{Base, Strand};

/// Runtime handle for a target index.
///
/// The mapped file holds the format described by `index::archive`; this
/// type owns the mmap lifetime and carries the native offset directory needed by
/// `TargetView`.
pub struct TargetRegistry {
    mmap: Mmap,
    offsets: Vec<usize>,
}

impl TargetRegistry {
    /// Build a target index over normalized sequences.
    ///
    /// The index is held in an anonymous mapping carrying the same image
    /// [`save`](Self::save) writes, so a built and a reopened index are the same
    /// value. `threads` is the OpenMP team for suffix-array construction with
    /// the `openmp` feature and is ignored without it; `None` uses the OpenMP
    /// default (`OMP_NUM_THREADS`, else every core).
    pub fn build(targets: Vec<(String, Sequence)>, threads: Option<NonZeroUsize>) -> Result<Self> {
        if targets.is_empty() {
            return Err(Error::Input("No target sequences to index".into()));
        }

        let mut records = Vec::with_capacity(targets.len());
        let mut combined_bases = Vec::new();

        for (name, sequence) in targets {
            records.push(TargetRecord {
                name,
                len: sequence.len() as u64,
            });

            combined_bases.reserve(2 * sequence.len() + 2);
            combined_bases.extend(sequence[..].iter().rev().copied());
            combined_bases.push(Base::Gap);
            combined_bases.extend(sequence.iter().copied().map(Base::complement));
            combined_bases.push(Base::Gap);
        }

        let suffix_index = SuffixIndex::build(combined_bases, threads)?;

        let store = TargetStore {
            targets: records,
            sequence: suffix_index.sequence,
            suffix_array: suffix_index.suffix_array,
        };

        let mmap = archive::map_anon(&store)?;
        let source = "a freshly built target index";
        let root = archive::access(archive::payload(mmap.as_ref()), &source)?;
        let offsets = archive::validate(root, &source)?;

        Ok(Self { mmap, offsets })
    }

    /// Write this index to `output` through a temporary file.
    pub fn save(&self, output: &Path) -> Result<()> {
        archive::write(output, self.mmap.as_ref())
    }

    /// Map an index file written by [`save`](Self::save), validating its header
    /// and offset directory.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        // memmap2 drops the path that fs-err would have carried.
        let mmap = unsafe { Mmap::map(&file) }.map_err(|err| {
            std::io::Error::new(
                err.kind(),
                format!("failed to memory-map file `{}`: {err}", path.display()),
            )
        })?;

        archive::read_header(mmap.as_ref(), path)?;
        let source = path.display();
        let root = archive::access(archive::payload(mmap.as_ref()), &source)?;
        let offsets = archive::validate(root, &source)?;

        Ok(Self { mmap, offsets })
    }

    /// Number of indexed targets.
    #[inline]
    pub fn len(&self) -> usize {
        self.root().targets.len()
    }

    /// Whether the index holds no targets.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.root().targets.is_empty()
    }

    /// Name of the target at `idx`.
    #[inline]
    pub fn get_name(&self, idx: usize) -> &str {
        self.root().targets[idx].name.as_str()
    }

    /// Get a seed-search view of the combined SA and sequence.
    pub fn view(&self) -> TargetView<'_> {
        TargetView::new(self.suffixes(), &self.offsets)
    }

    /// Return the selected physical target strand in duplex-column order.
    ///
    /// For input `T` written 5' to 3', Forward is `R(T)` and Reverse is `C(T)`.
    #[inline]
    pub fn target(&self, target_idx: usize, strand: Strand) -> &[Base] {
        self.view().target(target_idx, strand)
    }

    #[inline]
    fn root(&self) -> &ArchivedTargetStore {
        // SAFETY: `open` validates the archive before constructing `TargetRegistry`.
        unsafe { archive::root_unchecked(archive::payload(self.mmap.as_ref())) }
    }

    #[inline]
    fn suffixes(&self) -> SuffixIndexView<'_> {
        let root = self.root();
        let bytes = root.sequence.as_slice();
        let suffix_array = archive::archived_u64_as_native(root.suffix_array.as_slice());
        // SAFETY: `open` validates every stored byte as a Base discriminant.
        unsafe { SuffixIndexView::from_bytes_unchecked(bytes, suffix_array) }
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use tempfile::tempdir;

    use crate::index::sa::SuffixIndex;
    use crate::seq::Sequence;
    use crate::types::{Base, Strand};

    use super::TargetRegistry;

    fn targets(records: &[(&str, &str)]) -> Vec<(String, Sequence)> {
        records
            .iter()
            .map(|(name, seq)| {
                let (sequence, _) = Sequence::normalize(name, seq.as_bytes()).unwrap();
                (name.to_string(), sequence)
            })
            .collect()
    }

    #[test]
    fn roundtrip_build_open_target_view() {
        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");

        let built = TargetRegistry::build(targets(&[("chrA", "ACGUGA"), ("chrB", "UUUGCA")]), None)
            .unwrap();
        built.save(&index_path).unwrap();
        let store = TargetRegistry::open(&index_path).unwrap();
        assert_eq!(store.mmap.as_ref(), built.mmap.as_ref());
        let expected = [("chrA", 6usize), ("chrB", 6usize)];
        assert_eq!(store.len(), expected.len());

        let target = store.view();
        assert!(!target.suffixes().is_empty());
        assert_eq!(store.offsets.len(), store.len());
        let forward = target.map_target_range(0, Strand::Forward, 1..4);
        assert_eq!(forward, 2..5);
        assert_eq!(target.map_target_range(0, Strand::Forward, forward), 1..4);
        assert_eq!(target.map_target_range(0, Strand::Reverse, 1..4), 1..4);

        for (i, (expected_name, expected_seq_len)) in expected.iter().enumerate().take(store.len())
        {
            assert_eq!(store.get_name(i), *expected_name);
            for strand in [Strand::Forward, Strand::Reverse] {
                assert_eq!(store.target(i, strand).len(), *expected_seq_len);
            }
        }

        assert_eq!(
            store.target(0, Strand::Forward),
            &[Base::A, Base::G, Base::U, Base::G, Base::C, Base::A]
        );
        assert_eq!(
            store.target(0, Strand::Reverse),
            &[Base::U, Base::G, Base::C, Base::A, Base::C, Base::U]
        );
    }

    #[test]
    fn map_seed_pos_resolves_the_target_block() {
        let sequence = [Base::Gap; 24];
        let offsets = [0, 12];
        let suffix_index = SuffixIndex::build(sequence.to_vec(), Some(NonZeroUsize::MIN)).unwrap();
        let target = crate::index::view::TargetView::new(suffix_index.view(), &offsets);
        let seed_len = 2;

        for (block_idx, block_start) in offsets.iter().copied().enumerate() {
            assert_eq!(
                target.map_seed_pos(block_start, seed_len),
                Some((block_idx, Strand::Forward, 0))
            );
            assert_eq!(
                target.map_seed_pos(block_start + 6, seed_len),
                Some((block_idx, Strand::Reverse, 0))
            );
            assert_eq!(
                target.map_seed_pos(block_start + 6 + 2, seed_len),
                Some((block_idx, Strand::Reverse, 2))
            );
            assert_eq!(target.map_seed_pos(block_start + 5, seed_len), None);
            assert_eq!(target.map_seed_pos(block_start + 11, seed_len), None);
        }
    }

    #[test]
    fn target_view_maps_positions_in_both_strands() {
        let sequence = [Base::Gap; 12];
        let offsets = [0];
        let suffix_index = SuffixIndex::build(sequence.to_vec(), Some(NonZeroUsize::MIN)).unwrap();
        let target = crate::index::view::TargetView::new(suffix_index.view(), &offsets);
        let seq_len = 5;
        let seed_len = 2;

        assert_eq!(
            target.map_seed_pos(1, seed_len),
            Some((0, Strand::Forward, 1))
        );
        assert_eq!(
            target.map_seed_pos(seq_len + 1 + 2, seed_len),
            Some((0, Strand::Reverse, 2))
        );
        assert_eq!(target.map_seed_pos(seq_len - 1, seed_len), None);
        assert_eq!(target.map_seed_pos(2 * seq_len, seed_len), None);
        assert_eq!(target.map_seed_pos(seq_len, seed_len), None);
        assert_eq!(target.map_seed_pos(2 * seq_len + 1, seed_len), None);
    }

    #[test]
    fn global_offsets_delimit_target_blocks() {
        let store = TargetRegistry::build(
            targets(&[("t1", "ACGU"), ("t2", "GGCC"), ("t3", "AA")]),
            None,
        )
        .unwrap();

        assert_eq!(store.offsets, [0, 10, 20]);
    }
}
