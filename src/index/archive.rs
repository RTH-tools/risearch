//! On-disk format for the target index.
//!
//! A file is a fixed [`HEADER_LEN`] byte header followed by an rkyv archive of
//! [`TargetStore`]. The header carries everything needed to reject a file before
//! its bytes are interpreted as an archive, which an in-archive version field
//! cannot do: a layout change breaks validation before any field is reachable.
//!
//! Targets are stored as lengths rather than offsets. Block offsets are the
//! prefix sum of `2 * len + 2`, so a monotonic, zero-based, gap-free directory
//! is the only representable one.

use std::fmt::Display;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use memmap2::{Mmap, MmapOptions};
use rayon::prelude::*;
use rkyv::ser::writer::{Buffer, IoWriter};
use rkyv::ser::Positional;

use crate::error::{Error, Result};
use crate::types::Base;

/// Bytes preceding the rkyv payload.
///
/// A multiple of 16, so a page-aligned mapping leaves the payload aligned for
/// rkyv.
pub(super) const HEADER_LEN: usize = 32;

const MAGIC: [u8; 8] = *b"RIS3IDX\0";
const FORMAT_VERSION: u32 = 4;
const LITTLE_ENDIAN: u8 = 0;

/// Width in bytes of one stored suffix position.
const SA_WIDTH: u8 = 8;

#[derive(rkyv::Archive, rkyv::Serialize)]
pub(super) struct TargetStore {
    pub(super) targets: Vec<TargetRecord>,
    pub(super) sequence: Vec<u8>,
    pub(super) suffix_array: Vec<u64>,
}

#[derive(rkyv::Archive, rkyv::Serialize)]
pub(super) struct TargetRecord {
    pub(super) name: String,
    pub(super) len: u64,
}

/// Reject `bytes` unless its header describes an archive this build can read.
pub(super) fn read_header(bytes: &[u8], path: &Path) -> Result<()> {
    let rebuild = || format!("{}; rebuild it with `risearch index`", path.display());

    if bytes.len() < HEADER_LEN {
        return Err(Error::Index(format!(
            "Target index is too small to hold a header: {}",
            rebuild()
        )));
    }
    if bytes[..MAGIC.len()] != MAGIC {
        return Err(Error::Index(format!(
            "Not a risearch target index: {}",
            rebuild()
        )));
    }

    let version = u32::from_le_bytes(bytes[8..12].try_into().expect("4 bytes"));
    if version != FORMAT_VERSION {
        return Err(Error::Index(format!(
            "Unsupported target index version {version} in {}",
            rebuild()
        )));
    }
    if bytes[12] != SA_WIDTH {
        return Err(Error::Index(format!(
            "Target index stores {}-byte suffix positions, expected {SA_WIDTH}: {}",
            bytes[12],
            rebuild()
        )));
    }
    if bytes[13] != LITTLE_ENDIAN {
        return Err(Error::Index(format!(
            "Target index is not little-endian: {}",
            rebuild()
        )));
    }

    let payload_len = u64::from_le_bytes(bytes[16..24].try_into().expect("8 bytes"));
    let expected = HEADER_LEN as u64 + payload_len;
    if expected != bytes.len() as u64 {
        return Err(Error::Index(format!(
            "Target index is {} bytes, its header declares {expected}: {}",
            bytes.len(),
            rebuild()
        )));
    }

    Ok(())
}

/// Borrow the archive payload that follows a validated header.
#[inline]
pub(super) fn payload(bytes: &[u8]) -> &[u8] {
    &bytes[HEADER_LEN..]
}

/// Validate `payload` as an archive and borrow its root.
pub(super) fn access<'a>(
    payload: &'a [u8],
    source: &dyn Display,
) -> Result<&'a ArchivedTargetStore> {
    rkyv::access::<ArchivedTargetStore, rkyv::rancor::BoxedError>(payload)
        .map_err(|err| Error::Index(format!("Invalid target index archive: {source}: {err}")))
}

/// Borrow the root of a payload already validated by [`access`].
///
/// # Safety
///
/// `payload` must have passed [`access`].
#[inline]
pub(super) unsafe fn root_unchecked(payload: &[u8]) -> &ArchivedTargetStore {
    unsafe { rkyv::access_unchecked::<ArchivedTargetStore>(payload) }
}

/// Check the invariants the archive layout cannot express, returning the block
/// offset of each target.
pub(super) fn validate(root: &ArchivedTargetStore, source: &dyn Display) -> Result<Vec<usize>> {
    let seq = root.sequence.as_slice();
    let sa = root.suffix_array.as_slice();
    if seq.len() != sa.len() {
        return Err(Error::Index(format!(
            "Target index sequence/SA length mismatch (seq={}, sa={}): {}",
            seq.len(),
            sa.len(),
            source
        )));
    }

    // Reinterpreting these bytes as `Base` requires every discriminant to be
    // valid, and no byte pattern is invalid for the `u8` bytecheck sees. Reduced
    // as a per-chunk max rather than a short-circuiting `any`, which cannot
    // vectorize; at index scale that is ~0.1 versus ~0.8 G instructions.
    let max_rank = seq
        .par_chunks(1 << 20)
        .map(|chunk| chunk.iter().copied().fold(0u8, u8::max))
        .max()
        .unwrap_or(0);
    if max_rank > Base::U.as_u8() {
        return Err(Error::Index(format!(
            "Target index contains an invalid base rank: {source}"
        )));
    }

    let gap = Base::Gap.as_u8();
    let mut offsets = Vec::with_capacity(root.targets.len());
    let mut cursor = 0usize;
    for target in root.targets.iter() {
        let name = target.name.as_str();
        let len = usize::try_from(target.len.to_native())
            .map_err(|_| Error::Index(format!("Target '{name}' length does not fit in usize")))?;
        let block_end = len
            .checked_mul(2)
            .and_then(|doubled| doubled.checked_add(1))
            .and_then(|last| cursor.checked_add(last))
            .ok_or_else(|| {
                Error::Index(format!("Target '{name}' block overflows the address space"))
            })?;
        if block_end >= seq.len() {
            return Err(Error::Index(format!(
                "Target '{name}' block exceeds the indexed sequence (end={block_end}, len={}): {}",
                seq.len(),
                source
            )));
        }
        if seq[cursor + len] != gap || seq[block_end] != gap {
            return Err(Error::Index(format!(
                "Target '{name}' block is missing a separator Gap in {source}"
            )));
        }
        offsets.push(cursor);
        cursor = block_end + 1;
    }

    if cursor != seq.len() {
        return Err(Error::Index(format!(
            "Target blocks cover {cursor} bytes of a {}-byte indexed sequence: {}",
            seq.len(),
            source
        )));
    }

    Ok(offsets)
}

/// Encode the fixed header that precedes `payload_len` payload bytes.
pub(super) fn header(version: u32, sa_width: u8, payload_len: usize) -> [u8; HEADER_LEN] {
    let mut header = [0u8; HEADER_LEN];
    header[..MAGIC.len()].copy_from_slice(&MAGIC);
    header[8..12].copy_from_slice(&version.to_le_bytes());
    header[12] = sa_width;
    header[13] = LITTLE_ENDIAN;
    header[16..24].copy_from_slice(&(payload_len as u64).to_le_bytes());
    header
}

/// Serialize `store` into an anonymous read-only mapping holding the same
/// header-and-payload image a file would.
pub(super) fn map_anon(store: &TargetStore) -> Result<Mmap> {
    let payload_len = serialize(store, IoWriter::new(io::sink()))?.pos();

    let mut image = MmapOptions::new()
        .len(HEADER_LEN + payload_len)
        .map_anon()?;
    image[..HEADER_LEN].copy_from_slice(&header(FORMAT_VERSION, SA_WIDTH, payload_len));
    serialize(store, Buffer::from(&mut image[HEADER_LEN..]))?;

    Ok(image.make_read_only()?)
}

/// Write `store` into `file`, header last once the payload length is known.
pub(super) fn write_store(file: impl Write + Seek, store: &TargetStore) -> Result<()> {
    let mut file = BufWriter::new(file);
    file.write_all(&[0; HEADER_LEN])?;

    let writer = serialize(store, IoWriter::new(file))?;
    let payload_len = writer.pos();

    let mut file = writer.into_inner();
    file.seek(SeekFrom::Start(0))?;
    file.write_all(&header(FORMAT_VERSION, SA_WIDTH, payload_len))?;
    Ok(file.flush()?)
}

fn serialize<W: rkyv::ser::Writer<rkyv::rancor::BoxedError>>(
    store: &TargetStore,
    writer: W,
) -> Result<W> {
    // Only the writer can fail while serializing a `TargetStore`.
    Ok(rkyv::api::high::to_bytes_in(store, writer).map_err(io::Error::other)?)
}

#[cfg(target_endian = "little")]
#[inline]
pub(super) fn archived_u64_as_native(values: &[rkyv::primitive::ArchivedU64]) -> &[u64] {
    // SAFETY: rkyv's aligned little-endian archived u64 is a transparent-sized,
    // 8-aligned wrapper over native u64 bytes on little-endian targets.
    unsafe { std::slice::from_raw_parts(values.as_ptr().cast::<u64>(), values.len()) }
}

#[cfg(not(target_endian = "little"))]
#[inline]
pub(super) fn archived_u64_as_native(_values: &[rkyv::primitive::ArchivedU64]) -> &[u64] {
    panic!("mmap-backed target indexes currently require a little-endian target")
}

#[cfg(test)]
mod tests {
    use rkyv::ser::writer::IoWriter;
    use tempfile::tempdir;

    use crate::error::Error;
    use crate::index::store::TargetRegistry;
    use crate::types::Base;

    use super::{
        header, serialize, TargetRecord, TargetStore, FORMAT_VERSION, HEADER_LEN, SA_WIDTH,
    };

    /// One target of one base: `A Gap A Gap`.
    fn tiny_store() -> TargetStore {
        TargetStore {
            targets: vec![TargetRecord {
                name: "t1".to_string(),
                len: 1,
            }],
            sequence: vec![1, 0, 1, 0],
            suffix_array: vec![0, 1, 2, 3],
        }
    }

    fn write_image(path: &std::path::Path, version: u32, payload: &[u8]) {
        let mut image = header(version, SA_WIDTH, payload.len()).to_vec();
        image.extend_from_slice(payload);
        fs_err::write(path, &image).unwrap();
    }

    fn write_forged(path: &std::path::Path, version: u32, store: &TargetStore) {
        let payload = rkyv::to_bytes::<rkyv::rancor::Error>(store).unwrap();
        write_image(path, version, payload.as_slice());
    }

    fn open_err(path: &std::path::Path) -> String {
        let Err(err) = TargetRegistry::open(path) else {
            panic!("expected {} to be rejected", path.display());
        };
        format!("{err:#}")
    }

    #[test]
    fn a_valid_index_opens() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        write_forged(&path, FORMAT_VERSION, &tiny_store());

        TargetRegistry::open(&path).expect("tiny_store is a valid index");
    }

    #[test]
    fn open_rejects_file_shorter_than_the_header() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        fs_err::write(&path, vec![0u8; HEADER_LEN - 1]).unwrap();

        assert!(open_err(&path).contains("too small to hold a header"));
    }

    /// `HEADER_LEN` bytes is a complete header, so the size check must not
    /// reject it. The empty payload still fails, but later and for its own reason.
    #[test]
    fn a_header_sized_file_clears_the_size_check() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        write_image(&path, FORMAT_VERSION, &[]);

        let err = open_err(&path);
        assert!(!err.contains("too small to hold a header"), "{err}");
    }

    #[test]
    fn open_rejects_a_foreign_magic() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        fs_err::write(&path, vec![b'x'; HEADER_LEN * 2]).unwrap();

        assert!(open_err(&path).contains("Not a risearch target index"));
    }

    #[test]
    fn open_rejects_older_versions() {
        for version in [1, 2, 3] {
            let dir = tempdir().unwrap();
            let path = dir.path().join("t.idx");
            write_forged(&path, version, &tiny_store());

            let message = open_err(&path);
            assert!(message.contains(&format!("Unsupported target index version {version}")));
            assert!(message.contains("rebuild it with `risearch index"));
        }
    }

    #[test]
    fn open_rejects_a_declared_length_that_does_not_match() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        write_forged(&path, FORMAT_VERSION, &tiny_store());

        let mut bytes = fs_err::read(&path).unwrap();
        bytes.truncate(bytes.len() - 8);
        fs_err::write(&path, &bytes).unwrap();

        assert!(open_err(&path).contains("its header declares"));
    }

    #[test]
    fn open_rejects_a_payload_that_is_not_an_archive() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        write_image(&path, FORMAT_VERSION, b"not-a-valid-rkyv-archive");

        assert!(open_err(&path).contains("Invalid target index archive"));
    }

    /// A file from before the header existed starts with archive bytes, so it is
    /// turned away by the magic check instead of failing bytecheck.
    #[test]
    fn open_rejects_a_headerless_legacy_archive() {
        #[derive(rkyv::Archive, rkyv::Serialize)]
        struct TargetRecordV3 {
            name: String,
            offset: u64,
        }

        #[derive(rkyv::Archive, rkyv::Serialize)]
        struct SuffixIndexV3 {
            sequence: Vec<u8>,
            suffix_array: Vec<u64>,
        }

        #[derive(rkyv::Archive, rkyv::Serialize)]
        struct TargetStoreV3 {
            version: u32,
            targets: Vec<TargetRecordV3>,
            suffix_index: SuffixIndexV3,
        }

        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        let legacy = TargetStoreV3 {
            version: 3,
            targets: vec![TargetRecordV3 {
                name: "ENSG00000139618.15_transcript".to_string(),
                offset: 0,
            }],
            suffix_index: SuffixIndexV3 {
                sequence: vec![1, 0, 1, 0],
                suffix_array: vec![0, 1, 2, 3],
            },
        };
        let bytes = rkyv::to_bytes::<rkyv::rancor::Error>(&legacy).unwrap();
        fs_err::write(&path, bytes.as_slice()).unwrap();

        let message = open_err(&path);
        assert!(message.contains("Not a risearch target index"));
        assert!(message.contains("rebuild it with `risearch index`"));
    }

    #[test]
    fn open_rejects_sequence_sa_length_mismatch() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        let mut store = tiny_store();
        store.suffix_array.pop();
        write_forged(&path, FORMAT_VERSION, &store);

        assert!(open_err(&path).contains("sequence/SA length mismatch"));
    }

    #[test]
    fn open_rejects_invalid_base_rank() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        let mut store = tiny_store();
        store.sequence[0] = Base::U.as_u8() + 1;
        write_forged(&path, FORMAT_VERSION, &store);

        assert!(open_err(&path).contains("invalid base rank"));
    }

    #[test]
    fn open_rejects_missing_separator_gap() {
        // Both halves separately: the mid separator, then the trailing one.
        for corrupt_idx in [1, 3] {
            let dir = tempdir().unwrap();
            let path = dir.path().join("t.idx");
            let mut store = tiny_store();
            store.sequence[corrupt_idx] = Base::A.as_u8();
            write_forged(&path, FORMAT_VERSION, &store);

            assert!(
                open_err(&path).contains("missing a separator Gap"),
                "byte {corrupt_idx} was accepted"
            );
        }
    }

    #[test]
    fn open_rejects_blocks_that_do_not_cover_the_sequence() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        let mut store = tiny_store();
        store.sequence = vec![1, 0, 1, 0, 1, 0];
        store.suffix_array = vec![0, 1, 2, 3, 4, 5];
        write_forged(&path, FORMAT_VERSION, &store);

        assert!(open_err(&path).contains("cover 4 bytes of a 6-byte indexed sequence"));
    }

    #[test]
    fn open_rejects_a_length_that_overruns_the_sequence() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        let mut store = tiny_store();
        store.targets[0].len = 4;
        write_forged(&path, FORMAT_VERSION, &store);

        assert!(open_err(&path).contains("exceeds the indexed sequence"));
    }

    /// Bounds are enforced at read time by the clamp in `base_unchecked`, so
    /// opening must not pay an O(n) scan of the suffix array to prove it.
    #[test]
    fn open_tolerates_out_of_range_suffix_position() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.idx");
        let mut store = tiny_store();
        let past_end = store.sequence.len() as u64;
        store.suffix_array[0] = past_end;
        write_forged(&path, FORMAT_VERSION, &store);

        let registry = TargetRegistry::open(&path).expect("open tolerates a stale position");
        let suffixes = registry.view().suffixes();
        assert_eq!(suffixes.suffix_positions(0..1), &[past_end]);
        // SAFETY: index 0 is inside the suffix array.
        assert_eq!(unsafe { suffixes.base_unchecked(0, 0) }, Base::Gap);
    }

    /// Disk-full and similar failures surface from the writer inside rkyv and
    /// must stay `Error::Io`, which the Python binding raises as `OSError`.
    #[test]
    fn serialization_write_errors_are_io_errors() {
        let mut full = [0u8; 8];
        let Err(err) = serialize(&tiny_store(), IoWriter::new(&mut full[..])) else {
            panic!("an 8-byte buffer cannot hold the archive");
        };
        assert!(matches!(err, Error::Io(_)), "{err}");
    }
}
