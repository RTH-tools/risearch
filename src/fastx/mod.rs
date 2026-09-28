//! FASTA/FASTQ reading and normalization into [`Sequence`].

use needletail::parse_fastx_reader;
use std::collections::HashSet;
use std::io::Read;
use std::path::Path;

use crate::error::{Error, Result};
use crate::seq::Sequence;

/// Read a FASTA/FASTQ file into normalized sequences.
///
/// Record ids must be non-empty and unique. Records that normalize to nothing
/// are dropped; an input where every record does is an error.
pub fn read_sequences(filename: impl AsRef<Path>) -> Result<Vec<(String, Sequence)>> {
    let filename = filename.as_ref();

    let md = fs_err::metadata(filename)?;
    if md.is_dir() {
        return Err(Error::Input(format!(
            "Input path is not a file: {}",
            filename.display()
        )));
    }
    if md.is_file() && md.len() == 0 {
        return Err(Error::Input(format!(
            "No sequences found in input file: {}",
            filename.display()
        )));
    }

    let file = fs_err::File::open(filename)
        .map_err(|err| Error::Input(format!("Failed to open FASTA/FASTQ file: {err}")))?;
    read_sequences_from(file, &filename.display().to_string())
}

/// [`read_sequences`] over an already open `reader`, such as stdin; `source`
/// names it in errors.
pub fn read_sequences_from(
    reader: impl Read + Send,
    source: &str,
) -> Result<Vec<(String, Sequence)>> {
    let mut reader = parse_fastx_reader(reader)
        .map_err(|err| Error::Input(format!("Failed to read FASTA/FASTQ from {source}: {err}")))?;

    let mut sequences = Vec::new();
    let mut seen = HashSet::new();

    while let Some(record) = reader.next() {
        let rec = record.map_err(|err| {
            Error::Input(format!(
                "Failed to parse FASTA/FASTQ record from {source}: {err}"
            ))
        })?;

        let id = String::from_utf8_lossy(rec.id()).into_owned();
        if id.trim().is_empty() {
            return Err(Error::Input(format!(
                "Encountered empty FASTA record id in {source}"
            )));
        }
        if !seen.insert(id.clone()) {
            return Err(Error::Input(format!(
                "Duplicate FASTA record id '{id}' in {source}"
            )));
        }

        if let Some(sequence) = normalize_record(&id, &rec.seq())? {
            sequences.push((id, sequence));
        }
    }

    if sequences.is_empty() {
        return Err(Error::Input(format!(
            "All sequences were empty after normalization in {source}"
        )));
    }

    Ok(sequences)
}

/// Normalizes a raw sequence and handles logging for gaps/N-conversions.
/// Returns `Ok(None)` if the sequence is entirely empty after normalization.
fn normalize_record(id: &str, raw_seq: &[u8]) -> Result<Option<Sequence>> {
    let (sequence, stats) = Sequence::normalize(id, raw_seq)?;

    if sequence.is_empty() {
        log::warn!(
            "Skipping empty sequence after normalization: '{}' (removed_gaps={}, converted_to_n={})",
            id,
            stats.removed_gaps,
            stats.converted_to_n
        );
        return Ok(None);
    }

    if stats.removed_gaps > 0 || stats.converted_to_n > 0 {
        log::debug!(
            "Normalized sequence '{}': removed_gaps={}, converted_to_n={}",
            id,
            stats.removed_gaps,
            stats.converted_to_n
        );
    }

    Ok(Some(sequence))
}
