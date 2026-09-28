"""The public functions, re-exported by the package."""

from __future__ import annotations

import os
from os import PathLike

import polars as pl

from ._native import TargetRegistry
from ._native import build_index as _build_index
from ._native import search as _native_search


def index(
    fasta: str | PathLike[str],
    output: str | PathLike[str],
    *,
    threads: int | None = None,
) -> None:
    """Build a reusable target index from a FASTA file.

    Parameters
    ----------
    fasta
        Target sequences, FASTA or FASTQ, plain or gzip-compressed. T and U
        are the same base, so DNA and RNA input give the same index.
    output
        Where to write the index.
    threads
        Threads for building the index; only builds with the ``openmp``
        feature use it, and the published wheels do not.

    Raises
    ------
    ValueError
        Input that cannot be indexed: unparsable FASTA/FASTQ, empty or
        duplicate IDs, or a file whose sequences are all empty.
    FileNotFoundError, PermissionError, OSError
        A file that is missing, unreadable, or cannot be written.

    Notes
    -----
    An index records its format version. Opening one written by an
    incompatible release raises ``ValueError``; build it again. Building
    releases the GIL, so other Python threads keep running.
    """
    _build_index(fasta, output, threads=threads)


def search(
    query: str | PathLike[str] | list[str | PathLike[str]],
    target: TargetRegistry,
    *,
    seed_length: int | None = None,
    seed_start: int | None = None,
    seed_end: int | None = None,
    energy_threshold: float = -20.0,
    mismatches: int = 0,
    mismatch_prefix: int = 1,
    mismatch_suffix: int = 0,
    seed_wobble: bool = False,
    matrix: str | PathLike[str] = "t04",
    penalty: float = 0.0,
    temperature: int | None = None,
    max_extension: int = 20,
    seed_energy: float = 0.0,
    no_max_prune: bool = False,
    no_dedup: bool = False,
    alignment: bool = True,
    threads: int | None = None,
) -> pl.DataFrame:
    """Search for RNA and DNA interactions and return one row per hit.

    Parameters
    ----------
    query
        One query file or a list of them, FASTA or FASTQ. IDs must be unique
        across all files.
    target
        An index opened with `TargetRegistry.open`.
    seed_length
        Minimum length of a complementary seed; ``None`` means 6, or the whole
        seed interval when ``seed_start`` and ``seed_end`` are set.
    seed_start, seed_end
        Only look for seeds inside these query positions (1-based; ``-1`` is
        the last base). Set both or neither.
    energy_threshold
        Report hits with free energy at or below this, in kcal/mol.
    mismatches
        Mismatches allowed inside a seed.
    mismatch_prefix, mismatch_suffix
        Matches required at the 5' and 3' end of a mismatched seed.
    seed_wobble
        Allow G-U pairs inside seeds.
    matrix
        Energy model id (``t04``, ``slh04``, ``s95-rna-dna``, ``s95-dna-rna``)
        or path to a custom table.
    penalty
        kcal/mol added per extended nucleotide; included in the reported
        energy.
    temperature
        Temperature in °C, 0 to 50; ``None`` means 37. A custom table ignores
        it and logs a warning.
    max_extension
        Extend at most this many nucleotides on each side of the seed; ``-1``
        extends over the whole query.
    seed_energy
        Accepted for compatibility; the search does not use it yet.
    no_max_prune
        Keep seeds that could be extended by one more pair.
    no_dedup
        Report every seed's hit instead of one best hit per duplex.
    alignment
        Fill the ``alignment`` column; ``False`` leaves it null.
    threads
        Worker threads for this search; ``None`` uses ``RAYON_NUM_THREADS``,
        else every core. Must be at least 1.

    Returns
    -------
    polars.DataFrame
        One row per hit; see Notes for the columns.

    Raises
    ------
    ValueError
        Invalid options, including ``threads`` below 1; an empty query list;
        query input that cannot be searched (unparsable FASTA/FASTQ, empty or
        duplicate IDs, files whose sequences are all empty); an energy model
        that cannot be loaded.
    FileNotFoundError, PermissionError, OSError
        A query file that is missing or cannot be read.
    RuntimeError
        The worker threads cannot be started.

    Notes
    -----
    Columns of the result:

    - ``query_idx``, ``target_idx`` (UInt64): position of the query in the input
      and of the target in the index.
    - ``query_name``, ``target_name`` (String): FASTA identifiers.
    - ``q_start``, ``q_end``, ``t_start``, ``t_end`` (UInt64): 0-based, inclusive
      coordinates; target coordinates refer to the target as written.
    - ``strand`` (String): ``+`` or ``-``.
    - ``energy`` (Float64): free energy in kcal/mol.
    - ``alignment`` (String, nullable): pairing string, one character per
      alignment column, query 5' to 3'.

    The pairing string uses ``P`` for a Watson-Crick pair, ``W`` for a G-U
    wobble pair, ``U`` for a mismatch, ``T`` for an extra base on the target
    side and ``Q`` for an extra base on the query side.

    The row order can change between runs; sort when you need a stable order.
    The search releases the GIL, so other Python threads keep running.
    """
    if isinstance(query, (str, os.PathLike)):
        query = [query]

    result = _native_search(
        query,
        target,
        seed_length=seed_length,
        seed_start=seed_start,
        seed_end=seed_end,
        energy_threshold=energy_threshold,
        mismatches=mismatches,
        mismatch_prefix=mismatch_prefix,
        mismatch_suffix=mismatch_suffix,
        seed_wobble=seed_wobble,
        matrix=matrix,
        penalty=penalty,
        temperature=temperature,
        max_extension=max_extension,
        seed_energy=seed_energy,
        no_max_prune=no_max_prune,
        no_dedup=no_dedup,
        alignment=alignment,
        threads=threads,
    )
    return pl.DataFrame(result)
