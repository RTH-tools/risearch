"""Python interface to RIsearch."""

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
    """Build a reusable target index from a FASTA file."""
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
    """Search for RNA-RNA interactions and return one row per hit.

    ``query`` accepts one path or a list of paths. Result order is not
    stable across runs; sort the returned DataFrame when deterministic ordering
    matters. ``matrix`` is a bundled model id or a path to a custom TSV table;
    ``temperature`` defaults to 37 °C and has no effect on a custom table, where
    passing it logs a warning.
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


__all__ = [
    "TargetRegistry",
    "index",
    "search",
]
