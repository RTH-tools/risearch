# RIsearch for Python

[![PyPI](https://img.shields.io/pypi/v/risearch)](https://pypi.org/project/risearch/)
[![Python](https://img.shields.io/pypi/pyversions/risearch)](https://pypi.org/project/risearch/)
[![CI](https://github.com/saiden89/risearch/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/saiden89/risearch/actions/workflows/ci.yml?query=branch%3Amain)
[![License](https://img.shields.io/badge/license-BUSL--1.1-blue)](https://github.com/saiden89/risearch/blob/main/LICENSE)
[![Docs](https://img.shields.io/badge/docs-online-blue)](https://saiden89.github.io/risearch/)

RIsearch predicts RNA and DNA interactions: RNA-RNA, DNA-DNA and RNA-DNA hybrids.
You give it query sequences, such as miRNAs or siRNAs, and a set of targets, such as a transcriptome.
It reports every duplex whose predicted free energy is at or below a threshold, with coordinates, strand, energy and the pairing pattern, as a [Polars](https://pola.rs) DataFrame.

This package wraps the same Rust library as the [`risearch` command-line tool](https://github.com/saiden89/risearch), a rewrite of [RIsearch2](https://doi.org/10.1093/nar/gkw1325).
It is a pre-release: options and output may still change before 3.0.0.

## Installation

```bash
pip install risearch
```

or `uv add risearch`.
Wheels are built for Linux x86_64 (manylinux) and macOS on Apple Silicon, for every supported Python version.
On other platforms pip builds from the source distribution, which needs Rust 1.88 or newer.

On older x86-64 CPUs, or x86-64 Python running under Rosetta, install the compatibility Polars runtime; the module is still imported as `polars`:

```bash
pip install "risearch[lts-cpu]"
```

The package does not install the `risearch` command-line program; its [releases](https://github.com/saiden89/risearch/releases) ship prebuilt binaries.

## Quick start

```python
import risearch

risearch.index("targets.fa", "targets.idx")  # build once
targets = risearch.TargetRegistry.open("targets.idx")  # reuse for every search

hits = risearch.search("mirnas.fa", targets)
```

`hits` has one row per interaction.
With the example in the repository, two human miRNAs against the 19 transcripts of RHOC in [`tests/data`](https://github.com/saiden89/risearch/tree/main/tests/data):

```python
>>> hits.sort("energy", "target_name", "t_start").select(
...     "query_name", "target_name", "t_start", "t_end", "strand", "energy"
... ).head(3)
shape: (3, 6)
┌────────────────────────────┬─────────────────────────────────┬─────────┬───────┬────────┬──────────┐
│ query_name                 ┆ target_name                     ┆ t_start ┆ t_end ┆ strand ┆ energy   │
│ ---                        ┆ ---                             ┆ ---     ┆ ---   ┆ ---    ┆ ---      │
│ str                        ┆ str                             ┆ u64     ┆ u64   ┆ str    ┆ f64      │
╞════════════════════════════╪═════════════════════════════════╪═════════╪═══════╪════════╪══════════╡
│ hsa-miR-24-3p MIMAT0000080 ┆ ENSG00000155366|ENST00000527563 ┆ 481     ┆ 505   ┆ +      ┆ -25.3257 │
│ hsa-miR-24-3p MIMAT0000080 ┆ ENSG00000155366|ENST00000285735 ┆ 305     ┆ 322   ┆ +      ┆ -23.0914 │
│ hsa-miR-24-3p MIMAT0000080 ┆ ENSG00000155366|ENST00000285735 ┆ 1373    ┆ 1398  ┆ -      ┆ -21.7137 │
└────────────────────────────┴─────────────────────────────────┴─────────┴───────┴────────┴──────────┘
```

Coordinates are 0-based and inclusive, unlike the 1-based CLI output.
The row order can change from run to run, so sort when you need a stable order.

The [documentation](https://saiden89.github.io/risearch/) covers every function and option, the result columns, the energy models and logging.

## Citation

If you use RIsearch, please cite:

> Roncelli S, Favaro L, Anthon C, Gorodkin J.
> RIsearch and siOFF: An integrated, high-performance framework for RNA-RNA interaction and siRNA off-target prediction.
> *Bioinformatics*.

## License

[Business Source License 1.1](https://github.com/saiden89/risearch/blob/main/LICENSE).
