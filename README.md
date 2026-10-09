# RIsearch

[![CI](https://github.com/RTH-tools/risearch/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/RTH-tools/risearch/actions/workflows/ci.yml?query=branch%3Amain)
[![Release](https://img.shields.io/github/v/release/RTH-tools/risearch?include_prereleases&sort=semver)](https://github.com/RTH-tools/risearch/releases)
[![License](https://img.shields.io/badge/license-BUSL--1.1-blue)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.88%2B-orange)](https://www.rust-lang.org)
[![Python](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2FRTH-tools%2Frisearch%2Fmain%2Fbindings%2Fpython%2Fpyproject.toml)](bindings/python)

RIsearch predicts RNA and DNA interactions: RNA-RNA, DNA-DNA and RNA-DNA hybrids.
You give it query sequences, such as miRNAs or siRNAs, and a set of targets, such as a transcriptome.
It reports every duplex whose predicted free energy is at or below a threshold, with coordinates, strand, energy and optionally the full alignment.

The targets are indexed once into a suffix array.
A search looks up short complementary seeds in that index and extends each one in both directions with a dynamic-programming alignment scored by nearest-neighbour stacking energies.

Version 3 is a rewrite of [RIsearch2](https://doi.org/10.1093/nar/gkw1325) in Rust, with Python bindings.
It is a pre-release.
Options and output may still change before 3.0.0.

## Installation

Prebuilt binaries for Linux (x86_64, aarch64) and macOS (Apple Silicon, Intel) are attached to every [release](https://github.com/RTH-tools/risearch/releases), each with a SHA-256 checksum.

Or build it with Cargo (Rust 1.88 or newer):

```bash
cargo install --locked risearch --version 3.0.0-alpha.4
```

Add `--features openmp` to let `risearch index` build the suffix array on several threads.
This needs an OpenMP runtime the C compiler can find.
On macOS that means installing `libomp` and making it visible to the compiler and linker.

## Quick start

The repository includes a small example: two human miRNAs, miR-24-3p and miR-876-5p, as queries and the 19 transcripts of RHOC as targets.

```bash
risearch index tests/data/target.fa target.idx
risearch search -q tests/data/query.fa -t target.idx -o hits.tsv
```

`hits.tsv` holds one tab-separated line per interaction:

```text
hsa-miR-24-3p MIMAT0000080	2	22	ENSG00000155366|ENST00000285735	306	323	+	-23.09
hsa-miR-24-3p MIMAT0000080	1	22	ENSG00000155366|ENST00000285735	1374	1399	-	-21.71
hsa-miR-24-3p MIMAT0000080	1	22	ENSG00000155366|ENST00000369642	528	553	-	-21.71
```

The columns are query name, query start and end, target name, target start and end, strand, and free energy in kcal/mol.
Coordinates are 1-based and inclusive, and target coordinates always refer to the target as written in the FASTA file.
The order of hits can change from run to run, so sort the output if you need a stable order.

`risearch index --help` and `risearch search --help` list every option; `risearch search --help` groups them by stage and ends with examples, and `-h` gives a short summary.

## Input

`risearch index` reads target sequences and `risearch search` reads queries, both from FASTA or FASTQ.
Files can be plain or compressed with gzip, bzip2, xz or zstd; the format and the compression are detected from the file contents.
Either input can be `-` to read from stdin.
T and U are treated as the same base, so DNA and RNA input give the same result.
Record IDs must be unique within an input.

Index files record their format version.
If a release changes the format, `search` stops with an error naming the file, and you rebuild it with `risearch index`.
Indexes built by RIsearch2 cannot be read.

## Output formats

Choose a format with `-f/--format`.

| Format | Columns after the 8 standard ones |
| --- | --- |
| `minimal` (default) | none |
| `cigar` | pairing string |
| `bindingsite` | pairing string, aligned target, 5' and 3' target flanks (up to 20 nt each) |
| `detailed` | none, but each hit is preceded by a three-line alignment |

The pairing string has one character per alignment column, following the query 5' to 3':

| Code | Meaning |
| --- | --- |
| `P` | Watson-Crick pair (A-U, C-G) |
| `W` | G-U wobble pair |
| `U` | mismatch |
| `T` | extra base on the target side (gap in the query) |
| `Q` | extra base on the query side (gap in the target) |

In `detailed` output the query is on top (5' to 3'), the target underneath (3' to 5'), and the middle row marks Watson-Crick pairs with `|` and wobble pairs with `:`:

```text
uggcucaguu-----cagcaggaacag
|||||| ||      |||||||  |||
accgagacacccugugucgucc-cguc
hsa-miR-24-3p MIMAT0000080	1	22	ENSG00000155366|ENST00000369642	528	553	-	-21.71
```

Output goes to stdout unless you pass `-o`.
A file name ending in `.gz` or `.zst` is compressed with gzip or zstd.
`--compress` overrides the choice made from the file name.
`--multifile` writes one file per query into the directory given by `-o`.

## Energy parameters

Pick a bundled parameter set with `-P/--params` (default `t04`), or a table of your own with `--params-file`.

| Set | Parameters | Query / target |
| --- | --- | --- |
| `t04` | Turner 2004 | RNA / RNA |
| `slh04` | SantaLucia and Hicks 2004 | DNA / DNA |
| `s95-rna-dna` | Sugimoto 1995 | RNA / DNA |
| `s95-dna-rna` | Sugimoto 1995 | DNA / RNA |

Each set ships tables for 0, 25, 37, 42 and 50 °C, chosen with `-T`.
Temperatures between those are interpolated; temperatures outside 0 to 50 °C are rejected.

`-d` adds an energy penalty in kcal/mol for every nucleotide of the duplex, on both strands and including the seed, and the reported energy includes it.

### Custom energy tables

Pass a tab-separated file to `--params-file`, with this header:

```text
q1	q2	t1	t2	delta_g_kcal_per_mol
```

Each row gives the stacking free energy of two adjacent base pairs, `q1` with `t1` and `q2` with `t2`:

```text
query   5'─ q1 ─ q2 ─ 3'
             |    |
target  3'─ t1 ─ t2 ─ 5'
```

The query dinucleotide is read 5' to 3' and the target dinucleotide 3' to 5', so the target side appears reversed compared with its FASTA sequence.
For example, the row `A C U G` scores the query 5'-AC-3' paired with the target 3'-UG-5', which reads GU in the target's FASTA file.
Bases are `A C G U N`, and `-` stands for a gap.
Energies are in kcal/mol, negative is favourable.

Rows that are missing, or have a value of 20 or more, count as unobserved and score +20 kcal/mol.
Rows of the form `X - X -` or `- X - X` are the helix initiation terms.
The file needs at least one of them, and the initiation offset is derived from them.
A custom table is used as-is at any temperature; setting a temperature with one prints a warning.

## Python

The Python package wraps the same library and returns results as a [Polars](https://pola.rs) DataFrame; see [bindings/python/README.md](bindings/python/README.md).

## Troubleshooting

No hits: raise the energy threshold (`-e -15` instead of `-20`), shorten the seed with `--seed-length`, allow seed mismatches with `--mismatch-max`, or add `--seed-wobble`.

`omp.h file not found` when building with `--features openmp`: the compiler cannot find an OpenMP installation.
Install one for your platform, or build without the feature.

## License

[Business Source License 1.1](LICENSE)
