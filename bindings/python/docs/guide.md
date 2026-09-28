# Guide

## Energy models

`search(matrix=...)` picks the nearest-neighbour energy model that scores each duplex.

| `matrix` | Parameters | Query / target |
| --- | --- | --- |
| `t04` | Turner 2004 | RNA / RNA |
| `slh04` | SantaLucia and Hicks 2004 | DNA / DNA |
| `s95-rna-dna` | Sugimoto 1995 | RNA / DNA |
| `s95-dna-rna` | Sugimoto 1995 | DNA / RNA |

Each model ships tables for 0, 25, 37, 42 and 50 °C.
Temperatures between those are interpolated; temperatures outside 0 to 50 °C are rejected.

### Custom energy tables

`matrix` also accepts a path to a tab-separated file with this header:

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
Energies are in kcal/mol; negative is favourable.

Rows that are missing, or have a value of 20 or more, count as unobserved and score +20 kcal/mol.
Rows of the form `X - X -` or `- X - X` are the helix initiation terms; the file needs at least one of them, and the initiation offset is derived from them.
A custom table is used as-is at any temperature; passing `temperature` with one logs a warning.

## Logging

Log records go to Python's `logging` under the `risearch` logger:

```python
import logging

logging.basicConfig(level=logging.INFO)
```

Silence them with `logging.getLogger("risearch").setLevel(logging.ERROR)`.
Configure logging before the first call, since the level is read once per module on first use.
Release builds log at `INFO` and above.
