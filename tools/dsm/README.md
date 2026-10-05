# DSM table generator

Generates dinucleotide stacking energy matrices (DSM) for risearch.
Each table is a 36×36 matrix of nearest-neighbor thermodynamic parameters derived from ViennaRNA duplexfold energies via linear regression.

## Tables

Three parameter sets, five temperatures each (0, 25, 37, 42, 50 °C):

| ID     | Parameters                         | Strand symmetry |
|--------|------------------------------------|-----------------|
| t04    | Turner 2004 (RNA/RNA)              | Yes             |
| slh04  | SantaLucia-Hicks 2004 (DNA/DNA)    | Yes             |
| s95    | Sugimoto 1995 (RNA/DNA hybrid)     | No              |

## How it works

1. **Generate** random query/target duplex pairs (600K total per temperature)
2. **Fold** each pair with ViennaRNA `duplexfold` using the appropriate energy parameters
3. **Count** dinucleotide stacks from the folded structure
4. **Regress** ΔG ~ stack counts via OLS on symmetry-class features
5. **Fill** the 36×36 matrix from regression coefficients
6. **Postprocess**: symmetrize homoduplex matrices, compute initiation offset

## Setup

```bash
uv sync
```

## Usage

Generate all 15 TSV tables:

```bash
uv run python pipeline.py all
```

Each TSV loads directly in risearch via `--params-file path/to/37.tsv`.
`s95` tables load in RNA-query/DNA-target orientation; the transposed form exists only as the bundled `s95-dna-rna`.

Generate the Rust tables and their `mod.rs` for the runtime (writes to `../../src/dsm/tables/`):

```bash
uv run python pipeline.py rust
```

Generate one table:

```bash
uv run python pipeline.py one t04 37
```

## Files

- `dtypes.py` — types: `Base`, `ParamSet`, `DsmId`, `DsmTable`, `DsmData`
- `pipeline.py` — generation, regression, matrix construction, CLI
- `rna_dna_sugimoto1995.par` — Sugimoto 1995 RNA/DNA hybrid parameters for ViennaRNA

## Notes

- Tables are deterministic (seed 19328471) but differ from historical tables due to ViennaRNA version differences and the s95 strand orientation fix
- The s95 orientation follows Sugimoto 1995: query = RNA strand, target = DNA strand
- Generated Rust files are checked into `src/dsm/tables/` and should be regenerated when pipeline parameters change
