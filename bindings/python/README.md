# risearch Python bindings

## Quick start

From the repo root:

```bash
cd bindings/python
uv sync --locked
uv run --locked maturin develop
uv run --locked python
```

`uv` manages the locked development environment; Maturin builds and installs
the local extension into that environment. Run `maturin develop` again after
changing Rust binding code.

Import the public package as `risearch`. The compiled module is installed as
`risearch._native` and is a private implementation detail.

This Python package does not install the `risearch` command-line program. From
the repository root, install only the CLI with
`cargo install --locked --path .`.

## Legacy CPU / Rosetta install

For older x86-64 CPUs or x86-64 Python running under Rosetta on Apple Silicon,
install the compatibility Polars runtime through the `lts-cpu` extra:

```bash
uv sync --locked --extra lts-cpu
uv run --locked maturin develop
```

For published wheels, the equivalent pip form is:

```bash
pip install "risearch[lts-cpu]"
```

This extra uses Polars' `rtcompat` runtime, so the Python module is still
imported as `polars`.

## What the binding exposes

The Python package is very small:

- `risearch.index(fasta, output, threads=None)` builds a binary target index
  (`threads` needs the `openmp` build feature; wheels are single-threaded here)
- `risearch.TargetRegistry.open(path)` opens that index for reuse
- `risearch.search(query, target, **kwargs)` runs the search and returns a Polars `DataFrame`

The result schema is:

| Column | Polars type | Meaning |
| --- | --- | --- |
| `query_idx`, `target_idx` | `UInt64` | Registry positions |
| `query_name`, `target_name` | `String` | FASTA identifiers |
| `q_start`, `q_end`, `t_start`, `t_end` | `UInt64` | Zero-based, inclusive coordinates |
| `strand` | `String` | `+` or `-` |
| `energy` | `Float64` | Free energy in kcal/mol |
| `alignment` | `String` (nullable) | Pairing fingerprint |

The `search()` kwargs map to the canonical Rust-facing options:

- `seed_length`, `seed_start`, `seed_end`
- `mismatches`, `mismatch_prefix`, `mismatch_suffix`
- `seed_wobble` (off by default; set `True` to allow G-U pairs in seeds)
- `matrix` (bundled model id or path to a custom TSV table), `penalty`, `temperature`
  (defaults to 37 °C; has no effect on a custom table, and a warning is logged if both are given)
- `max_extension`
- `energy_threshold`, `seed_energy`, `no_max_prune`, `no_dedup`
- `alignment` — set `False` to skip DP traceback; the `alignment` column becomes all-null
- `threads` — rayon worker width for query parsing and search; defaults to rayon's
  choice, which honours `RAYON_NUM_THREADS`. Must be at least 1; `0` raises `ValueError`.

## Logging

Rust log records are forwarded into Python's `logging` under the `risearch`
logger, so the application decides where they go:

```python
import logging

logging.basicConfig(level=logging.INFO)
```

Unconfigured, Python's own fallback prints warnings to stderr. Silence the
library with `logging.getLogger("risearch").setLevel(logging.ERROR)`.
Configure logging before the first call: pyo3-log caches each module's
effective level on first use, so later `setLevel` calls may not take effect.
`TRACE` never reaches Python, and release wheels also drop `DEBUG`.

## Where defaults live

`risearch.search()` declares every default; the compiled `_native.search`
declares none. `_native._default_options()` reports the Rust config defaults
and the test suite holds the wrapper to them. This covers the Python API only
— the CLI resolves some options its own way, so the two can still differ
(`alignment` is one: the CLI derives it from the output format).

## Smoke test

Run this from `bindings/python/`:

```bash
uv run --locked python - <<'PY'
from pathlib import Path
import tempfile
import risearch

root = Path.cwd().parents[1]
target_fa = root / "tests" / "data" / "target.fa"
query_fa = root / "tests" / "data" / "query.fa"

with tempfile.TemporaryDirectory() as tmp:
    idx = Path(tmp) / "RHOC.idx"
    risearch.index(target_fa, idx)
    target = risearch.TargetRegistry.open(idx)
    df = risearch.search(query_fa, target, seed_length=8, energy_threshold=-10.0)
    print(df.shape)
    print(df.columns)
PY
```

## Development check

Rust-side build check:

```bash
cargo test -p risearch-python
```

Build the extension and run the Python suite:

```bash
uv sync --locked
uv run --locked maturin develop
uv run --locked pytest -q
```

Build a wheel directly:

```bash
uv run --python 3.10 --locked maturin build --out ../../dist
```

Build the source distribution and verify it by rebuilding a wheel from the
unpacked archive:

```bash
uv run --python 3.10 --locked maturin build --sdist --out ../../dist
```
