# Developing the Python bindings

Run everything from `bindings/python`.
`uv` manages the locked development environment; Maturin builds the extension and installs it into it.

```bash
uv sync --locked
uv run --locked maturin develop --generate-stubs
```

Run `maturin develop --generate-stubs` again after changing Rust code.
Import the public package as `risearch`; the compiled module is installed as `risearch._native` and is a private implementation detail.

`--generate-stubs` writes the type stub `risearch/_native.pyi` from the compiled module.
It is not committed; ty and pyrefly, including the hk hooks, need it on disk, and the published wheels generate their own.

For older x86-64 CPUs or x86-64 Python under Rosetta, sync with the compatibility Polars runtime instead:

```bash
uv sync --locked --extra lts-cpu
uv run --locked maturin develop --generate-stubs
```

## Tests

```bash
cargo test -p risearch-python
uv run --locked pytest -q
```

The Python/CLI parity test needs a built CLI: run `cargo build --bin risearch` first, or set `CARGO_BIN_EXE_risearch`.

## Where defaults live

`risearch.search()` declares every default; the compiled `_native.search` declares none.
`_native._default_options()` reports the Rust config defaults and the test suite holds the wrapper to them.
This covers the Python API only: the CLI resolves some options its own way, so the two can still differ (`alignment` is one: the CLI derives it from the output format).

## Wheels and the source distribution

Build a wheel:

```bash
uv run --python 3.10 --locked maturin build --out ../../dist
```

Build the source distribution and verify it by rebuilding a wheel from the unpacked archive:

```bash
uv run --python 3.10 --locked maturin build --sdist --out ../../dist
```
