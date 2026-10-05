"""Smoke tests for the risearch Python bindings.

Requires the extension to be built first:
    cd bindings/python && uv run --locked maturin develop
"""

import gzip
import inspect
import logging
import os
import subprocess
from pathlib import Path
from typing import Any

import polars as pl
import pytest
from polars.testing import assert_frame_equal

import risearch
import risearch._native as native

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DATA = Path(__file__).resolve().parents[3] / "tests" / "data"
TARGET_FA = DATA / "target.fa"
QUERY_FA = DATA / "query.fa"


@pytest.fixture(scope="session")
def target_index(tmp_path_factory):
    idx = tmp_path_factory.mktemp("idx") / "RHOC.idx"
    risearch.index(TARGET_FA, idx)
    return idx


@pytest.fixture(scope="session")
def store(target_index):
    return risearch.TargetRegistry.open(target_index)


# ---------------------------------------------------------------------------
# index / TargetRegistry.open
# ---------------------------------------------------------------------------


def test_index_creates_file(tmp_path):
    idx = tmp_path / "out.idx"
    risearch.index(TARGET_FA, idx)
    assert idx.exists() and idx.stat().st_size > 0


def test_target_registry_repr_shows_count(store):
    assert "targets=" in repr(store)
    count = int(repr(store).split("targets=")[1].rstrip(")"))
    assert count > 0


# ---------------------------------------------------------------------------
# search — DataFrame shape and schema
# ---------------------------------------------------------------------------

EXPECTED_COLUMNS = {
    "query_idx": pl.UInt64,
    "query_name": pl.String,
    "target_idx": pl.UInt64,
    "target_name": pl.String,
    "q_start": pl.UInt64,
    "q_end": pl.UInt64,
    "t_start": pl.UInt64,
    "t_end": pl.UInt64,
    "strand": pl.String,
    "energy": pl.Float64,
    "alignment": pl.String,
}


def test_search_schema(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    for col, dtype in EXPECTED_COLUMNS.items():
        assert col in df.columns, f"missing column: {col}"
        assert df[col].dtype == dtype, f"{col}: expected {dtype}, got {df[col].dtype}"
    assert df["query_name"].is_not_null().all()
    assert df["target_name"].is_not_null().all()
    assert (df["query_name"].str.len_chars() > 0).all()
    assert (df["target_name"].str.len_chars() > 0).all()


def test_empty_result_has_correct_schema(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-1000.0)
    assert isinstance(df, pl.DataFrame)
    assert len(df) == 0
    for col, dtype in EXPECTED_COLUMNS.items():
        assert col in df.columns
        assert df[col].dtype == dtype


# ---------------------------------------------------------------------------
# search — value correctness
# ---------------------------------------------------------------------------


def test_all_energies_below_threshold(store):
    threshold = -15.0
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=threshold)
    assert (df["energy"] <= threshold).all()


def test_search_value_invariants(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert df["strand"].is_in(["+", "-"]).all()
    assert (df["q_start"] <= df["q_end"]).all()
    assert (df["t_start"] <= df["t_end"]).all()


# ---------------------------------------------------------------------------
# search — filtering behaviour
# ---------------------------------------------------------------------------


def test_stricter_cutoff_returns_exactly_the_eligible_rows(store):
    options = dict(seed_length=8, alignment=False)
    loose = risearch.search(QUERY_FA, store, energy_threshold=-5.0, **options)
    energies = sorted(set(loose["energy"]))
    assert len(energies) >= 2, "fixture must exercise both acceptance and rejection"
    cutoff = energies[len(energies) // 2]
    expected = loose.filter(pl.col("energy") <= cutoff)
    assert 0 < len(expected) < len(loose)
    actual = risearch.search(QUERY_FA, store, energy_threshold=cutoff, **options)
    assert_frame_equal(expected.sort(SORT_COLS), actual.sort(SORT_COLS))


def test_seed_length_and_wobble_select_the_specified_spans(tmp_path):
    """AGG/UUU has no two-base WC seed; wobble permits all five spans of length >= 2."""
    query, target, index = (tmp_path / name for name in ("q.fa", "t.fa", "t.idx"))
    query.write_text(">q\nAGG\n")
    target.write_text(">t\nUUU\n")
    risearch.index(target, index)
    store = risearch.TargetRegistry.open(index)
    all_spans = [
        (0, 1, 0, 1, "+"),
        (0, 1, 1, 2, "+"),
        (0, 2, 0, 2, "+"),
        (1, 2, 0, 1, "+"),
        (1, 2, 1, 2, "+"),
    ]
    for length, wobble, expected in [
        (2, False, []),
        (2, True, all_spans),
        (3, True, [(0, 2, 0, 2, "+")]),
    ]:
        actual = risearch.search(
            query,
            store,
            seed_length=length,
            seed_wobble=wobble,
            max_extension=0,
            no_max_prune=True,
            no_dedup=True,
            energy_threshold=100.0,
        )
        assert (
            sorted(
                actual.select("q_start", "q_end", "t_start", "t_end", "strand").rows()
            )
            == expected
        )


# ---------------------------------------------------------------------------
# search — multi-file input
# ---------------------------------------------------------------------------


def test_single_path_and_list_equivalent(store):
    df_str = risearch.search(
        str(QUERY_FA), store, seed_length=8, energy_threshold=-10.0
    )
    df_list = risearch.search([QUERY_FA], store, seed_length=8, energy_threshold=-10.0)
    assert len(df_str) > 0
    assert_frame_equal(df_str.sort(SORT_COLS), df_list.sort(SORT_COLS))


def test_split_files_match_full_file(store, tmp_path):
    lines = QUERY_FA.read_text().splitlines()
    records, current = [], []
    for line in lines:
        if line.startswith(">") and current:
            records.append("\n".join(current))
            current = []
        current.append(line)
    if current:
        records.append("\n".join(current))

    mid = len(records) // 2
    f1, f2 = tmp_path / "part1.fa", tmp_path / "part2.fa"
    f1.write_text("\n".join(records[:mid]) + "\n")
    f2.write_text("\n".join(records[mid:]) + "\n")

    df_split = risearch.search([f1, f2], store, seed_length=8, energy_threshold=-10.0)
    df_full = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert len(df_full) > 0
    assert_frame_equal(df_split.sort(SORT_COLS), df_full.sort(SORT_COLS))


@pytest.mark.parametrize(
    "fastq,compressed", [(False, True), (True, False), (True, True)]
)
def test_fastq_and_gzip_preserve_fasta_search_results(tmp_path, fastq, compressed):
    """Exercise both input boundaries; FASTQ quality values do not change bases."""

    def write_records(name, records, *, fastq=False, compressed=False):
        text = "".join(
            f"@{name}\n{seq}\n+\n{'!' * len(seq)}\n" if fastq else f">{name}\n{seq}\n"
            for name, seq in records
        )
        path = tmp_path / name
        path.write_bytes(gzip.compress(text.encode()) if compressed else text.encode())
        return path

    queries = [("q-a", "AAAA"), ("q-c", "CCCC")]
    targets = [("t-u", "UUUUUU"), ("t-g", "GGGGGG")]
    plain_query = write_records("query.fa", queries)
    plain_target = write_records("target.fa", targets)
    suffix = (".fq" if fastq else ".fa") + (".gz" if compressed else "")
    variant_query = write_records(
        "query-variant" + suffix, queries, fastq=fastq, compressed=compressed
    )
    variant_target = write_records(
        "target-variant" + suffix, targets, fastq=fastq, compressed=compressed
    )
    plain_index, variant_index = tmp_path / "plain.idx", tmp_path / "variant.idx"
    risearch.index(plain_target, plain_index)
    risearch.index(variant_target, variant_index)
    options = dict(seed_length=2, max_extension=2, energy_threshold=100.0)
    expected = risearch.search(
        plain_query, risearch.TargetRegistry.open(plain_index), **options
    )
    actual = risearch.search(
        variant_query, risearch.TargetRegistry.open(variant_index), **options
    )
    assert set(expected["query_name"]) == {"q-a", "q-c"}
    assert_frame_equal(expected.sort(SORT_COLS), actual.sort(SORT_COLS))


# ---------------------------------------------------------------------------
# search — alignment and threads
# ---------------------------------------------------------------------------

SORT_COLS = [c for c in EXPECTED_COLUMNS if c != "alignment"]


def canonical_rows(df: pl.DataFrame) -> list[tuple]:
    """Return a stable named multiset, retaining coordinates and alignment."""
    return sorted(
        tuple(
            row[col] if col != "energy" else round(float(row[col]), 4)
            for col in EXPECTED_COLUMNS
        )
        for row in df.to_dicts()
    )


def cli_path() -> Path:
    """Locate the Cargo-built CLI used by the cross-frontend fixture."""
    configured = os.environ.get("CARGO_BIN_EXE_risearch")  # noqa: SIM112
    candidates = [Path(configured)] if configured else []
    candidates.append(DATA.parents[1] / "target" / "debug" / "risearch")
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate
    pytest.fail(
        "Python/CLI parity test needs a built risearch binary; set "
        "CARGO_BIN_EXE_risearch or run `cargo build --bin risearch` first"
    )


def test_alignment_toggle(store):
    with_aln = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    without = risearch.search(
        QUERY_FA, store, seed_length=8, energy_threshold=-10.0, alignment=False
    )
    assert without["alignment"].is_null().all()
    for col, dtype in EXPECTED_COLUMNS.items():
        assert col in without.columns
        assert without[col].dtype == dtype
    assert_frame_equal(
        with_aln.drop("alignment").sort(SORT_COLS),
        without.drop("alignment").sort(SORT_COLS),
    )


@pytest.mark.parametrize("threads", [1, 2])
def test_threads_does_not_change_results(store, threads):
    default = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    pinned = risearch.search(
        QUERY_FA, store, seed_length=8, energy_threshold=-10.0, threads=threads
    )
    columns = list(EXPECTED_COLUMNS)
    assert_frame_equal(default.sort(columns), pinned.sort(columns))


def test_index_accepts_threads(tmp_path):
    idx = tmp_path / "threaded.idx"
    risearch.index(TARGET_FA, idx, threads=2)
    assert idx.exists() and idx.stat().st_size > 0


@pytest.mark.parametrize("threads", [0, -1])
def test_threads_below_one_are_rejected(tmp_path, store, threads):
    with pytest.raises(ValueError, match="threads must be >= 1"):
        risearch.index(TARGET_FA, tmp_path / "rejected.idx", threads=threads)
    with pytest.raises(ValueError, match="threads must be >= 1"):
        risearch.search(QUERY_FA, store, threads=threads)


def test_fasta_case_wrapping_and_tu_spelling_preserve_hits(tmp_path):
    query_plain = tmp_path / "query-plain.fa"
    query_variant = tmp_path / "query-variant.fa"
    target_plain = tmp_path / "target-plain.fa"
    target_variant = tmp_path / "target-variant.fa"
    query_plain.write_text(">q\nAAAA\n")
    query_variant.write_text(">q\naA\naA\n")
    target_plain.write_text(">t\nUUUUUU\n")
    target_variant.write_text(">t\nttTT\ntt\n")

    plain_index = tmp_path / "plain.idx"
    variant_index = tmp_path / "variant.idx"
    risearch.index(target_plain, plain_index, threads=1)
    risearch.index(target_variant, variant_index, threads=1)
    plain = risearch.search(
        query_plain,
        risearch.TargetRegistry.open(plain_index),
        seed_length=2,
        max_extension=2,
        energy_threshold=100.0,
    )
    variant = risearch.search(
        query_variant,
        risearch.TargetRegistry.open(variant_index),
        seed_length=2,
        max_extension=2,
        energy_threshold=100.0,
    )
    assert len(plain) > 0, "normalization fixture must produce hits"
    assert canonical_rows(plain) == canonical_rows(variant)


# ---------------------------------------------------------------------------
# search — invalid arguments
# ---------------------------------------------------------------------------


def test_seed_start_without_end_raises(store):
    with pytest.raises(ValueError, match="seed_start"):
        risearch.search(QUERY_FA, store, seed_start=2)


def test_python_search_matches_cli_minimal_fixture(tmp_path):
    """Cross-check the Python Arrow surface against CLI TSV at 2dp energy."""
    target = tmp_path / "targets.fa"
    query = tmp_path / "queries.fa"
    index = tmp_path / "targets.idx"
    target.write_text(">t-a\nUUUUUU\n>t-c\nGGGGGG\n")
    query.write_text(">q-a\nAAAA\n>q-c\nCCCC\n")

    cli = cli_path()
    subprocess.run(
        [str(cli), "index", str(target), str(index)],
        check=True,
        capture_output=True,
        text=True,
    )
    completed = subprocess.run(
        [
            str(cli),
            "search",
            "-q",
            str(query),
            "-t",
            str(index),
            "-o",
            "-",
            "--seed-length",
            "2",
            "-l",
            "2",
            "-e",
            "100.0",
            "--format",
            "minimal",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    cli_rows = []
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        assert len(fields) == 8
        cli_rows.append(
            (
                fields[0],
                fields[3],
                int(fields[1]) - 1,
                int(fields[2]) - 1,
                int(fields[4]) - 1,
                int(fields[5]) - 1,
                fields[6],
                round(float(fields[7]), 2),
            )
        )

    python = risearch.search(
        query,
        risearch.TargetRegistry.open(index),
        seed_length=2,
        max_extension=2,
        energy_threshold=100.0,
        alignment=False,
    )
    python_rows = [
        (
            row["query_name"],
            row["target_name"],
            row["q_start"],
            row["q_end"],
            row["t_start"],
            row["t_end"],
            row["strand"],
            round(float(row["energy"]), 2),
        )
        for row in python.to_dicts()
    ]
    assert cli_rows
    assert sorted(cli_rows) == sorted(python_rows)


def test_invalid_matrix_raises(store, caplog):
    """An unknown id fails validation without first being warned about as a table."""
    with (
        caplog.at_level(logging.WARNING, logger="risearch"),
        pytest.raises(ValueError, match="unknown energy parameter set"),
    ):
        risearch.search(QUERY_FA, store, matrix="t05", temperature=50)
    assert not [r for r in caplog.records if "temperature" in r.getMessage()]


def test_excessive_max_extension_raises(store):
    with pytest.raises(ValueError, match="max extension"):
        risearch.search(QUERY_FA, store, max_extension=257)


# ---------------------------------------------------------------------------
# Public/private package boundary
# ---------------------------------------------------------------------------


def test_public_objects_identify_as_risearch():
    assert risearch.TargetRegistry.__module__ == "risearch"
    assert risearch.index.__name__ == "index"


def test_public_and_native_search_params_match():
    """The wrapper must forward every native keyword, in the same order."""
    public = inspect.signature(risearch.search)
    private = inspect.signature(native.search)
    assert list(public.parameters) == list(private.parameters)


def test_native_search_declares_no_defaults():
    """Defaults belong to the wrapper alone, so the native side must have none."""
    private = inspect.signature(native.search)
    defaulted = [
        name
        for name, p in private.parameters.items()
        if p.default is not inspect.Parameter.empty
    ]
    assert defaulted == []


def test_wrapper_defaults_match_rust_config_defaults():
    """The wrapper's defaults are the Rust config's defaults, not a copy of them."""
    public = inspect.signature(risearch.search).parameters
    exempt = {"query", "target", "threads", "alignment"}
    assert {n: p.default for n, p in public.items() if n not in exempt} == (
        native._default_options()
    )


def test_documented_defaults_are_pinned():
    """Absolute pins: the relative check above moves with the Rust config."""
    public = inspect.signature(risearch.search).parameters
    assert public["energy_threshold"].default == -20.0
    assert public["penalty"].default == 0.0
    assert public["seed_wobble"].default is False
    assert public["max_extension"].default == 20
    assert public["temperature"].default is None
    assert public["matrix"].default == "t04"


def test_native_search_is_callable_with_every_kwarg(store):
    """Keeps a runtime call on the native entry point the stub documents."""
    options: dict[str, Any] = {
        **native._default_options(),
        "seed_length": 8,
        "alignment": True,
    }
    result = native.search([QUERY_FA], store, **options, threads=None)
    assert pl.DataFrame(result).height > 0


def test_explicit_temperature_with_custom_table_warns(store, tmp_path, caplog):
    """`temperature` cannot apply to a user TSV; saying it anyway is flagged once."""
    tsv = tmp_path / "custom.tsv"
    tsv.write_text(
        "q1\tq2\tt1\tt2\tdelta_g_kcal_per_mol\n"
        "A\t-\tU\t-\t3.0\n"
        "-\tA\t-\tU\t3.0\n"
        "G\tG\tC\tC\t-3.3\n"
    )
    # pyo3-log caches the "risearch" level at first use; no test may raise it.
    with caplog.at_level(logging.WARNING, logger="risearch"):
        risearch.search(QUERY_FA, store, matrix=tsv)
        assert not [r for r in caplog.records if "temperature" in r.getMessage()]
        risearch.search(QUERY_FA, store, matrix=tsv, temperature=50)
        assert len([r for r in caplog.records if "temperature" in r.getMessage()]) == 1


# ---------------------------------------------------------------------------
# Error mapping
# ---------------------------------------------------------------------------


def test_missing_index_raises_file_not_found(tmp_path):
    """An absent index is an OS-level miss, not a risearch-specific failure."""
    with pytest.raises(FileNotFoundError):
        risearch.TargetRegistry.open(tmp_path / "absent.idx")


def test_corrupt_index_raises_value_error(tmp_path):
    """A file that is not an index is rejected by the header, before the archive."""
    bogus = tmp_path / "bogus.idx"
    bogus.write_bytes(b"definitely not a risearch index")
    with pytest.raises(ValueError):
        risearch.TargetRegistry.open(bogus)


def test_missing_query_file_raises_file_not_found(store):
    """A query path that does not exist surfaces as FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        risearch.search("no-such-query.fa", store)


def test_empty_query_list_raises_value_error(store):
    """An empty query list is an error, not an empty result."""
    with pytest.raises(ValueError):
        risearch.search([], store)


@pytest.mark.parametrize(
    "contents,reason",
    [
        (">duplicate\nAAAA\n>duplicate\nCCCC\n", "Duplicate"),
        (">empty\n----\n", "empty"),
        ("@query\nAAAA\n+\n!!\n", "parse"),
    ],
)
@pytest.mark.parametrize("operation", ["index", "search"])
def test_invalid_sequence_input_raises_value_error(
    store, tmp_path, contents, reason, operation
):
    """Bad biological input raises ValueError at both public entry points."""
    source = tmp_path / "invalid.fastx"
    source.write_text(contents)
    with pytest.raises(ValueError, match=reason):
        if operation == "index":
            risearch.index(source, tmp_path / "invalid.idx")
        else:
            risearch.search(source, store, seed_length=2)
