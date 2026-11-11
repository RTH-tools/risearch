#!/usr/bin/env python3

import random
from dataclasses import replace
from functools import cache
from itertools import product
from multiprocessing import Pool
from pathlib import Path
from typing import Annotated

import numpy as np
import RNA
import typer
from loguru import logger
from numpy.typing import NDArray
from sklearn.linear_model import LinearRegression

from dtypes import (
    ALL_DSM_IDS,
    COMPLEMENT_TABLES,
    DSM_N,
    RNA_ALPHABET,
    SEED,
    STRATEGIES,
    Base,
    ComplementStrategy,
    DsmData,
    DsmId,
    DsmTable,
    Duplex,
    ParamSet,
    StackIndex,
    StructureChar,
    TrainingStrategy,
)


def reverse_complement(
    seq: str, strategy: ComplementStrategy = ComplementStrategy.STRICT
) -> str:
    """Reverse complement a sequence. WOBBLE strategy allows G-U pairing."""
    table = COMPLEMENT_TABLES[strategy]
    return "".join(random.choice(table[b]) for b in reversed(seq))


def random_seq(length: int) -> str:
    return "".join(random.choices(RNA_ALPHABET, k=length))


_MIN_STACK_COUNT = 4
_FLANK_OVERLAP = 4
_FLANK_MIN = 10
_FLANK_MAX = 20
_G = Base.Gap


def _val(s: StackIndex) -> tuple[int, ...]:
    return tuple(6 if b is _G else b.value for b in s)


def _regression_group(s: StackIndex, *, strand_symmetric: bool) -> object:
    """Group stacks that share one regression coefficient.

    Priority: strand canon → init pairs → WC-anchored → gap-diagonal → gap-normalized → raw.
    """
    if strand_symmetric:
        s = min(s, (s[3], s[2], s[1], s[0]), key=_val)

    i, j, k, l = s
    cross = ((i, k), (j, l))
    wc = [p for p in cross if p[0].pairs_with(p[1])]
    gap = [p for p in cross if p == (_G, _G)]

    if len(wc) == 1 and len(gap) == 1:
        p = wc[0]
        return ("init", min(p, p[::-1], key=lambda t: (t[0].value, t[1].value)))

    if len(wc) == 1:
        other = cross[1] if wc[0] is cross[0] else cross[0]
        return (wc[0], other)

    if s.count(_G) == 1:
        diag = [((j, k), (i, l)), ((i, l), (j, k))]
        for pair, ctx in diag:
            if pair[0].pairs_with(pair[1]):
                if strand_symmetric and ctx[0] is _G:
                    ctx = ctx[::-1]
                return (pair, ctx)
        if j is _G:
            s = (j, i, k, l)
        elif l is _G:
            s = (i, j, l, k)

    return s


def _build_table(strand_symmetric: bool) -> dict[StackIndex, StackIndex]:
    """Precompute StackIndex → canonical representative for all 1296 stacks."""
    table: dict[StackIndex, StackIndex] = {}
    groups: dict[object, StackIndex] = {}
    for stack in product(Base, Base, Base, Base):
        rep = groups.setdefault(
            _regression_group(stack, strand_symmetric=strand_symmetric), stack
        )
        table[stack] = rep
    return table


_CLASS_TABLE = {sym: _build_table(sym) for sym in (True, False)}


def symmetry_class(key: StackIndex, *, strand_symmetric: bool = True) -> StackIndex:
    """Map a stack to its canonical representative for regression feature pooling."""
    return _CLASS_TABLE[strand_symmetric][key]


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------


def generate_duplex_pairs(
    strategy: TrainingStrategy, temperature: float
) -> list[Duplex]:
    OV = _FLANK_OVERLAP
    flank = strategy.flank
    result: list[Duplex] = []
    for x in range(1, strategy.n_samples + 1):
        left = random_seq(random.randint(_FLANK_MIN, _FLANK_MAX))
        right = random_seq(random.randint(_FLANK_MIN, _FLANK_MAX))
        core = random_seq(strategy.core_size.sample())
        core_rc = (
            reverse_complement(core, strategy.complement)
            if strategy.complement is not None
            else random_seq(strategy.core_size.sample())
        )

        query = left + core + right
        target = (
            reverse_complement(right[-OV:], flank)
            + reverse_complement(right[:-OV])
            + core_rc
            + reverse_complement(left[OV:])
            + reverse_complement(left[:OV], flank)
        )

        result.append(Duplex(query, target, temperature, x))
        result.append(Duplex(target, query, temperature, x))
    return result


_HYBRID_PARAMS = Path(__file__).with_name("rna_dna_sugimoto1995.par")


def _init_worker(dsm_id: DsmId, temperature: float) -> None:
    """Per-process ViennaRNA setup: load the right energy parameter set."""
    RNA.cvar.temperature = temperature
    match dsm_id.params:
        case ParamSet.HYBRID:
            RNA.params_load_from_string(_HYBRID_PARAMS.read_text())
            RNA.update_fold_params()
        case ParamSet.DNA:
            RNA.params_load_DNA_Mathews2004()
        case ParamSet.RNA:
            RNA.params_load_RNA_Turner2004()


_FC = Base.from_char


def score_duplex(pair: Duplex) -> Duplex:
    """Fold with ViennaRNA duplexfold, count dinucleotide stacks on the fly."""
    query, target = pair.query, pair.target
    duplex = RNA.duplexfold(query, target)

    O, C, U = StructureChar.OPEN, StructureChar.CLOSE, StructureChar.UNPAIRED

    raw_q, sep, raw_t = duplex.structure.partition("&")
    if not sep or len(raw_q) != len(query) or len(raw_t) != len(target):
        return pair
    qs = [StructureChar(c) for c in raw_q]
    ts = [StructureChar(c) for c in raw_t]

    if qs[-1] is U and ts[0] is U:
        qs[-1], ts[0] = O, C
    if qs[0] is U and ts[-1] is U:
        qs[0], ts[-1] = O, C
    if U in (qs[0], qs[-1], ts[0], ts[-1]):
        return pair

    stacks: dict[StackIndex, int] = {}
    prev_q, prev_t = _G, _G
    qi, tj = 0, len(target) - 1
    while qi < len(query) and tj >= 0:
        match qs[qi], ts[tj]:
            case (StructureChar.UNPAIRED, StructureChar.CLOSE):
                cur_q, cur_t = _FC(query[qi]), _G
                qi += 1
            case (StructureChar.OPEN, StructureChar.UNPAIRED):
                cur_q, cur_t = _G, _FC(target[tj])
                tj -= 1
            case _:
                cur_q, cur_t = _FC(query[qi]), _FC(target[tj])
                qi += 1
                tj -= 1
        idx: StackIndex = (prev_q, cur_q, prev_t, cur_t)
        stacks[idx] = stacks.get(idx, 0) + 1
        prev_q, prev_t = cur_q, cur_t
    terminal: StackIndex = (prev_q, _G, prev_t, _G)
    stacks[terminal] = stacks.get(terminal, 0) + 1

    return replace(pair, energy=duplex.energy, stacks=stacks)


# ---------------------------------------------------------------------------
# Pipeline phases
# ---------------------------------------------------------------------------


def generate_corpus(
    dsm_id: DsmId, temperature: float, workers: int | None = None
) -> list[Duplex]:
    """Generate and score a corpus at `temperature` in degrees Celsius."""
    random.seed(SEED)
    logger.info("generating corpus for {} at {}°C", dsm_id.id, temperature)
    scored: list[Duplex] = []
    with Pool(
        workers, initializer=_init_worker, initargs=(dsm_id, temperature)
    ) as pool:
        for strategy in STRATEGIES:
            pairs = generate_duplex_pairs(strategy, temperature)
            scored.extend(pool.map(score_duplex, pairs))
            logger.debug("{}: {} pairs scored", strategy.name, len(pairs))
    logger.info("corpus: {} duplexes", len(scored))
    return scored


def _class_stacks(
    stacks: dict[StackIndex, int], *, strand_symmetric: bool
) -> dict[StackIndex, int]:
    table = _CLASS_TABLE[strand_symmetric]
    classes: dict[StackIndex, int] = {}
    for k, v in stacks.items():
        classes[table[k]] = classes.get(table[k], 0) + v
    return classes


def _regress(
    energies: list[float], stacks_list: list[dict[StackIndex, int]]
) -> dict[StackIndex, float]:
    """Fit energy ~ stack counts via OLS, drop weak features, return coefficients."""
    keys = sorted({k for s in stacks_list for k in s}, key=repr)
    key_idx = {k: i for i, k in enumerate(keys)}

    X = np.zeros((len(energies), len(keys)))
    y = np.array(energies)
    for i, s in enumerate(stacks_list):
        for k, v in s.items():
            X[i, key_idx[k]] = v

    col_sums = X.sum(axis=0)
    strong = [j for j in range(len(keys)) if col_sums[j] >= _MIN_STACK_COUNT]
    weak_mask = col_sums < _MIN_STACK_COUNT
    row_mask: NDArray[np.bool_] = (X[:, weak_mask] == 0).all(axis=1)
    X_fit: NDArray[np.float64] = X[row_mask][:, strong]
    y_fit: NDArray[np.float64] = y[row_mask]

    regr = LinearRegression(fit_intercept=False).fit(X_fit, y_fit)
    logger.info(
        "regress: {} rows, {} features, rank {}", len(y_fit), X_fit.shape[1], regr.rank_
    )
    return {keys[j]: round(float(regr.coef_[i]), 4) for i, j in enumerate(strong)}


def fit(scored: list[Duplex], *, strand_symmetric: bool) -> dict[StackIndex, float]:
    """Linear regression: energy ~ symmetry-class stack counts, no intercept."""
    active = [
        (d.energy, d.stacks) for d in scored if d.energy is not None and d.energy < 0
    ]
    energies = [e for e, _ in active]

    class_stacks = [
        _class_stacks(s, strand_symmetric=strand_symmetric) for _, s in active
    ]
    return _regress(energies, class_stacks)


_N_SUBS = (Base.A, Base.C, Base.G, Base.U)


@cache
def _expand_n(stack: StackIndex) -> tuple[StackIndex, ...]:
    """Expand N bases to all ACGU combinations. No-N stacks return (self,)."""
    if Base.N not in stack:
        return (stack,)
    result: list[StackIndex] = []
    for combo in product(*(_N_SUBS if b is Base.N else (b,) for b in stack)):
        result.append((combo[0], combo[1], combo[2], combo[3]))
    return tuple(result)


def build_matrix(coefs: dict[StackIndex, float], *, strand_symmetric: bool) -> DsmData:
    """Fill the 36x36 matrix from regression coefficients. N entries get max over ACGU."""
    table = _CLASS_TABLE[strand_symmetric]
    dsm = DsmData(np.full((DSM_N, DSM_N), DsmData.UNOBSERVED))
    observed = 0
    for stack in product(Base, Base, Base, Base):
        val = max(coefs.get(table[k], DsmData.UNOBSERVED) for k in _expand_n(stack))
        dsm.set(stack, val)
        if val < DsmData.UNOBSERVED:
            observed += 1
    logger.info("matrix: {}/{} cells observed", observed, DSM_N * DSM_N)
    return dsm


def postprocess(raw: DsmData, dsm_id: DsmId) -> DsmData:
    """Symmetrize equal-chemistry strands, then adjust the initiation offset."""
    mat = raw.symmetrize() if dsm_id.params.strand_symmetric else raw
    return mat.reoffset(mat.auto_offset())


def build_table(
    dsm_id: DsmId, temperature: float, workers: int | None = None
) -> DsmTable:
    """Full pipeline: generate corpus → fit → build matrix → postprocess → DsmTable."""
    strand_sym = dsm_id.params.strand_symmetric
    corpus = generate_corpus(dsm_id, temperature, workers)
    coefs = fit(corpus, strand_symmetric=strand_sym)
    raw = build_matrix(coefs, strand_symmetric=strand_sym)
    dsm = postprocess(raw, dsm_id)
    return DsmTable(dsm_id.id, dsm_id.params, temperature, dsm)


def run_pipeline(
    dsm_id: DsmId, temperature: float, output: Path, workers: int | None = None
) -> DsmTable:
    table = build_table(dsm_id, temperature, workers)
    table.save(output)
    logger.success("{}/{}/{}.tsv", output, dsm_id.id, int(temperature))
    return table


def generate_all(output: Path, workers: int | None = None) -> list[DsmTable]:
    return [
        run_pipeline(dsm_id, temperature, output, workers)
        for dsm_id in ALL_DSM_IDS
        for temperature in dsm_id.temperatures
    ]


_RUST_TABLES = Path(__file__).parents[2] / "src" / "dsm" / "tables"


def emit_rust_mod(path: Path) -> None:
    """Write `mod.rs`: one `mod` line per table file plus the `TABLES` registry."""
    mods = [(d.id, int(t)) for d in ALL_DSM_IDS for t in d.temperatures]
    lines = [
        "//! Bundled DSM tables.",
        "//! Generated by tools/dsm/pipeline.py",
        "",
        "use super::DsmTable;",
        "use crate::types::Energy;",
        "",
        *(f"mod {n};" for n in sorted(f"{i}_{t}" for i, t in mods)),
        "",
        "pub(crate) static TABLES: &[(&str, i32, Energy, &DsmTable)] = &[",
        *(f'    ("{i}", {t}, {i}_{t}::INITIATION, &{i}_{t}::TABLE),' for i, t in mods),
        "];",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

app = typer.Typer()

_DSM_BY_ID = {d.id: d for d in ALL_DSM_IDS}


@app.command("all")
def cmd_all(
    output: Annotated[Path, typer.Option(help="Output directory")] = Path("."),
    workers: Annotated[int | None, typer.Option(help="Pool size")] = None,
) -> None:
    _ = generate_all(output, workers)


@app.command("one")
def cmd_one(
    dsm_id: Annotated[str, typer.Argument(help="DSM ID")],
    temperature: Annotated[float, typer.Argument(help="Temperature in Celsius")],
    output: Annotated[Path, typer.Option(help="Output directory")] = Path("."),
    workers: Annotated[int | None, typer.Option(help="Pool size")] = None,
) -> None:
    """Generate one DSM at `temperature` in degrees Celsius."""
    did = _DSM_BY_ID.get(dsm_id)
    if did is None:
        logger.error("unknown DSM ID '{}', valid: {}", dsm_id, ", ".join(_DSM_BY_ID))
        raise SystemExit(1)
    _ = run_pipeline(did, temperature, output, workers)


@app.command("rust")
def cmd_rust(
    output: Annotated[Path, typer.Option(help="Rust output directory")] = _RUST_TABLES,
    workers: Annotated[int | None, typer.Option(help="Pool size")] = None,
) -> None:
    """Generate every `{id}_{temperature}.rs` table for the runtime, then `mod.rs`."""
    for dsm_id in ALL_DSM_IDS:
        for temperature in dsm_id.temperatures:
            path = output / f"{dsm_id.id}_{int(temperature)}.rs"
            build_table(dsm_id, temperature, workers).emit_rust(path)
            logger.success("{}", path)
    emit_rust_mod(output / "mod.rs")
    logger.success("{}", output / "mod.rs")


if __name__ == "__main__":
    app()
