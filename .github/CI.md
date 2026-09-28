# CI operator guide

`ci.yml` runs on pushes and pull requests and gates merges.

## ci.yml — build and test

Runs the workspace build and the test suite on every push and PR.
Uses `dtolnay/rust-toolchain@stable`, `Swatinem/rust-cache@v2`, and `taiki-e/install-action@nextest` (cargo-nextest as the test runner).
A concurrency group cancels superseded runs on the same ref.
Default job permissions are `contents: read`, `actions/checkout` runs with `persist-credentials: false`, and every job sets a `timeout-minutes`.

Cargo validation uses `--locked`.
The Python package builds through `uv run --locked`, and Maturin is configured to use the committed Cargo lockfile.
A manifest change that is not accompanied by its lockfile update therefore fails the build instead of silently resolving.
`Cargo.lock` and `bindings/python/uv.lock` are both committed for this reason.

The `supply-chain` job runs `cargo deny --locked --workspace check`, which covers all four of cargo-deny's lanes (advisories, licenses, bans, sources) against `deny.toml` at the repository root.
`[graph]` sets `all-features = true` and both release targets, so an advisory reachable only through `openmp` or only on one platform is still seen.
Licence checking sets `include-dev` and `include-build`, so a copyleft build-dependency cannot slip in unexamined; the allowlist is permissive-only, with explicit `[[licenses.exceptions]]` entries admitting GPL-3.0-only for `risearch` and `risearch-python` themselves.
`[sources]` denies unknown registries and git sources outright, while `[bans] multiple-versions` is `warn` rather than `deny`, because duplicate transitive versions are common and not by themselves a defect.
Prune `[advisories] ignore` entries when the dependency that justified them leaves the tree: a stale entry emits `warning[advisory-not-detected]` but exits zero, so it silently widens the allowlist instead of failing.

The `lint` job runs rustfmt and clippy (`-D warnings`), then two extra guards.
It builds the `openmp` feature (`cargo build --locked -p risearch --features openmp`): it is a documented, user-facing feature that the default-feature test jobs never exercise, so without this it could rot unnoticed.
`build`, not `check`, is deliberate; only a full build links `libgomp`, which is where the native OpenMP failures this lane targets actually surface.
This runs only on Ubuntu because Apple clang ships no OpenMP, so the feature cannot link there with the default toolchain.
It then builds the docs with `cargo doc --locked --no-deps -p risearch --all-features` under `RUSTDOCFLAGS: -D warnings`, so a refactor that breaks an intra-doc link fails here rather than silently degrading the published docs.rs page.

The `actions-lint` job lints the workflow files with `actionlint` (shellcheck-backed), pinned to its `rhysd/actionlint` Docker image by digest.

The `msrv` job pins the declared minimum supported Rust version (`rust-version = "1.88"` in the root `Cargo.toml`, inherited by both crates) and runs `cargo check --locked --workspace --all-targets --all-features` on it.
The floor is 1.88 because `libsais 0.2` uses `let`-chains (stabilised in Rust 1.88) and declares no `rust-version` of its own, so `cargo msrv find` is the source of truth.
Bump both the `rust-version` field and this job's pinned toolchain together when the floor moves.

The `minimal-versions` job proves the declared dependency floors of the published `risearch` crate are real: it resolves every direct dependency down to its declared minimum and builds against that.
`--locked` is intentionally absent because the tool rewrites the lock down to those minimums.
`--ignore-private` drops `risearch-python` (`publish = false`) from the graph: its `pyo3-log` dependency requires `log ~0.4.21`, which cannot unify with the `log = "0.4.8"` floor this job exists to verify, so without it the resolution fails outright and the published crate's floors are never tested.
Where `msrv` fixes the compiler (1.88) and uses the latest deps, this fixes the deps to their floors and uses the current compiler; together they bound the support envelope.

The `semver` job gates the public API of `risearch` against the pull request's base commit with `cargo-semver-checks`, so a breaking change cannot land without being noticed.
It runs on pull requests only, because `github.event.pull_request.base.sha` is what supplies the baseline; the checkout uses `fetch-depth: 0` since `--baseline-rev` resolves that commit out of the local `.git`.

`--release-type minor` is load-bearing and must not be dropped.
With no explicit release type the tool derives one from the version numbers, and when the baseline version carries a prerelease suffix it classifies the change as major; a major bump satisfies every lint's requirement, so every lint is filtered out and the job passes having checked nothing.
Both sides of a pull request read `3.0.0-alpha.1`, so that is the case here, and it will recur at every future prerelease.
Pinning the level to minor runs the breaking-change lints and skips the additive ones, which is what a stability gate wants.
`--default-features` is also deliberate: the default heuristic enables every feature that is not obviously unstable, which would pull in `openmp` and require an OpenMP toolchain on the runner for no gain, as no public item is gated on a non-default feature.
There is no `--locked`; the tool has no such flag and builds both sides through manifests it generates itself.

The gate covers only the documented tier.
`cargo-semver-checks` excludes `#[doc(hidden)]` items from the public API, which matches the three tiers the crate already maintains: plain `pub` is supported and documented, `#[doc(hidden)] pub` is reachable only because bench targets and the binary compile as separate crates and carries no stability guarantee, and everything else is `pub(crate)`.
Engine internals under `dp`, `dsm`, `seed` and `adapter` are therefore free to change without tripping this job.
To land an intentional breaking change, apply the `semver-breaking` label to the pull request, which skips the job and leaves the decision recorded on the PR.
Breaking changes are expected while the crate is in alpha, so the job exists to force that decision to be explicit rather than to forbid it.
Do not reach for `continue-on-error`: a job that reports green is a job nobody reads.

The `test-python` job builds the extension with `maturin develop --generate-stubs` into the locked `uv` environment and runs the Python suite from `bindings/python/`.
`--generate-stubs` writes `risearch/_native.pyi` from the compiled module; the stub is not committed, so ty and pyrefly run after that step, while ruff runs before it with `risearch` declared first-party.
`uv sync` prunes packages it does not track, so it must run before `maturin develop`; `uv run` does not prune, so the suite step leaves the installed extension in place.
The job pins CPython 3.10, the declared floor, because that is where the Python layer's syntax and typing assumptions break first — the `cp310-abi3` extension itself is identical on every supported interpreter.
The Rust `test` job excludes `risearch-python`, so this is the only lane that exercises the binding.
Its last step, `maturin build --sdist`, builds the source distribution and then a wheel from it, so a manifest or `include` list that breaks the sdist fails here.

The `packaging` job builds and tests the wheels with cibuildwheel on Linux and macOS; `config-settings` passes `--generate-stubs`, so every wheel ships the stub.
Nothing is uploaded to an index yet.

The `Docs` workflow (`docs.yml`) builds the Python API site with Zensical on every push to `main` and deploys it to GitHub Pages.
It generates the stub first, because mkdocstrings reads the API from it.

### OS matrix

The `test` job is a matrix over both runners, reported as `Test ubuntu` and `Test macOS`:

- `ubuntu-latest` (x86_64)
- `macos-latest` (arm64)

Both runners execute the complete suite with identical steps and no test filter.
`fail-fast: false` keeps one OS from cancelling the other, so architecture-specific indexing, scoring, and output failures are visible on both platforms.

Correctness tests use direct seed enumeration, exhaustive short alignment paths, independent scalar path scoring, and explicit output expectations.
The native dependencies of the Rust crate remain; the separate Ubuntu OpenMP build still verifies the optional indexing feature.

## Action pinning and Dependabot

Official actions from well-governed orgs (GitHub `actions/*`) use a major tag (e.g. `@v7`); every other third-party action (`dtolnay/rust-toolchain`, `Swatinem/rust-cache`, `taiki-e/install-action`, `astral-sh/setup-uv`) is pinned to a full commit SHA with a `# version` comment, and `rhysd/actionlint` to its image digest.
SHA/digest pins are immutable (a moved tag cannot inject code); Dependabot bumps them and keeps the comment current.
`astral-sh/setup-uv` additionally pins the installed `uv` version via its `version:` input, since an unpinned `setup-uv` installs the latest `uv` at run time.

Dependabot tracks three ecosystems: `cargo` (root workspace), `github-actions` (workflows), and `uv` (`bindings/python`, so `uv.lock` and the PEP 735 `[dependency-groups]` stay current).


## Extended verification

The test matrix also executes the complete Rust suite with `--release` on both platforms; doctests run once per platform in debug.
A separate Ubuntu job runs tests with `openmp` enabled and a bounded OpenMP worker count.
The Python lane builds the CLI before comparing its rows against Arrow binding results.

`verification.yml` compiles Kani harnesses on relevant PRs, runs generated correctness, Miri, native AddressSanitizer and bounded proofs on scheduled/manual runs, and produces coverage only on manual dispatch.
Miri selections fail if a named test does not exist.
These jobs add no mutation campaign.
