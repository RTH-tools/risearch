//! The CLI and the library must report the same hits for the same run.
//!
//! Output order is deliberately unspecified, so both sides are reduced to a
//! sorted named hit key before comparison.
#![cfg(feature = "cli")]

use assert_cmd::cargo::cargo_bin_cmd;
use risearch::fastx::read_sequences;
use risearch::{run_search, Energy, SearchConfig, SearchHit, Sequence, TargetRegistry, VecSink};
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

// Mirrors the CLI arguments below; the two must stay in step or the test
// compares two different searches.
fn config() -> SearchConfig {
    let mut config = SearchConfig::default();
    config.seed.seed_length = Some(6);
    // `--mismatch-max` alone expands to prefix = suffix = max on the CLI side.
    config.seed.max_mismatches = 2;
    config.seed.min_prefix_matches = 2;
    config.seed.min_suffix_matches = 2;
    config.extend.max_extension = 30;
    config.extend.build_alignment = false;
    config.filter.delta_g = Energy::from_kcal(-5.0);
    config
}

fn collect(
    query: &[(String, Sequence)],
    target: &TargetRegistry,
    config: &SearchConfig,
) -> Vec<SearchHit> {
    let sink = VecSink::default();
    run_search(query, target, config, &sink).unwrap();
    sink.into_hits()
}

fn minimal_cli_keys(text: &str) -> Vec<String> {
    let mut keys: Vec<_> = text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            assert_eq!(fields.len(), 8, "unexpected minimal row: {line}");
            let energy = (fields[7].parse::<f64>().unwrap() * 100.0).round_ties_even() as i64;
            format!(
                "{}\0{}\0{}:{}\0{}:{}\0{}\0{}",
                fields[0],
                fields[3],
                fields[1].parse::<usize>().unwrap() - 1,
                fields[2].parse::<usize>().unwrap() - 1,
                fields[4].parse::<usize>().unwrap() - 1,
                fields[5].parse::<usize>().unwrap() - 1,
                fields[6],
                energy
            )
        })
        .collect();
    keys.sort();
    keys
}

fn build_cli_index(target: &Path, index: &Path) {
    cargo_bin_cmd!("risearch")
        .arg("index")
        .arg(target)
        .arg(index)
        .assert()
        .success();
}

fn cli_search_keys(query: &OsStr, index: &Path, stdin: Option<Vec<u8>>) -> Vec<String> {
    let mut cmd = cargo_bin_cmd!("risearch");
    cmd.arg("search")
        .arg("-q")
        .arg(query)
        .arg("-t")
        .arg(index)
        .arg("-o")
        .arg("-")
        .arg("--seed-length")
        .arg("6")
        .arg("--mismatch-max")
        .arg("2")
        .arg("-l")
        .arg("30")
        .arg("-e")
        .arg("-5.0")
        .arg("--format")
        .arg("minimal");
    if let Some(input) = stdin {
        cmd.write_stdin(input);
    }
    let output = cmd.assert().success().get_output().stdout.clone();
    minimal_cli_keys(std::str::from_utf8(&output).unwrap())
}

#[test]
fn cli_minimal_output_matches_library_hits() {
    let query = fixture("query.fa");
    let dir = TempDir::new().unwrap();
    let index = dir.path().join("targets.idx");
    build_cli_index(&fixture("target.fa"), &index);

    let cli = cli_search_keys(query.as_os_str(), &index, None);
    assert!(!cli.is_empty(), "CLI fixture must produce hits");

    let config = config();
    let queries = read_sequences(&query).unwrap();
    let targets = TargetRegistry::open(&index).unwrap();
    let mut library: Vec<_> = collect(&queries, &targets, &config)
        .into_iter()
        .map(|hit| {
            // Same tie rule as the output formatter, or exact halves split.
            let energy = (hit.energy.to_kcal() * 100.0).round_ties_even() as i64;
            format!(
                "{}\0{}\0{}:{}\0{}:{}\0{}\0{}",
                queries[hit.query_idx as usize].0,
                targets.get_name(hit.target_index()),
                hit.q_start,
                hit.q_end,
                hit.t_start,
                hit.t_end,
                hit.strand,
                energy
            )
        })
        .collect();
    library.sort();
    assert!(!library.is_empty(), "library fixture must produce hits");
    assert_eq!(cli, library, "CLI and library disagree on minimal fixture");
}

#[test]
fn cli_stdin_input_matches_file() {
    let query = fixture("query.fa");
    let dir = TempDir::new().unwrap();
    let index = dir.path().join("targets.idx");
    build_cli_index(&fixture("target.fa"), &index);

    let stdin_index = dir.path().join("stdin.idx");
    cargo_bin_cmd!("risearch")
        .arg("index")
        .arg("-")
        .arg(&stdin_index)
        .write_stdin(fs_err::read(fixture("target.fa")).unwrap())
        .assert()
        .success();
    assert_eq!(
        fs_err::read(&index).unwrap(),
        fs_err::read(&stdin_index).unwrap(),
        "stdin targets changed the index"
    );

    let from_file = cli_search_keys(query.as_os_str(), &index, None);
    let from_stdin = cli_search_keys(OsStr::new("-"), &index, Some(fs_err::read(&query).unwrap()));
    assert!(!from_file.is_empty(), "CLI fixture must produce hits");
    assert_eq!(from_file, from_stdin, "stdin query changed the hits");
}
