//! Application entry point - handles CLI dispatch and orchestration.

use std::io;
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, Result};
use log::{debug, info, trace};

use crate::cli::args::validate_output_parent;
use risearch::config::STDIO;
use risearch::fastx::{read_sequences, read_sequences_from};
use risearch::{output, search, Sequence, TargetRegistry};

use crate::cli::{Cli, Commands, SearchArgs};

pub(crate) fn run(cli: Cli) -> Result<()> {
    match &cli.command {
        Commands::Index(cmd) => {
            init_runtime(cli.verbose, cli.jobs)?;
            cmd_index(&cmd.input, &cmd.output, cli.jobs)
        }
        Commands::Search(cmd) => {
            init_runtime(cli.verbose, cli.jobs)?;
            cmd_search(cmd)
        }
    }
}

fn cmd_index(input: &Path, output: &Path, threads: Option<NonZeroUsize>) -> Result<()> {
    info!("Creating index: {:?} -> {:?}", input, output);
    validate_output_parent(output)?;
    let targets = read_input(input).context("Failed to read targets")?;
    TargetRegistry::build_to(targets, threads, output).context("Failed to create index")?;
    info!("Index saved to {:?}", output);
    Ok(())
}

fn cmd_search(cmd: &SearchArgs) -> Result<()> {
    let cmd = cmd.clone().resolve_legacy();
    let query_path = &cmd.input.query;
    let target_path = cmd.input.target.as_deref().expect("clap requires a target");
    let output_path = &cmd.output.path;

    let (opts, output) = cmd.clone().try_into_configs()?;

    debug!("Loading queries from {:?}", query_path);
    let queries = read_input(query_path).context("Failed to load queries")?;
    for (id, sequence) in &queries {
        opts.seed
            .resolve(sequence.len())
            .with_context(|| format!("Invalid seed spec for query '{id}'"))
            .context("Failed to load queries")?;
    }
    info!("Loaded {} queries", queries.len());

    debug!("Loading target index from {:?}", target_path);
    let targets = TargetRegistry::open(target_path).context("Failed to load index")?;
    trace!("Index loaded: {} targets", targets.len());

    debug!("Starting search...");

    let sink = output::TextSink::new(&queries, &targets, &output, output_path)?;
    search::run_search(&queries, &targets, &opts, &sink)?;
    info!("Done");
    Ok(())
}

fn read_input(path: &Path) -> Result<Vec<(String, Sequence)>> {
    let records = if path == Path::new(STDIO) {
        read_sequences_from(io::stdin(), "<stdin>")?
    } else {
        read_sequences(path)?
    };
    Ok(records)
}

fn init_runtime(verbosity: u8, threads: Option<NonZeroUsize>) -> Result<()> {
    init_logging(verbosity);
    init_thread_pool(threads)
}

fn init_thread_pool(threads: Option<NonZeroUsize>) -> Result<()> {
    let mut builder = rayon::ThreadPoolBuilder::new();
    // Only pin the count when the user passed one; otherwise let rayon read
    // RAYON_NUM_THREADS and fall back to all cores on its own.
    if let Some(n) = threads {
        builder = builder.num_threads(n.get());
    }
    builder
        .build_global()
        .context("Failed to initialize thread pool")?;
    debug!("Thread pool: {} threads", rayon::current_num_threads());
    Ok(())
}

pub(crate) fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };

    env_logger::Builder::new()
        .filter_level(level)
        .format_timestamp(None)
        .format_target(false)
        .format_module_path(false)
        .init();
}
