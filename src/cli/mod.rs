pub(crate) mod app;
pub(crate) mod args;

use clap::{Parser, Subcommand};
use std::num::NonZeroUsize;
use std::path::PathBuf;

pub(crate) use args::SearchArgs;

pub(crate) fn main() -> anyhow::Result<()> {
    app::run(Cli::parse())
}

#[derive(Parser, Debug)]
#[command(name = "risearch")]
#[command(
    author,
    version,
    about = "Energy based RNA-RNA interaction predictions",
    long_about = None,
    subcommand_required = true,
    arg_required_else_help = true,
    after_help = "Subcommand help:\n  risearch index --help\n  risearch search --help"
)]
pub(crate) struct Cli {
    /// Increase verbosity (-v = info, -vv = debug, -vvv = trace)
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, global = true)]
    pub(crate) verbose: u8,

    /// Number of threads (global): search workers, and OpenMP threads for
    /// `index` in openmp builds. Overrides RAYON_NUM_THREADS / OMP_NUM_THREADS.
    #[arg(
        short = 'j',
        long = "jobs",
        alias = "threads",
        value_name = "N",
        global = true
    )]
    pub(crate) jobs: Option<NonZeroUsize>,

    #[command(subcommand)]
    pub(crate) command: Commands,
}

#[derive(clap::Args, Debug)]
pub(crate) struct IndexCommand {
    /// Input file in FASTA format.
    #[arg(value_name = "INPUT")]
    pub(crate) input: PathBuf,

    /// Save index to given index file path
    #[arg(value_name = "OUTPUT")]
    pub(crate) output: PathBuf,
}

#[derive(Subcommand, Debug)]
#[allow(clippy::large_enum_variant)] // CLI parsing - allocation overhead is negligible
pub(crate) enum Commands {
    /// Create index for target sequence(s)
    Index(IndexCommand),

    /// Search for interactions in the given sequence(s)
    Search(SearchArgs),
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::Parser;

    #[test]
    fn jobs_must_be_at_least_one() {
        assert!(Cli::try_parse_from(["risearch", "-j", "0", "index", "in.fa", "out.idx"]).is_err());
        assert!(Cli::try_parse_from(["risearch", "-j", "1", "index", "in.fa", "out.idx"]).is_ok());
    }
}
