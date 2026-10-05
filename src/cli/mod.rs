pub(crate) mod app;
pub(crate) mod args;

use clap::builder::styling::{AnsiColor, Styles};
use clap::{Parser, Subcommand};
use std::num::NonZeroUsize;
use std::path::PathBuf;

pub(crate) use args::SearchArgs;

pub(crate) fn main() -> anyhow::Result<()> {
    app::run(Cli::parse())
}

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().bold())
    .usage(AnsiColor::Green.on_default().bold())
    .literal(AnsiColor::Cyan.on_default().bold())
    .placeholder(AnsiColor::Cyan.on_default());

const EXAMPLES: &[&str] = &[
    "  risearch index targets.fa targets.idx",
    "  risearch search -q queries.fa -t targets.idx -o hits.tsv",
];

const SEARCH_EXAMPLES: &[&str] = &[
    "  Search with the defaults, hits to a file:",
    "    risearch search -q queries.fa -t targets.idx -o hits.tsv",
    "  miRNA-style seeds at query positions 2-7, weaker duplexes too:",
    "    risearch search -q mirnas.fa -t utrs.idx --seed-start 2 --seed-end 7 -e -15",
    "  Up to two mismatches per seed, pairing strings, gzip output:",
    "    risearch search -q sirnas.fa -t transcripts.idx --mismatch-max 2 -f cigar -o hits.tsv.gz",
    "  RNA queries against a DNA target:",
    "    risearch search -q guides.fa -t genome.idx -P s95-rna-dna",
    "  One zstd-compressed file per query in results/:",
    "    risearch search -q queries.fa -t targets.idx -o results --multifile --compress zstd",
];

fn examples(lines: &[&str]) -> String {
    let header = STYLES.get_header();
    format!("{header}Examples:{header:#}\n{}", lines.join("\n"))
}

#[derive(Parser, Debug)]
#[command(name = "risearch")]
#[command(styles = STYLES)]
#[command(
    version,
    about = "Energy-based RNA and DNA interaction prediction",
    max_term_width = 100,
    after_help = examples(EXAMPLES)
)]
pub(crate) struct Cli {
    /// Increase verbosity (-v = info, -vv = debug, -vvv = trace; debug and trace need a debug build)
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, global = true, help_heading = "Global options")]
    pub(crate) verbose: u8,

    /// Number of worker threads [default: all cores]
    ///
    /// For index, the suffix array is built on one thread unless risearch was compiled with the
    /// openmp feature, in which case -j applies to it too.
    #[arg(
        short = 'j',
        long = "jobs",
        alias = "threads",
        value_name = "N",
        global = true,
        help_heading = "Global options"
    )]
    pub(crate) jobs: Option<NonZeroUsize>,

    #[command(subcommand)]
    pub(crate) command: Commands,
}

#[derive(clap::Args, Debug)]
pub(crate) struct IndexCommand {
    /// Target sequences: FASTA or FASTQ file, optionally compressed, or '-' for stdin
    ///
    /// gzip, bzip2, xz and zstd are detected from the file contents. Letters are case-insensitive
    /// and T and U are the same base; '-' and '.' are dropped and any other letter becomes N.
    /// Record IDs must be unique.
    #[arg(value_name = "SEQUENCES")]
    pub(crate) input: PathBuf,

    /// Path of the index file to write
    #[arg(value_name = "TARGET")]
    pub(crate) output: PathBuf,
}

#[derive(Subcommand, Debug)]
#[allow(clippy::large_enum_variant)] // CLI parsing - allocation overhead is negligible
pub(crate) enum Commands {
    /// Build a target index from RNA or DNA sequences
    #[command(arg_required_else_help = true)]
    Index(IndexCommand),

    /// Find RNA and DNA interactions between query sequences and a target index
    #[command(
        arg_required_else_help = true,
        after_help = examples(&SEARCH_EXAMPLES[..2]),
        after_long_help = examples(SEARCH_EXAMPLES)
    )]
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
