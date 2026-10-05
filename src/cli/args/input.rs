use std::path::PathBuf;

#[derive(clap::Args, Debug, Clone)]
pub(crate) struct InputArgs {
    /// Query sequences: FASTA or FASTQ file, optionally compressed, or '-' for stdin
    ///
    /// gzip, bzip2, xz and zstd are detected from the file contents. Letters are case-insensitive
    /// and T and U are the same base; '-' and '.' are dropped and any other letter becomes N.
    /// Record IDs must be unique.
    #[arg(short = 'q', long = "query", value_name = "FILE")]
    pub(crate) query: PathBuf,

    /// Target index built by `risearch index`
    #[arg(
        short = 't',
        long = "target",
        value_name = "TARGET",
        required = true,
        conflicts_with = "legacy_target"
    )]
    pub(crate) target: Option<PathBuf>,
}
