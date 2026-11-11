//! Command-line entry point for `risearch`.

mod cli;

use clap::Parser;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> anyhow::Result<()> {
    cli::app::run(cli::Cli::parse())
}
