//! Entry point for the `pathotypr` application.
//!
//! This binary orchestrates the different subcommands available in the tool,
//! such as `train`, `predict`, `classify`, and `SplitFastq`. It uses the `clap`
//! crate to parse command-line arguments and dispatches control to the
//! corresponding module.

use clap::{Parser, Subcommand};
use std::sync::Once;

// Internal module declarations.
mod classify;
mod classify_split_fastq;
mod predict;
mod split_kmer;
mod train;

/* ----------------- Command-Line Interface Definition ----------------- */

#[derive(Parser)]
#[command(
    name = "pathotypr",
    version = "0.1.0",
    author = "Paula Ruiz Rodriguez",
    about = "A tool to classify genomes using machine learning and marker-based approaches."
)]
struct Cli {
    #[command(subcommand)]
    cmd: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Train a new Random Forest model on a set of genomes.
    Train(train::Args),
    /// Predict the lineage of genomes using a pre-trained model.
    Predict(predict::PredictArgs),
    /// Classify genomes based on a set of known marker k-mers.
    Classify(classify::Args),
    /// Perform split-k-mer typing directly from FASTQ files.
    SplitFastq(classify_split_fastq::SplitFastqArgs),
}

/* ----------------- Logger Initialization ----------------- */

// A static guard to ensure the logger is initialized only once.
static INIT: Once = Once::new();

/// Initializes the global logger using `env_logger`.
///
/// This function is wrapped in a `Once` block to prevent multiple initializations,
/// which would cause a panic. The default log level is "info" but can be
/// overridden by the `RUST_LOG` environment variable.
fn init_logger() {
    INIT.call_once(|| {
        let _ = env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
            .try_init();
    });
}

/* ----------------- Main Function ----------------- */

fn main() -> anyhow::Result<()> {
    // Ensure the logger is ready.
    init_logger();

    // Parse the command-line arguments and execute the corresponding subcommand.
    match Cli::parse().cmd {
        Commands::Train(a) => train::run(a).map_err(|e| anyhow::anyhow!("{}", e))?,
        Commands::Predict(a) => predict::run(a).map_err(|e| anyhow::anyhow!("{}", e))?,
        Commands::Classify(a) => classify::run(a).map_err(|e| anyhow::anyhow!("{}", e))?,
        Commands::SplitFastq(a) => classify_split_fastq::run(a)?,
    }
    Ok(())
}
