//! Entry point for the `pathotypr` application.
//!
//! This binary orchestrates the different subcommands available in the tool,
//! such as `train`, `predict`, `classify`, and `SplitFastq`. It uses the `clap`
//! crate to parse command-line arguments and dispatches control to the
//! corresponding module.

use clap::{Parser, Subcommand};
use log::LevelFilter;
use std::sync::Once;
// --- NEW IMPORTS ---
// Import chrono to get the local time
use chrono::Local;
// Import the Write trait to format the log message
use std::io::Write;


// Internal module declarations.
mod classify;
mod classify_split_fastq;
mod common;
mod errors;
// MODIFIED: Add the new 'match' module. Use r# to avoid keyword conflict.
mod r#match;
mod predict;
mod split_kmer;
mod train;

// Use the AppResult type from the errors module.
use errors::AppResult;

/* ----------------- Command-Line Interface Definition ----------------- */

#[derive(Parser)]
#[command(
    name = "pathotypr",
    version = "0.2.0",
    author = "Paula Ruiz Rodriguez",
    about = "A versatile toolkit for genome classification and variant genotyping."
)]
struct Cli {
    #[command(subcommand)]
    cmd: Commands,

    /// Set the verbosity level. Use -v for debug, -vv for trace.
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    verbose: u8,
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
    // MODIFIED: Add the 'match' subcommand.
    /// Find the best matching reference for a set of FASTQ reads.
    Match(r#match::Args),
}

/* ----------------- Logger Initialization ----------------- */

static INIT: Once = Once::new();

/// Initializes the global logger based on the verbosity flag.
fn init_logger(verbosity: u8) {
    INIT.call_once(|| {
        let level = match verbosity {
            0 => LevelFilter::Info,
            1 => LevelFilter::Debug,
            _ => LevelFilter::Trace,
        };
        
        // --- MODIFIED LOGGER INITIALIZATION ---
        // We now use a custom format to include the local timestamp.
        env_logger::Builder::new()
            .filter_level(level)
            .format(|buf, record| {
                writeln!(
                    buf,
                    "[{}] [{}] - {}",
                    // Get the current local time and format it
                    Local::now().format("%Y-%m-%d %H:%M:%S"),
                    record.level(),
                    record.args()
                )
            })
            .init();
    });
}

/* ----------------- Main Function ----------------- */

fn main() -> AppResult<()> {
    let cli = Cli::parse();
    init_logger(cli.verbose);

    // Execute the corresponding subcommand.
    match cli.cmd {
        Commands::Train(a) => train::run(a)?,
        Commands::Predict(a) => predict::run(a)?,
        Commands::Classify(a) => classify::run(a)?,
        Commands::SplitFastq(a) => classify_split_fastq::run(a)?,
        // MODIFIED: Add the handler for the 'match' subcommand.
        Commands::Match(a) => r#match::run(a)?,
    }
    Ok(())
}
