//! Entry point for the `pathotypr` application.
//!
//! This binary orchestrates the different subcommands available in the tool,
//! such as `train`, `predict`, `classify`, and `SplitFastq`. It uses the `clap`
//! crate to parse command-line arguments and dispatches control to the
//! corresponding module.

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use clap::{Parser, Subcommand};
use log::LevelFilter;
use std::sync::Once;
// --- NEW IMPORTS ---
// Import chrono to get the local time
use chrono::Local;
// Import the Write trait to format the log message
use std::io::Write;

// Internal module declarations.
// Items marked `pub` in these modules are part of the library API consumed by
// Tauri — they appear unused from the binary's perspective, so we suppress
// dead_code per-module rather than with a blanket crate-level allow.
#[allow(dead_code, unused_imports)] mod classify;
#[allow(dead_code)] mod classify_split_fastq;
#[allow(dead_code)] mod common;
#[allow(dead_code)] mod errors;
#[allow(dead_code)] mod excel;
#[allow(dead_code)] mod fasta_io;
#[allow(dead_code)] mod lineage;
#[allow(dead_code, unused_imports)] mod r#match;
#[allow(dead_code)] mod model;
#[allow(dead_code)] mod paired_end;
#[allow(dead_code)] mod predict;
#[allow(dead_code)] mod sparse_tree;
#[allow(dead_code)] mod split_kmer;
#[allow(dead_code)] mod train;
#[allow(dead_code)] mod vectorizer;

// Use the AppResult type from the errors module.
use errors::AppResult;

/* ----------------- Command-Line Interface Definition ----------------- */

#[derive(Parser)]
#[command(
    name = "pathotypr",
    version = env!("CARGO_PKG_VERSION"),
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
        Commands::Train(a) => { train::run(a)?; }
        Commands::Predict(a) => predict::run(a)?,
        Commands::Classify(a) => classify::run(a)?,
        Commands::SplitFastq(a) => classify_split_fastq::run(a)?,
        Commands::Match(a) => r#match::run(a)?,
    }
    Ok(())
}
