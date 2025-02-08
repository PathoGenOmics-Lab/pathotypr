use clap::{Parser, Subcommand};

/// Pathotypr tool: classify genomes using three subcommands: train, predict, and classify.
#[derive(Parser)]
#[command(name = "pathotypr")]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Train a model using an input FASTA file
    Train(train::Args),
    /// Predict genome classifications using a saved model and an input FASTA file
    Predict(predict::PredictArgs),
    /// Classify genomes (with markers) using TSV/FASTA inputs
    Classify(classify::Args),
}

mod train;
mod predict;
mod classify;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Train(args) => {
            train::run(args)?;
        }
        Commands::Predict(args) => {
            predict::run(args)?;
        }
        Commands::Classify(args) => {
            classify::run(args)?;
        }
    }
    Ok(())
}
