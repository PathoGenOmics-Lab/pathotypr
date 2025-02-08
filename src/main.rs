use clap::{Parser, Subcommand};

/// Pathotypr: A tool to classify genomes.
///
/// This tool supports three subcommands:
/// - **train**: Train a machine learning model using an input FASTA file.
/// - **predict**: Predict genome classifications using a saved model.
/// - **classify**: Classify genomes using marker data from TSV/FASTA inputs.
///
/// Author: Paula Ruiz Rodriguez <paula.ruiz.rodriguez@csic.es>
/// Version: 0.1.0
#[derive(Parser)]
#[command(
    name = "pathotypr",
    author = "Paula Ruiz Rodriguez <paula.ruiz.rodriguez@csic.es>",
    version = "0.1.0",
    about = "A tool to classify genomes using the train, predict, and classify subcommands.",
    long_about = "Pathotypr is a command-line tool developed to classify genomes. \
                  It provides various functionalities:\n\n\
                  * train   - Train a model using an input FASTA file.\n\
                  * predict - Predict genome classifications using a saved model.\n\
                  * classify - Classify genomes based on SNPs from TSV/FASTA files.\n\n\
                  Use the respective subcommand to access the desired functionality."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Train a model using an input FASTA file.
    Train(train::Args),
    /// Predict genome classifications using a saved model.
    Predict(predict::PredictArgs),
    /// Classify genomes using marker data.
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
