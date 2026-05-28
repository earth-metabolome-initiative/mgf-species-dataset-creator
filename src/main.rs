//! Command-line entrypoint for the MGF species dataset creator.

use std::path::PathBuf;

use clap::Parser;
use mgf_species_dataset_creator::{PipelineConfig, Polarity, run_pipeline};
use tracing_subscriber::EnvFilter;

/// Top-level command-line parser.
#[derive(Debug, Parser)]
#[command(author, version, about)]
struct Cli {
    /// Selected command.
    #[command(subcommand)]
    command: Command,
}

/// Supported CLI commands.
#[derive(Debug, clap::Subcommand)]
enum Command {
    /// Download/cache the Zenodo dataset and write taxonomically enriched MGF files.
    Run(RunArgs),
}

/// Arguments for the `run` command.
#[derive(Debug, Parser)]
struct RunArgs {
    /// Zenodo record ID to fetch.
    #[arg(long, default_value_t = mgf_species_dataset_creator::DEFAULT_RECORD_ID)]
    record_id: u64,
    /// Archive file name to fetch from the Zenodo record.
    #[arg(long, default_value = mgf_species_dataset_creator::DEFAULT_ZENODO_ARCHIVE)]
    archive_name: String,
    /// Directory where the archive is cached.
    #[arg(long, default_value = ".cache/zenodo")]
    cache_dir: PathBuf,
    /// Directory where the archive is extracted.
    #[arg(long, default_value = ".cache/extracted")]
    extraction_dir: PathBuf,
    /// Use an existing extracted dataset root instead of downloading/extracting.
    #[arg(long)]
    extracted_dataset_dir: Option<PathBuf>,
    /// Directory containing NCBI taxdump nodes.dmp, names.dmp, and optionally merged.dmp.
    #[arg(long)]
    taxdump_dir: PathBuf,
    /// Directory where enriched MGF files and reports are written.
    #[arg(long, default_value = "output")]
    output_dir: PathBuf,
    /// Which polarity to process.
    #[arg(long, default_value_t = Polarity::Both)]
    polarity: Polarity,
    /// Replace existing enriched MGF files.
    #[arg(long)]
    overwrite: bool,
}

/// Parses CLI arguments, runs the selected command, and reports completion.
#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    tracing_subscriber::fmt()
        .with_env_filter(EnvFilter::from_default_env())
        .init();
    let cli = Cli::parse();
    match cli.command {
        Command::Run(args) => {
            let output_dir = args.output_dir.clone();
            let report = run_pipeline(PipelineConfig {
                record_id: args.record_id,
                archive_name: args.archive_name,
                cache_dir: args.cache_dir,
                extraction_dir: args.extraction_dir,
                extracted_dataset_dir: args.extracted_dataset_dir,
                taxdump_dir: args.taxdump_dir,
                output_dir: args.output_dir,
                polarity: args.polarity,
                overwrite: args.overwrite,
            })
            .await?;
            println!(
                "processed {} MGF file(s); report written to {}/pipeline_report.json",
                report.processed_mgfs.len(),
                output_dir.display()
            );
        }
    }
    Ok(())
}
