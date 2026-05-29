//! Command-line entrypoint for the MGF species dataset creator.

use std::path::PathBuf;

use clap::Parser;
use mgf_species_dataset_creator::repair::{PublishConfig, RepairConfig, repair_dataset};
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
    /// Prepare a corrected PF1600 archive and optionally publish it to Zenodo.
    RepairDataset(RepairDatasetArgs),
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
    /// Keep only the top K most intense peaks per spectrum; use 0 to keep all peaks.
    #[arg(long, default_value_t = 128)]
    top_k_peaks: usize,
    /// Replace existing enriched MGF files.
    #[arg(long)]
    overwrite: bool,
}

/// Arguments for the `repair-dataset` command.
#[derive(Debug, Parser)]
struct RepairDatasetArgs {
    /// Root of the extracted `pf1600_raw` dataset.
    #[arg(long, default_value = ".cache/extracted/pf1600_raw")]
    source_dataset_dir: PathBuf,
    /// Working directory where corrected files are substituted.
    #[arg(long, default_value = ".cache/corrected/pf1600_raw")]
    work_dataset_dir: PathBuf,
    /// Corrected archive path to create.
    #[arg(long, default_value = ".cache/corrected/pf1600_raw_corrected.tar.gz")]
    archive_path: PathBuf,
    /// Replace existing working copy and archive.
    #[arg(long)]
    overwrite: bool,
    /// Gzip level for archive creation, 0 fastest/no compression through 9 smallest/slowest.
    #[arg(long, default_value_t = 1)]
    gzip_level: u32,
    /// Publish the corrected archive to Zenodo after preparing it.
    #[arg(long)]
    publish: bool,
    /// Publish to Zenodo sandbox instead of production Zenodo.
    #[arg(long)]
    sandbox: bool,
    /// Environment variable containing the Zenodo token.
    #[arg(long, default_value = "ZENODO_TOKEN")]
    token_env: String,
    /// Zenodo dataset title.
    #[arg(
        long,
        default_value = "PF1600 raw dataset with corrected feature MS2 MGF files"
    )]
    title: String,
    /// Zenodo creator display name.
    #[arg(long, default_value = "The Earth Metabolome Initiative")]
    creator: String,
    /// Zenodo HTML description.
    #[arg(
        long,
        default_value = "<p>Corrected PF1600 raw dataset archive with replacement MassIVE MGF files for VGF138_A01 and VGF151_E05 positive-mode feature MS2 data.</p>"
    )]
    description_html: String,
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
                top_k_peaks: args.top_k_peaks,
                overwrite: args.overwrite,
            })
            .await?;
            println!(
                "processed {} MGF file(s); report written to {}/pipeline_report.json",
                report.processed_mgfs.len(),
                output_dir.display()
            );
        }
        Command::RepairDataset(args) => {
            let publish = if args.publish {
                Some(PublishConfig {
                    token: std::env::var(&args.token_env)?,
                    sandbox: args.sandbox,
                    title: args.title,
                    creator: args.creator,
                    description_html: args.description_html,
                })
            } else {
                None
            };
            let report = repair_dataset(
                RepairConfig {
                    source_dataset_dir: args.source_dataset_dir,
                    work_dataset_dir: args.work_dataset_dir,
                    archive_path: args.archive_path,
                    overwrite: args.overwrite,
                    gzip_level: args.gzip_level,
                },
                publish,
                &mgf_species_dataset_creator::progress::ProgressReporter::enabled(),
            )
            .await?;
            println!(
                "created corrected archive at {}",
                report.archive_path.display()
            );
            if let Some(record_id) = report.zenodo_record_id {
                println!("published Zenodo record {record_id}");
            }
        }
    }
    Ok(())
}
