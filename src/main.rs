//! Command-line entrypoint for the MGF species dataset creator.

use std::path::PathBuf;

use clap::Parser;
use mgf_species_dataset_creator::plot::{PlotConfig, plot_input_stats};
use mgf_species_dataset_creator::repair::{
    PublishConfig, RepairConfig, publish_archive_version, repair_dataset,
    restore_legacy_sirius_files,
};
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
    /// Plot input MGF distributions split by positive and negative mode.
    PlotInputStats(PlotInputStatsArgs),
    /// Prepare a corrected PF1600 archive and optionally publish it to Zenodo.
    RepairDataset(RepairDatasetArgs),
    /// Publish an existing corrected archive as a new Zenodo record version.
    PublishArchive(PublishArchiveArgs),
    /// Restore the two legacy Sirius artifact files from the original Zenodo archive.
    RestoreLegacySirius(RestoreLegacySiriusArgs),
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
    /// Disable mascot-rs intensity normalization of output peaks.
    #[arg(long = "no-normalize-peak-intensities", action = clap::ArgAction::SetFalse, default_value_t = true)]
    normalize_peak_intensities: bool,
    /// Replace existing enriched MGF files.
    #[arg(long)]
    overwrite: bool,
}

/// Arguments for the `plot-input-stats` command.
#[derive(Debug, Parser)]
struct PlotInputStatsArgs {
    /// Root of the extracted input dataset.
    #[arg(long, default_value = ".cache/extracted/pf1600_raw")]
    dataset_dir: PathBuf,
    /// Directory where SVG plots and statistics files are written.
    #[arg(long, default_value = "plots")]
    output_dir: PathBuf,
    /// Which polarity to include in the plots.
    #[arg(long, default_value_t = Polarity::Both)]
    polarity: Polarity,
    /// Number of histogram bins.
    #[arg(long, default_value_t = 60)]
    bins: usize,
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
    /// Create a separate new dataset instead of publishing a new record version.
    #[arg(long)]
    new_dataset: bool,
    /// Existing Zenodo record/deposition ID to publish a new version of.
    #[arg(long, default_value_t = mgf_species_dataset_creator::DEFAULT_RECORD_ID)]
    record_id: u64,
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

/// Arguments for the `publish-archive` command.
#[derive(Debug, Parser)]
struct PublishArchiveArgs {
    /// Existing corrected archive path to publish.
    #[arg(long, default_value = ".cache/corrected/pf1600_raw_corrected.tar.gz")]
    archive_path: PathBuf,
    /// Publish to Zenodo sandbox instead of production Zenodo.
    #[arg(long)]
    sandbox: bool,
    /// Create a separate new dataset instead of publishing a new record version.
    #[arg(long)]
    new_dataset: bool,
    /// Existing Zenodo record/deposition ID to publish a new version of.
    #[arg(long, default_value_t = mgf_species_dataset_creator::DEFAULT_RECORD_ID)]
    record_id: u64,
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

/// Arguments for the `restore-legacy-sirius` command.
#[derive(Debug, Parser)]
struct RestoreLegacySiriusArgs {
    /// Original Zenodo archive containing the legacy Sirius artifact files.
    #[arg(long, default_value = ".cache/zenodo/pf1600_raw.tar.gz")]
    original_archive_path: PathBuf,
    /// Existing corrected working dataset where Sirius artifact files are restored.
    #[arg(long, default_value = ".cache/corrected/pf1600_raw")]
    work_dataset_dir: PathBuf,
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
                normalize_peak_intensities: args.normalize_peak_intensities,
                overwrite: args.overwrite,
            })
            .await?;
            println!(
                "processed {} MGF file(s); report written to {}/pipeline_report.json",
                report.processed_mgfs.len(),
                output_dir.display()
            );
        }
        Command::PlotInputStats(args) => {
            let output_dir = args.output_dir.clone();
            let report = plot_input_stats(
                PlotConfig {
                    dataset_dir: args.dataset_dir,
                    output_dir: args.output_dir,
                    polarity: args.polarity,
                    bins: args.bins,
                },
                &mgf_species_dataset_creator::progress::ProgressReporter::enabled(),
            )?;
            println!(
                "parsed {} MGF file(s), {} spectra; plots written to {}",
                report.mgf_files,
                report.spectra,
                output_dir.display()
            );
        }
        Command::RepairDataset(args) => {
            let publish = if args.publish {
                Some(PublishConfig {
                    token: env_or_dotenv(&args.token_env)?,
                    sandbox: args.sandbox,
                    new_dataset: args.new_dataset,
                    record_id: args.record_id,
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
        Command::PublishArchive(args) => {
            let (record_id, doi) = publish_archive_version(
                &args.archive_path,
                PublishConfig {
                    token: env_or_dotenv(&args.token_env)?,
                    sandbox: args.sandbox,
                    new_dataset: args.new_dataset,
                    record_id: args.record_id,
                    title: args.title,
                    creator: args.creator,
                    description_html: args.description_html,
                },
                &mgf_species_dataset_creator::progress::ProgressReporter::enabled(),
            )
            .await?;
            if let Some(record_id) = record_id {
                println!("published Zenodo record {record_id}");
            }
            if let Some(doi) = doi {
                println!("doi {doi}");
            }
        }
        Command::RestoreLegacySirius(args) => {
            let restored = restore_legacy_sirius_files(
                &args.original_archive_path,
                &args.work_dataset_dir,
                &mgf_species_dataset_creator::progress::ProgressReporter::enabled(),
            )?;
            println!("restored {} legacy Sirius file(s)", restored.len());
            for path in restored {
                println!("{}", path.display());
            }
        }
    }
    Ok(())
}

/// Loads a configuration value from the process environment or `.env`.
fn env_or_dotenv(key: &str) -> Result<String, Box<dyn std::error::Error>> {
    if let Ok(value) = std::env::var(key) {
        return Ok(value);
    }
    let dotenv = std::fs::read_to_string(".env")?;
    for line in dotenv.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((line_key, value)) = line.split_once('=') else {
            continue;
        };
        if line_key.trim() == key {
            let value = unquote_env_value(value.trim()).to_string();
            if value.is_empty() || value == "replace_with_your_zenodo_token" {
                break;
            }
            return Ok(value);
        }
    }
    Err(format!("missing {key}; set it in the environment or in .env").into())
}

/// Removes simple shell-style quotes around an environment value.
fn unquote_env_value(value: &str) -> &str {
    value
        .strip_prefix('"')
        .and_then(|value| value.strip_suffix('"'))
        .or_else(|| {
            value
                .strip_prefix('\'')
                .and_then(|value| value.strip_suffix('\''))
        })
        .unwrap_or(value)
}
