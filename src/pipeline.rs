//! End-to-end orchestration for downloading, extracting, resolving, and enriching.

use std::path::{Path, PathBuf};

use serde::Serialize;
use tracing::info;

use crate::Error;
use crate::archive::{ensure_extracted, find_mgf_files, find_taxo_outputs};
use crate::config::PipelineConfig;
use crate::error::{IoResultExt, Result};
use crate::mgf::{EnrichmentStats, enrich_mgf_file};
use crate::progress::ProgressReporter;
use crate::taxo::TaxoDataset;
use crate::taxonomy::NcbiTaxonomyResolver;
use crate::zenodo::ensure_archive;

/// Report for one processed MGF file.
#[derive(Debug, Clone, Serialize)]
pub struct ProcessedMgf {
    /// Source MGF path.
    pub input_path: PathBuf,
    /// Enriched MGF path.
    pub output_path: PathBuf,
    /// Number of spectra parsed from the source file.
    pub spectra_total: usize,
    /// Number of spectra enriched with taxonomic metadata.
    pub spectra_enriched: usize,
    /// Number of spectra without matching taxonomic metadata.
    pub spectra_without_taxo: usize,
    /// Number of enriched spectra without NCBI resolution.
    pub spectra_without_ncbi_resolution: usize,
    /// Number of most intense peaks kept per spectrum; zero means all peaks were kept.
    pub top_k_peaks: usize,
}

/// Report for a completed pipeline run.
#[derive(Debug, Clone, Serialize)]
pub struct PipelineReport {
    /// Zenodo record ID used for the run.
    pub record_id: u64,
    /// Dataset root used for discovery.
    pub dataset_root: PathBuf,
    /// Number of indexed taxonomic records.
    pub taxo_records: usize,
    /// Per-MGF processing summaries.
    pub processed_mgfs: Vec<ProcessedMgf>,
}

/// Runs the full dataset enrichment pipeline with terminal progress enabled.
pub async fn run_pipeline(config: PipelineConfig) -> Result<PipelineReport> {
    run_pipeline_with_progress(config, &ProgressReporter::enabled()).await
}

/// Runs the full dataset enrichment pipeline with the provided progress reporter.
pub async fn run_pipeline_with_progress(
    config: PipelineConfig,
    progress: &ProgressReporter,
) -> Result<PipelineReport> {
    let dataset_root = resolve_dataset_root(&config, progress).await?;

    let spinner = progress.spinner(format!(
        "Loading NCBI taxonomy from {}",
        config.taxdump_dir.display()
    ));
    let resolver = NcbiTaxonomyResolver::from_config(&config)?;
    spinner.finish_with_message("Loaded NCBI taxonomy");

    let spinner = progress.spinner("Locating sample taxo_output directories");
    let taxo_outputs = find_taxo_outputs(&dataset_root)?;
    if taxo_outputs.is_empty() {
        if is_blank_or_qc_sample(&dataset_root)? {
            std::fs::create_dir_all(&config.output_dir).with_path(&config.output_dir)?;
            let report = PipelineReport {
                record_id: config.record_id,
                dataset_root,
                taxo_records: 0,
                processed_mgfs: Vec::new(),
            };
            write_report(&config.output_dir, &report)?;
            return Ok(report);
        }
        return Err(Error::MissingPath(dataset_root.join("*/taxo_output")));
    }
    spinner.finish_with_message(format!(
        "Found {} taxo_output directorie(s)",
        taxo_outputs.len()
    ));

    let mut sample_jobs = Vec::new();
    let spinner = progress.spinner(format!(
        "Loading sample metadata and discovering {} MGF files",
        config.polarity
    ));
    let mut taxo_records = 0;
    for taxo_output in taxo_outputs {
        info!(path = %taxo_output.display(), "loading taxonomic metadata");
        let sample_root = taxo_output
            .parent()
            .ok_or_else(|| Error::MissingPath(taxo_output.clone()))?
            .to_path_buf();
        if is_blank_or_qc_sample(&sample_root)? {
            continue;
        }
        let taxo = TaxoDataset::from_taxo_output(&taxo_output, &resolver)?;
        taxo_records += taxo.len();
        let mgf_paths = find_mgf_files(&sample_root, config.polarity)?;
        if !mgf_paths.is_empty() {
            sample_jobs.push((sample_root, taxo, mgf_paths));
        }
    }
    let total_mgfs: usize = sample_jobs
        .iter()
        .map(|(_, _, mgf_paths)| mgf_paths.len())
        .sum();
    if total_mgfs == 0 {
        return Err(Error::NoMgfFiles(
            config.polarity.to_string(),
            dataset_root.clone(),
        ));
    }
    spinner.finish_with_message(format!(
        "Found {total_mgfs} MGF file(s) across {} sample directorie(s)",
        sample_jobs.len()
    ));

    std::fs::create_dir_all(&config.output_dir).with_path(&config.output_dir)?;
    let bar = progress.bar(
        total_mgfs as u64,
        if config.top_k_peaks == 0 {
            "Enriching MGF files without peak filtering".to_string()
        } else {
            format!(
                "Enriching MGF files and keeping top {} peaks",
                config.top_k_peaks
            )
        },
    );
    let mut processed_mgfs = Vec::with_capacity(total_mgfs);
    for (sample_root, taxo, mgf_paths) in sample_jobs {
        let sample_output_dir = config.output_dir.join(
            sample_root
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("sample"),
        );
        for input_path in mgf_paths {
            bar.set_message(format!(
                "Enriching {}",
                ProgressReporter::file_label(&input_path)
            ));
            info!(path = %input_path.display(), "enriching MGF");
            let stats = enrich_mgf_file(
                &input_path,
                &sample_output_dir,
                &taxo,
                config.top_k_peaks,
                config.overwrite,
            )?;
            processed_mgfs.push(processed_from_stats(stats));
            bar.inc(1);
        }
    }
    bar.finish_with_message("MGF enrichment complete");

    let report = PipelineReport {
        record_id: config.record_id,
        dataset_root,
        taxo_records,
        processed_mgfs,
    };
    let spinner = progress.spinner("Writing pipeline report");
    write_report(&config.output_dir, &report)?;
    spinner.finish_with_message(format!(
        "Report written to {}",
        config.output_dir.join("pipeline_report.json").display()
    ));
    Ok(report)
}

/// Resolves the dataset root by using existing extracted data or downloading and extracting.
async fn resolve_dataset_root(
    config: &PipelineConfig,
    progress: &ProgressReporter,
) -> Result<PathBuf> {
    if let Some(path) = &config.extracted_dataset_dir {
        if !path.exists() {
            return Err(Error::MissingPath(path.clone()));
        }
        let spinner = progress.spinner(format!("Using extracted dataset at {}", path.display()));
        spinner.finish_with_message(format!("Using extracted dataset at {}", path.display()));
        return Ok(path.clone());
    }
    let archive = ensure_archive(
        config.record_id,
        &config.archive_name,
        &config.cache_dir,
        progress,
    )
    .await?;
    let spinner = progress.spinner(format!(
        "Extracting {} into {} (this can take several minutes)",
        archive.display(),
        config.extraction_dir.display()
    ));
    ensure_extracted(&archive, &config.extraction_dir)?;
    spinner.finish_with_message(format!(
        "Dataset ready at {}",
        config.extraction_dir.display()
    ));
    Ok(config.extraction_dir.clone())
}

/// Converts internal enrichment statistics into the public report type.
fn processed_from_stats(stats: EnrichmentStats) -> ProcessedMgf {
    ProcessedMgf {
        input_path: stats.input_path,
        output_path: stats.output_path,
        spectra_total: stats.spectra_total,
        spectra_enriched: stats.spectra_enriched,
        spectra_without_taxo: stats.spectra_without_taxo,
        spectra_without_ncbi_resolution: stats.spectra_without_ncbi_resolution,
        top_k_peaks: stats.top_k_peaks,
    }
}

/// Writes the JSON pipeline report to the output directory.
fn write_report(output_dir: &Path, report: &PipelineReport) -> Result<()> {
    let path = output_dir.join("pipeline_report.json");
    let content = serde_json::to_vec_pretty(report).map_err(|source| Error::Json {
        path: path.clone(),
        source,
    })?;
    std::fs::write(&path, content).with_path(&path)
}

/// Returns whether a sample metadata file marks the sample as blank or QC.
fn is_blank_or_qc_sample(sample_root: &Path) -> Result<bool> {
    let Some(sample_id) = sample_root.file_name().and_then(|name| name.to_str()) else {
        return Ok(false);
    };
    let metadata_path = sample_root.join(format!("{sample_id}_metadata.tsv"));
    if !metadata_path.exists() {
        return Ok(false);
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(&metadata_path)
        .map_err(|source| Error::Csv {
            path: metadata_path.clone(),
            source,
        })?;
    let headers = reader
        .headers()
        .map_err(|source| Error::Csv {
            path: metadata_path.clone(),
            source,
        })?
        .clone();
    let Some(sample_type_index) = headers.iter().position(|header| header == "sample_type") else {
        return Ok(false);
    };
    for row in reader.records() {
        let row = row.map_err(|source| Error::Csv {
            path: metadata_path.clone(),
            source,
        })?;
        let sample_type = row
            .get(sample_type_index)
            .unwrap_or_default()
            .trim()
            .to_ascii_lowercase();
        if matches!(sample_type.as_str(), "blank" | "qc") {
            return Ok(true);
        }
    }
    Ok(false)
}
