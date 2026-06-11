//! End-to-end orchestration for downloading, extracting, resolving, and enriching.

use std::path::{Path, PathBuf};

use mascot_rs::prelude::MGFVec;
use serde::Serialize;
use tracing::info;

use crate::Error;
use crate::archive::{ensure_extracted, find_mgf_files, find_taxo_outputs};
use crate::config::PipelineConfig;
use crate::error::{IoResultExt, Result};
use crate::mgf::{EnrichmentStats, enrich_mgf_records, write_mgf_file};
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
    /// Whether peak intensities were normalized.
    pub normalize_peak_intensities: bool,
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
    let aggregate_outputs = aggregate_outputs_for(&config);
    for output_path in aggregate_outputs.iter().flatten() {
        if output_path.exists() && !config.overwrite {
            return Err(Error::OutputExists(output_path.clone()));
        }
    }
    let bar = progress.bar(
        total_mgfs as u64,
        if config.top_k_peaks == 0 {
            "Enriching MGF files without peak filtering".to_string()
        } else {
            format!(
                "Enriching MGF files, keeping top {} peaks, normalizing intensities: {}",
                config.top_k_peaks, config.normalize_peak_intensities
            )
        },
    );
    let mut processed_mgfs = Vec::with_capacity(total_mgfs);
    let mut pos_spectra = MGFVec::default();
    let mut neg_spectra = MGFVec::default();
    for (sample_root, taxo, mgf_paths) in sample_jobs {
        for input_path in mgf_paths {
            bar.set_message(format!(
                "Enriching {}",
                ProgressReporter::file_label(&input_path)
            ));
            info!(path = %input_path.display(), "enriching MGF");
            let output_path = aggregate_output_for_path(&input_path, &aggregate_outputs)?;
            let sample_metadata = sample_metadata_headers_for_mgf(
                &sample_root,
                polarity_label_for_path(&input_path),
            )?;
            let enriched = enrich_mgf_records(
                &input_path,
                &output_path,
                &taxo,
                &sample_metadata,
                config.top_k_peaks,
                config.normalize_peak_intensities,
            )?;
            let mut spectra = enriched.spectra;
            match polarity_label_for_path(&input_path) {
                PolarityLabel::Pos => pos_spectra.append(&mut spectra),
                PolarityLabel::Neg => neg_spectra.append(&mut spectra),
            }
            processed_mgfs.push(processed_from_stats(enriched.stats));
            bar.inc(1);
        }
    }
    bar.finish_with_message("MGF enrichment complete");

    let spinner = progress.spinner("Writing aggregate polarity MGF files");
    if let Some(output_path) = &aggregate_outputs.pos {
        write_mgf_file(&pos_spectra, output_path, config.overwrite)?;
    }
    if let Some(output_path) = &aggregate_outputs.neg {
        write_mgf_file(&neg_spectra, output_path, config.overwrite)?;
    }
    spinner.finish_with_message("Aggregate polarity MGF files written");

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

/// Returns requested sample-level metadata headers for one MGF polarity.
fn sample_metadata_headers_for_mgf(
    sample_root: &Path,
    polarity: PolarityLabel,
) -> Result<Vec<(String, String)>> {
    let Some(sample_id) = sample_root.file_name().and_then(|name| name.to_str()) else {
        return Ok(Vec::new());
    };
    let metadata_path = sample_root.join(format!("{sample_id}_metadata.tsv"));
    if !metadata_path.exists() {
        return Ok(Vec::new());
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
    let Some(row) = reader
        .records()
        .next()
        .transpose()
        .map_err(|source| Error::Csv {
            path: metadata_path.clone(),
            source,
        })?
    else {
        return Ok(Vec::new());
    };

    let mut output = Vec::new();
    push_metadata_columns(
        &headers,
        &row,
        &mut output,
        &["sample_id", "sample_type", "sample_plate_id"],
    );
    match polarity {
        PolarityLabel::Pos => push_metadata_columns(
            &headers,
            &row,
            &mut output,
            &["sample_filename_pos", "pos_injection_date"],
        ),
        PolarityLabel::Neg => push_metadata_columns(
            &headers,
            &row,
            &mut output,
            &["sample_filename_neg", "neg_injection_date"],
        ),
    }
    push_metadata_columns(
        &headers,
        &row,
        &mut output,
        &[
            "bio_leish_donovani_10ugml_inhibition",
            "bio_leish_donovani_2ugml_inhibition",
            "bio_tryp_brucei_rhodesiense_10ugml_inhibition",
            "bio_tryp_brucei_rhodesiense_2ugml_inhibition",
            "bio_tryp_cruzi_10ugml_inhibition",
            "bio_l6_cytotoxicity_10ugml_inhibition",
        ],
    );
    if !output.iter().any(|(key, _)| key == "sample_id") {
        output.insert(0, ("sample_id".to_string(), sample_id.to_string()));
    }
    Ok(output)
}

/// Pushes non-empty metadata values for selected TSV columns.
fn push_metadata_columns(
    headers: &csv::StringRecord,
    row: &csv::StringRecord,
    output: &mut Vec<(String, String)>,
    columns: &[&str],
) {
    for column in columns {
        if let Some(value) = headers
            .iter()
            .position(|header| header == *column)
            .and_then(|index| row.get(index))
            .and_then(non_empty_metadata_value)
        {
            output.push(((*column).to_string(), value));
        }
    }
}

/// Returns non-empty metadata strings, excluding common missing-value markers.
fn non_empty_metadata_value(value: &str) -> Option<String> {
    let trimmed = value.trim();
    (!trimmed.is_empty() && !matches!(trimmed, "NA" | "N/A" | "NaN" | "nan"))
        .then(|| trimmed.to_string())
}

/// Aggregate output paths selected for a pipeline run.
#[derive(Debug, Clone)]
struct AggregateOutputs {
    /// Positive-mode output path, when requested.
    pos: Option<PathBuf>,
    /// Negative-mode output path, when requested.
    neg: Option<PathBuf>,
}

impl AggregateOutputs {
    /// Returns an iterator over selected output paths.
    fn iter(&self) -> impl Iterator<Item = &Option<PathBuf>> {
        [&self.pos, &self.neg].into_iter()
    }
}

/// Returns the aggregate output files requested by the polarity selector.
fn aggregate_outputs_for(config: &PipelineConfig) -> AggregateOutputs {
    AggregateOutputs {
        pos: matches!(
            config.polarity,
            crate::Polarity::Pos | crate::Polarity::Both
        )
        .then(|| config.output_dir.join("emi_taxo_enriched_pos.mgf")),
        neg: matches!(
            config.polarity,
            crate::Polarity::Neg | crate::Polarity::Both
        )
        .then(|| config.output_dir.join("emi_taxo_enriched_neg.mgf")),
    }
}

/// Returns the aggregate output file matching one source MGF path.
fn aggregate_output_for_path(path: &Path, outputs: &AggregateOutputs) -> Result<PathBuf> {
    match polarity_label_for_path(path) {
        PolarityLabel::Pos => outputs
            .pos
            .clone()
            .ok_or_else(|| Error::NoMgfFiles("pos".to_string(), path.to_path_buf())),
        PolarityLabel::Neg => outputs
            .neg
            .clone()
            .ok_or_else(|| Error::NoMgfFiles("neg".to_string(), path.to_path_buf())),
    }
}

/// Polarity labels inferred from MGF file names.
#[derive(Debug, Clone, Copy)]
enum PolarityLabel {
    /// Positive ionization mode.
    Pos,
    /// Negative ionization mode.
    Neg,
}

/// Infers the binary polarity from a feature MS2 MGF path.
fn polarity_label_for_path(path: &Path) -> PolarityLabel {
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default();
    if name.contains("_ms2_neg") || name.ends_with("_neg.mgf") {
        PolarityLabel::Neg
    } else {
        PolarityLabel::Pos
    }
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
        normalize_peak_intensities: stats.normalize_peak_intensities,
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

#[cfg(test)]
/// Tests for pipeline metadata helpers.
mod tests {
    use super::*;

    #[test]
    /// Verifies polarity-specific sample filename and injection date fields are filtered.
    fn sample_metadata_keeps_only_matching_polarity_fields() {
        let temp_dir = tempfile::tempdir().expect("temp dir");
        let sample_root = temp_dir.path().join("VGF159_A02");
        std::fs::create_dir_all(&sample_root).expect("sample dir");
        std::fs::write(
            sample_root.join("VGF159_A02_metadata.tsv"),
            concat!(
                "sample_id\tsample_type\tsample_plate_id\tsample_filename_pos\tpos_injection_date\t",
                "bio_leish_donovani_10ugml_inhibition\tsample_filename_neg\tneg_injection_date\n",
                "VGF159_A02\tsample\tplate-1\tpos.raw\t2021-01-01\t42\tneg.raw\t2021-01-02\n",
            ),
        )
        .expect("metadata");

        let pos_metadata =
            sample_metadata_headers_for_mgf(&sample_root, PolarityLabel::Pos).expect("pos");
        let neg_metadata =
            sample_metadata_headers_for_mgf(&sample_root, PolarityLabel::Neg).expect("neg");

        assert!(pos_metadata.contains(&("sample_filename_pos".to_string(), "pos.raw".to_string())));
        assert!(
            pos_metadata.contains(&("pos_injection_date".to_string(), "2021-01-01".to_string()))
        );
        assert!(
            !pos_metadata
                .iter()
                .any(|(key, _)| key == "sample_filename_neg")
        );
        assert!(
            !pos_metadata
                .iter()
                .any(|(key, _)| key == "neg_injection_date")
        );

        assert!(neg_metadata.contains(&("sample_filename_neg".to_string(), "neg.raw".to_string())));
        assert!(
            neg_metadata.contains(&("neg_injection_date".to_string(), "2021-01-02".to_string()))
        );
        assert!(
            !neg_metadata
                .iter()
                .any(|(key, _)| key == "sample_filename_pos")
        );
        assert!(
            !neg_metadata
                .iter()
                .any(|(key, _)| key == "pos_injection_date")
        );
    }
}
