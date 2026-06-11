//! Utilities for preparing and publishing a corrected PF1600 archive.

use std::fs::File;
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use serde::Serialize;
use tar::{Archive, Builder};
use tracing::info;
use zenodo_rs::{
    AccessRight, Auth, Creator, DepositMetadataUpdate, DepositionId, Endpoint, FileReplacePolicy,
    UploadSpec, UploadType, ZenodoClient,
};

use crate::Error;
use crate::error::{IoResultExt, Result};
use crate::progress::ProgressReporter;

/// Corrected MassIVE URL for `VGF138_A01_features_ms2_pos.mgf`.
pub const VGF138_A01_CORRECTED_URL: &str = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000087728%2Fupdates%2F2021-12-03_pmallard_cff635d5%2Fpeak%2Findividual_mgf_files%2FVGF138_A01_features_ms2_pos.mgf&forceDownload=true";

/// Corrected MassIVE URL for `VGF151_E05_features_ms2_pos.mgf`.
pub const VGF151_E05_CORRECTED_URL: &str = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000087728%2Fupdates%2F2021-12-03_pmallard_cff635d5%2Fpeak%2Findividual_mgf_files%2FVGF151_E05_features_ms2_pos.mgf&forceDownload=true";

/// Configuration for preparing a corrected PF1600 archive.
#[derive(Debug, Clone)]
pub struct RepairConfig {
    /// Root of the extracted `pf1600_raw` dataset.
    pub source_dataset_dir: PathBuf,
    /// Working directory where the corrected dataset copy is prepared.
    pub work_dataset_dir: PathBuf,
    /// Path of the corrected `.tar.gz` archive to create.
    pub archive_path: PathBuf,
    /// Whether to replace an existing working copy or archive.
    pub overwrite: bool,
    /// Gzip compression level from 0 to 9.
    pub gzip_level: u32,
}

/// Configuration for publishing the corrected archive to Zenodo.
#[derive(Debug, Clone)]
pub struct PublishConfig {
    /// Zenodo API token.
    pub token: String,
    /// Whether to use Zenodo sandbox rather than production Zenodo.
    pub sandbox: bool,
    /// Whether to create a separate new dataset instead of versioning an existing record.
    pub new_dataset: bool,
    /// Published Zenodo deposition/record ID to version.
    pub record_id: u64,
    /// Zenodo dataset title.
    pub title: String,
    /// Zenodo creator display name.
    pub creator: String,
    /// HTML description for the Zenodo record.
    pub description_html: String,
}

/// Result of preparing and optionally publishing a corrected archive.
#[derive(Debug, Clone, Serialize)]
pub struct RepairReport {
    /// Corrected dataset working directory.
    pub work_dataset_dir: PathBuf,
    /// Corrected archive path.
    pub archive_path: PathBuf,
    /// Paths replaced in the working copy.
    pub replaced_files: Vec<PathBuf>,
    /// Legacy Sirius artifact paths preserved in the working copy.
    pub preserved_legacy_files: Vec<PathBuf>,
    /// Published Zenodo record ID, when publishing was requested.
    pub zenodo_record_id: Option<u64>,
    /// Published Zenodo DOI, when available.
    pub zenodo_doi: Option<String>,
}

/// Prepares a corrected PF1600 archive and optionally publishes it to Zenodo.
pub async fn repair_dataset(
    config: RepairConfig,
    publish: Option<PublishConfig>,
    progress: &ProgressReporter,
) -> Result<RepairReport> {
    validate_repair_config(&config)?;
    prepare_working_copy(&config, progress)?;
    let (replaced_files, preserved_legacy_files) =
        download_replacements(&config.work_dataset_dir, progress).await?;
    create_archive(
        &config.work_dataset_dir,
        &config.archive_path,
        config.overwrite,
        config.gzip_level,
        progress,
    )?;
    let (zenodo_record_id, zenodo_doi) = if let Some(publish) = publish {
        publish_archive_version(&config.archive_path, publish, progress).await?
    } else {
        (None, None)
    };
    Ok(RepairReport {
        work_dataset_dir: config.work_dataset_dir,
        archive_path: config.archive_path,
        replaced_files,
        preserved_legacy_files,
        zenodo_record_id,
        zenodo_doi,
    })
}

/// Validates source and output paths before doing any expensive work.
fn validate_repair_config(config: &RepairConfig) -> Result<()> {
    if !config.source_dataset_dir.exists() {
        return Err(Error::MissingPath(config.source_dataset_dir.clone()));
    }
    if config.gzip_level > 9 {
        return Err(Error::MissingValue(format!(
            "gzip level must be between 0 and 9, got {}",
            config.gzip_level
        )));
    }
    Ok(())
}

/// Creates a hard-linked working copy of the extracted dataset.
fn prepare_working_copy(config: &RepairConfig, progress: &ProgressReporter) -> Result<()> {
    if config.work_dataset_dir.exists() {
        if !config.overwrite {
            return Err(Error::OutputExists(config.work_dataset_dir.clone()));
        }
        std::fs::remove_dir_all(&config.work_dataset_dir).with_path(&config.work_dataset_dir)?;
    }
    let spinner = progress.spinner(format!(
        "Counting files in {} before copy",
        config.source_dataset_dir.display()
    ));
    let file_count = count_files(&config.source_dataset_dir)?;
    spinner.finish_with_message(format!("Found {file_count} files to copy"));

    let bar = progress.bar(
        file_count,
        format!(
            "Copying dataset to {} with hard links",
            config.work_dataset_dir.display()
        ),
    );
    copy_tree_hardlinking(&config.source_dataset_dir, &config.work_dataset_dir, &bar)?;
    bar.finish_with_message(format!(
        "Prepared corrected dataset copy at {}",
        config.work_dataset_dir.display()
    ));
    Ok(())
}

/// Recursively copies a directory tree, using hard links for files when possible.
fn copy_tree_hardlinking(
    source: &Path,
    destination: &Path,
    progress: &indicatif::ProgressBar,
) -> Result<()> {
    std::fs::create_dir_all(destination).with_path(destination)?;
    for entry in std::fs::read_dir(source).with_path(source)? {
        let entry = entry.with_path(source)?;
        let source_path = entry.path();
        let destination_path = destination.join(entry.file_name());
        if source_path.is_dir() {
            copy_tree_hardlinking(&source_path, &destination_path, progress)?;
        } else if std::fs::hard_link(&source_path, &destination_path).is_err() {
            std::fs::copy(&source_path, &destination_path).with_path(&destination_path)?;
            progress.inc(1);
        } else {
            progress.inc(1);
        }
    }
    Ok(())
}

/// Counts regular files under a directory tree.
fn count_files(root: &Path) -> Result<u64> {
    let mut count = 0;
    for entry in std::fs::read_dir(root).with_path(root)? {
        let entry = entry.with_path(root)?;
        let path = entry.path();
        if path.is_dir() {
            count += count_files(&path)?;
        } else {
            count += 1;
        }
    }
    Ok(count)
}

/// Restores the two original Sirius artifact files from the Zenodo archive.
pub fn restore_legacy_sirius_files(
    original_archive_path: &Path,
    work_dataset_dir: &Path,
    progress: &ProgressReporter,
) -> Result<Vec<PathBuf>> {
    let targets = [
        (
            Path::new("VGF138_A01/pos/VGF138_A01_features_ms2_pos.mgf"),
            Path::new("VGF138_A01/pos/VGF138_A01_sirius_pos.mgf"),
        ),
        (
            Path::new("VGF151_E05/pos/VGF151_E05_features_ms2_pos.mgf"),
            Path::new("VGF151_E05/pos/VGF151_E05_sirius_pos.mgf"),
        ),
    ];
    let spinner = progress.spinner(format!(
        "Restoring legacy Sirius artifacts from {} (sequential archive scan)",
        original_archive_path.display()
    ));
    let file = File::open(original_archive_path).with_path(original_archive_path)?;
    let decoder = GzDecoder::new(file);
    let mut archive = Archive::new(decoder);
    let mut restored = Vec::new();

    for entry in archive.entries().with_path(original_archive_path)? {
        let mut entry = entry.with_path(original_archive_path)?;
        let entry_path = entry.path().with_path(original_archive_path)?.into_owned();
        let Some((_, legacy_relative_path)) = targets
            .iter()
            .find(|(source_relative_path, _)| entry_path.ends_with(source_relative_path))
        else {
            continue;
        };
        let destination = work_dataset_dir.join(legacy_relative_path);
        if let Some(parent) = destination.parent() {
            std::fs::create_dir_all(parent).with_path(parent)?;
        }
        entry.unpack(&destination).with_path(&destination)?;
        restored.push(destination);
        if restored.len() == targets.len() {
            break;
        }
    }

    if restored.len() != targets.len() {
        return Err(Error::MissingValue(format!(
            "restored {} of {} legacy Sirius artifact files from {}",
            restored.len(),
            targets.len(),
            original_archive_path.display()
        )));
    }

    spinner.finish_with_message(format!(
        "Restored {} legacy Sirius artifact file(s)",
        restored.len()
    ));
    Ok(restored)
}

/// Preserves originals as Sirius artifacts and downloads corrected MassIVE files.
async fn download_replacements(
    work_dataset_dir: &Path,
    progress: &ProgressReporter,
) -> Result<(Vec<PathBuf>, Vec<PathBuf>)> {
    let replacements = [
        (
            VGF138_A01_CORRECTED_URL,
            PathBuf::from("VGF138_A01/pos/VGF138_A01_features_ms2_pos.mgf"),
            PathBuf::from("VGF138_A01/pos/VGF138_A01_sirius_pos.mgf"),
        ),
        (
            VGF151_E05_CORRECTED_URL,
            PathBuf::from("VGF151_E05/pos/VGF151_E05_features_ms2_pos.mgf"),
            PathBuf::from("VGF151_E05/pos/VGF151_E05_sirius_pos.mgf"),
        ),
    ];
    let client = reqwest::Client::new();
    let mut replaced = Vec::new();
    let mut preserved = Vec::new();
    for (url, relative_path, legacy_relative_path) in replacements {
        let destination = work_dataset_dir.join(&relative_path);
        if !destination.exists() {
            return Err(Error::MissingPath(destination));
        }
        let legacy_destination = work_dataset_dir.join(&legacy_relative_path);
        std::fs::copy(&destination, &legacy_destination).with_path(&legacy_destination)?;
        preserved.push(legacy_destination);

        let spinner =
            progress.spinner(format!("Downloading corrected {}", relative_path.display()));
        let bytes = client
            .get(url)
            .send()
            .await?
            .error_for_status()?
            .bytes()
            .await?;
        std::fs::write(&destination, &bytes).with_path(&destination)?;
        spinner.finish_with_message(format!(
            "Downloaded corrected {} ({} bytes)",
            relative_path.display(),
            bytes.len()
        ));
        replaced.push(destination);
    }
    Ok((replaced, preserved))
}

/// Creates a `.tar.gz` archive containing the corrected dataset directory.
fn create_archive(
    work_dataset_dir: &Path,
    archive_path: &Path,
    overwrite: bool,
    gzip_level: u32,
    progress: &ProgressReporter,
) -> Result<()> {
    if archive_path.exists() {
        if !overwrite {
            return Err(Error::OutputExists(archive_path.to_path_buf()));
        }
        std::fs::remove_file(archive_path).with_path(archive_path)?;
    }
    if let Some(parent) = archive_path.parent() {
        std::fs::create_dir_all(parent).with_path(parent)?;
    }
    let spinner = progress.spinner(format!(
        "Counting files before archive creation for {}",
        work_dataset_dir.display()
    ));
    let file_count = count_files(work_dataset_dir)?;
    spinner.finish_with_message(format!("Archiving {file_count} files"));

    let bar = progress.bar(
        file_count,
        format!(
            "Creating corrected archive {} with gzip level {}",
            archive_path.display(),
            gzip_level
        ),
    );
    let file = File::create(archive_path).with_path(archive_path)?;
    let encoder = GzEncoder::new(file, Compression::new(gzip_level));
    let mut archive = Builder::new(encoder);
    let archive_root = work_dataset_dir
        .file_name()
        .ok_or_else(|| Error::MissingValue("work dataset directory basename".to_string()))?;
    append_tree_with_progress(
        &mut archive,
        work_dataset_dir,
        Path::new(archive_root),
        &bar,
    )?;
    archive.finish().with_path(archive_path)?;
    bar.finish_with_message(format!(
        "Created corrected archive {}",
        archive_path.display()
    ));
    Ok(())
}

/// Recursively appends a directory tree to a tar archive with file-level progress.
fn append_tree_with_progress<W: std::io::Write>(
    archive: &mut Builder<W>,
    source: &Path,
    archive_path: &Path,
    progress: &indicatif::ProgressBar,
) -> Result<()> {
    archive.append_dir(archive_path, source).with_path(source)?;
    for entry in std::fs::read_dir(source).with_path(source)? {
        let entry = entry.with_path(source)?;
        let source_path = entry.path();
        let child_archive_path = archive_path.join(entry.file_name());
        if source_path.is_dir() {
            append_tree_with_progress(archive, &source_path, &child_archive_path, progress)?;
        } else {
            archive
                .append_path_with_name(&source_path, &child_archive_path)
                .with_path(&source_path)?;
            progress.inc(1);
        }
    }
    Ok(())
}

/// Publishes an existing corrected archive as a new version of an existing Zenodo dataset.
pub async fn publish_archive_version(
    archive_path: &Path,
    config: PublishConfig,
    progress: &ProgressReporter,
) -> Result<(Option<u64>, Option<String>)> {
    if !archive_path.exists() {
        return Err(Error::MissingPath(archive_path.to_path_buf()));
    }
    let spinner = if config.new_dataset {
        progress.spinner("Publishing corrected archive as a new Zenodo dataset")
    } else {
        progress.spinner(format!(
            "Publishing corrected archive as a new version of Zenodo record {}",
            config.record_id
        ))
    };
    let mut builder = ZenodoClient::builder(Auth::new(config.token));
    if config.sandbox {
        builder = builder.endpoint(Endpoint::Sandbox);
    }
    let client = builder.build()?;
    let metadata = DepositMetadataUpdate::builder()
        .title(config.title)
        .upload_type(UploadType::Dataset)
        .description_html(config.description_html)
        .creator(Creator::named(config.creator))
        .access_right(AccessRight::Open)
        .build()
        .map_err(|error| Error::MissingValue(format!("invalid Zenodo metadata: {error}")))?;
    let uploaded_name = archive_path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| Error::MissingValue("archive filename".to_string()))?
        .to_string();
    let upload = UploadSpec::from_path_as(archive_path, uploaded_name).with_path(archive_path)?;
    let published = if config.new_dataset {
        client
            .create_and_publish_dataset(&metadata, vec![upload])
            .await?
    } else {
        client
            .publish_dataset_with_policy(
                DepositionId(config.record_id),
                &metadata,
                FileReplacePolicy::ReplaceAll,
                vec![upload],
            )
            .await?
    };
    let record_id = published.record.id.0;
    let doi = published.record.doi.map(|doi| doi.to_string());
    if config.new_dataset {
        spinner.finish_with_message(format!(
            "Published corrected dataset as Zenodo record {record_id}"
        ));
    } else {
        spinner.finish_with_message(format!(
            "Published corrected dataset version as Zenodo record {record_id}"
        ));
    }
    info!(
        record_id,
        ?doi,
        "published corrected Zenodo dataset version"
    );
    Ok((Some(record_id), doi))
}
