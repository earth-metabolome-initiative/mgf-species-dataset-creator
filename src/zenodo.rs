//! Zenodo archive download and cache management.

use std::path::{Path, PathBuf};

use tracing::info;
use zenodo_rs::{Auth, RecordId, ZenodoClient};

use crate::error::{IoResultExt, Result};
use crate::progress::ProgressReporter;
use crate::{DEFAULT_ZENODO_ARCHIVE, Error};

/// Ensures that the configured archive exists in the local cache.
pub async fn ensure_archive(
    record_id: u64,
    archive_name: &str,
    cache_dir: &Path,
    progress: &ProgressReporter,
) -> Result<PathBuf> {
    std::fs::create_dir_all(cache_dir).with_path(cache_dir)?;
    let destination = cache_dir.join(archive_name);
    if destination.exists() {
        info!(path = %destination.display(), "using cached Zenodo archive");
        let spinner = progress.spinner(format!(
            "Using cached Zenodo archive at {}",
            destination.display()
        ));
        spinner.finish_with_message(format!(
            "Using cached Zenodo archive at {}",
            destination.display()
        ));
        return Ok(destination);
    }

    let spinner = progress.spinner(format!("Reading Zenodo record {record_id}"));
    let client = ZenodoClient::new(Auth::new(""))?;
    let record = client.get_record(RecordId(record_id)).await?;
    spinner.finish_with_message(format!("Read Zenodo record {record_id}"));

    let file = record
        .files
        .iter()
        .find(|file| file.key == archive_name)
        .or_else(|| {
            record
                .files
                .iter()
                .find(|file| file.key == DEFAULT_ZENODO_ARCHIVE)
        })
        .ok_or_else(|| Error::MissingZenodoFile {
            record_id,
            file: archive_name.to_string(),
        })?;
    let key = file.key.clone();

    info!(record_id, file = %key, path = %destination.display(), "downloading Zenodo archive");
    let bar = progress.bytes_bar(
        file.size,
        format!("Downloading {key} to {}", destination.display()),
    );
    client
        .download_record_file_by_key_to_path_with_progress(
            RecordId(record_id),
            &key,
            &destination,
            bar.clone(),
        )
        .await?;
    bar.finish_with_message(format!("Downloaded {key} to {}", destination.display()));
    Ok(destination)
}
