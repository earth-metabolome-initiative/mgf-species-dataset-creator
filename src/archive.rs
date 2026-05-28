//! Archive extraction and recursive dataset file discovery.

use std::fs::File;
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use tar::Archive;
use tracing::info;

use crate::error::{IoResultExt, Result};

/// Extracts a gzipped tar archive unless a previous extraction marker exists.
pub fn ensure_extracted(archive_path: &Path, extraction_dir: &Path) -> Result<PathBuf> {
    let marker = extraction_dir.join(".mgf-species-dataset-creator-extracted");
    if marker.exists() {
        return Ok(extraction_dir.to_path_buf());
    }

    std::fs::create_dir_all(extraction_dir).with_path(extraction_dir)?;
    info!(
        archive = %archive_path.display(),
        destination = %extraction_dir.display(),
        "extracting archive"
    );
    let archive_file = File::open(archive_path).with_path(archive_path)?;
    let decoder = GzDecoder::new(archive_file);
    let mut archive = Archive::new(decoder);
    archive.unpack(extraction_dir).with_path(extraction_dir)?;
    std::fs::write(&marker, b"ok\n").with_path(&marker)?;
    Ok(extraction_dir.to_path_buf())
}

/// Finds the first `taxo_output` directory under a dataset root.
pub fn find_taxo_output(root: &Path) -> Option<PathBuf> {
    find_dir_named(root, "taxo_output")
}

/// Finds all `taxo_output` directories under a dataset root.
pub fn find_taxo_outputs(root: &Path) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    collect_dirs_named(root, "taxo_output", &mut paths)?;
    paths.sort();
    Ok(paths)
}

/// Finds MGF files under a dataset root matching the requested polarity.
pub fn find_mgf_files(root: &Path, polarity: crate::Polarity) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    collect_files(root, &mut paths)?;
    paths.retain(|path| {
        let name = path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or_default();
        path.extension().and_then(|ext| ext.to_str()) == Some("mgf")
            && name.contains("_features_ms2_")
            && polarity.matches_path(path)
    });
    paths.sort();
    Ok(paths)
}

/// Recursively finds a directory with the requested basename.
fn find_dir_named(root: &Path, target: &str) -> Option<PathBuf> {
    let entries = std::fs::read_dir(root).ok()?;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            if path.file_name().and_then(|name| name.to_str()) == Some(target) {
                return Some(path);
            }
            if let Some(found) = find_dir_named(&path, target) {
                return Some(found);
            }
        }
    }
    None
}

/// Recursively collects directories with the requested basename.
fn collect_dirs_named(root: &Path, target: &str, out: &mut Vec<PathBuf>) -> Result<()> {
    for entry in std::fs::read_dir(root).with_path(root)? {
        let entry = entry.with_path(root)?;
        let path = entry.path();
        if path.is_dir() {
            if path.file_name().and_then(|name| name.to_str()) == Some(target) {
                out.push(path);
            } else {
                collect_dirs_named(&path, target, out)?;
            }
        }
    }
    Ok(())
}

/// Recursively collects all regular files under a root path.
fn collect_files(root: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
    for entry in std::fs::read_dir(root).with_path(root)? {
        let entry = entry.with_path(root)?;
        let path = entry.path();
        if path.is_dir() {
            collect_files(&path, out)?;
        } else {
            out.push(path);
        }
    }
    Ok(())
}
