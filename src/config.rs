//! Configuration types for selecting inputs, outputs, and MGF polarity.

use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;

use crate::{DEFAULT_RECORD_ID, DEFAULT_ZENODO_ARCHIVE};

/// MGF ionization polarity selection for dataset processing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polarity {
    /// Process positive-mode MGF files.
    Pos,
    /// Process negative-mode MGF files.
    Neg,
    /// Process both positive- and negative-mode MGF files.
    Both,
}

impl Polarity {
    /// Returns whether a path matches this polarity selector.
    pub fn matches_path(self, path: &std::path::Path) -> bool {
        let name = path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or_default();
        match self {
            Self::Pos => name.ends_with("_pos.mgf") || name.contains("_ms2_pos"),
            Self::Neg => name.ends_with("_neg.mgf") || name.contains("_ms2_neg"),
            Self::Both => {
                name.ends_with("_pos.mgf")
                    || name.contains("_ms2_pos")
                    || name.ends_with("_neg.mgf")
                    || name.contains("_ms2_neg")
            }
        }
    }

    /// Returns the stable CLI string for this polarity.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pos => "pos",
            Self::Neg => "neg",
            Self::Both => "both",
        }
    }
}

impl fmt::Display for Polarity {
    /// Formats the polarity as its stable CLI string.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for Polarity {
    /// Error message returned for invalid polarity strings.
    type Err = String;

    /// Parses a polarity selector from a CLI string.
    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value {
            "pos" | "positive" => Ok(Self::Pos),
            "neg" | "negative" => Ok(Self::Neg),
            "both" => Ok(Self::Both),
            other => Err(format!("expected pos, neg, or both, got `{other}`")),
        }
    }
}

/// Configuration for one full dataset enrichment run.
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    /// Zenodo record ID to fetch.
    pub record_id: u64,
    /// Archive file name inside the Zenodo record.
    pub archive_name: String,
    /// Directory where the Zenodo archive is cached.
    pub cache_dir: PathBuf,
    /// Directory where the archive is extracted.
    pub extraction_dir: PathBuf,
    /// Optional existing extracted dataset root.
    pub extracted_dataset_dir: Option<PathBuf>,
    /// Directory containing NCBI taxdump files.
    pub taxdump_dir: PathBuf,
    /// Directory where enriched MGF files and reports are written.
    pub output_dir: PathBuf,
    /// MGF polarity selection.
    pub polarity: Polarity,
    /// Whether existing output files may be replaced.
    pub overwrite: bool,
}

impl Default for PipelineConfig {
    /// Returns the default PF1600 Zenodo pipeline configuration.
    fn default() -> Self {
        Self {
            record_id: DEFAULT_RECORD_ID,
            archive_name: DEFAULT_ZENODO_ARCHIVE.to_string(),
            cache_dir: PathBuf::from(".cache/zenodo"),
            extraction_dir: PathBuf::from(".cache/extracted"),
            extracted_dataset_dir: None,
            taxdump_dir: PathBuf::from("taxdump"),
            output_dir: PathBuf::from("output"),
            polarity: Polarity::Both,
            overwrite: false,
        }
    }
}
