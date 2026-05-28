//! Error and result types for the dataset enrichment pipeline.

use std::path::PathBuf;

/// Error type returned by this crate.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Filesystem I/O failed for a path.
    #[error("I/O error at {path}: {source}")]
    Io {
        /// Path involved in the failing I/O operation.
        path: PathBuf,
        /// Underlying I/O error.
        #[source]
        source: std::io::Error,
    },
    /// CSV or TSV parsing failed for a path.
    #[error("CSV error at {path}: {source}")]
    Csv {
        /// Path of the CSV or TSV file.
        path: PathBuf,
        /// Underlying CSV parser error.
        #[source]
        source: csv::Error,
    },
    /// JSON parsing or serialization failed for a path.
    #[error("JSON error at {path}: {source}")]
    Json {
        /// Path of the JSON file or report.
        path: PathBuf,
        /// Underlying JSON error.
        #[source]
        source: serde_json::Error,
    },
    /// Zenodo API or download operation failed.
    #[error("Zenodo error: {0}")]
    Zenodo(#[from] zenodo_rs::ZenodoError),
    /// MGF parsing, editing, or writing failed through mascot-rs.
    #[error("Mascot MGF error at {path}: {source}")]
    Mascot {
        /// Path of the MGF file being parsed or written.
        path: PathBuf,
        /// Underlying mascot-rs error.
        #[source]
        source: mascot_rs::prelude::MascotError,
    },
    /// A required input path was not found.
    #[error("missing required path: {0}")]
    MissingPath(PathBuf),
    /// The configured Zenodo file was not present in the record.
    #[error("missing Zenodo file `{file}` in record {record_id}")]
    MissingZenodoFile {
        /// Zenodo record ID that was queried.
        record_id: u64,
        /// File key expected in the Zenodo record.
        file: String,
    },
    /// No MGF files matched the selected polarity.
    #[error("no MGF files matched polarity `{0}` under {1}")]
    NoMgfFiles(String, PathBuf),
    /// NCBI taxdump parsing or taxon resolution failed.
    #[error("NCBI taxonomy error: {0}")]
    Taxonomy(String),
    /// Dataset taxonomic metadata was missing or inconsistent.
    #[error("taxonomic metadata error: {0}")]
    Taxo(String),
    /// An output path already exists and overwrite was disabled.
    #[error("output already exists: {0}; pass --overwrite to replace it")]
    OutputExists(PathBuf),
}

/// Crate-local result alias.
pub type Result<T> = std::result::Result<T, Error>;

/// Adds path context to standard I/O results.
pub(crate) trait IoResultExt<T> {
    /// Converts an I/O result into the crate error type with path context.
    fn with_path(self, path: impl Into<PathBuf>) -> Result<T>;
}

impl<T> IoResultExt<T> for std::io::Result<T> {
    /// Converts an I/O result into the crate error type with path context.
    fn with_path(self, path: impl Into<PathBuf>) -> Result<T> {
        let path = path.into();
        self.map_err(|source| Error::Io { path, source })
    }
}
