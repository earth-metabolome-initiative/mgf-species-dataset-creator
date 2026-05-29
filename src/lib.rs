//! Tools for creating taxonomically enriched MGF datasets from Zenodo records.

#![warn(missing_docs)]

/// Archive extraction and dataset file discovery helpers.
pub mod archive;
/// Pipeline configuration types.
pub mod config;
/// Error and result types used across the crate.
pub mod error;
/// MGF parsing, enrichment, and writing helpers.
pub mod mgf;
/// End-to-end pipeline orchestration.
pub mod pipeline;
/// Terminal progress reporting helpers.
pub mod progress;
/// Utilities for repairing and publishing corrected PF1600 archives.
pub mod repair;
/// Dataset-specific taxonomic metadata parsing.
pub mod taxo;
/// NCBI taxonomy dump parsing and resolution.
pub mod taxonomy;
/// Zenodo download and cache helpers.
pub mod zenodo;

pub use config::{PipelineConfig, Polarity};
pub use error::{Error, Result};
pub use pipeline::{PipelineReport, ProcessedMgf, run_pipeline};

/// Default Zenodo record ID for the PF1600 raw dataset.
pub const DEFAULT_RECORD_ID: u64 = 10_827_917;
/// Default archive file name in the Zenodo record.
pub const DEFAULT_ZENODO_ARCHIVE: &str = "pf1600_raw.tar.gz";
