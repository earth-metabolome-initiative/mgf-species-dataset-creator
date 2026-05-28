//! Parsing of dataset-specific taxonomic metadata files.

use std::collections::HashMap;
use std::path::Path;

use serde::Serialize;
use serde_json::Value;

use crate::Error;
use crate::error::{IoResultExt, Result};
use crate::taxonomy::{NcbiTaxonomyResolver, ResolvedTaxon};

/// Feature-level taxonomic metadata joined to an MGF spectrum.
#[derive(Debug, Clone, Serialize)]
pub struct TaxoRecord {
    /// Feature identifier used to join against MGF `FEATURE_ID` or `SCANS`.
    pub feature_id: String,
    /// Taxon ID as provided by source metadata, when present.
    pub source_taxon_id: Option<u64>,
    /// Taxon name as provided by source metadata, when present.
    pub source_taxon_name: Option<String>,
    /// Resolved NCBI taxonomy entry, when resolution succeeded.
    pub ncbi: Option<ResolvedTaxon>,
}

/// In-memory feature-indexed taxonomic metadata.
#[derive(Debug, Clone)]
pub struct TaxoDataset {
    by_feature_id: HashMap<String, TaxoRecord>,
    default_record: Option<TaxoRecord>,
}

impl TaxoDataset {
    /// Loads all supported taxonomic metadata files from `taxo_output`.
    pub fn from_taxo_output(path: &Path, resolver: &NcbiTaxonomyResolver) -> Result<Self> {
        let mut records = HashMap::new();
        for entry in std::fs::read_dir(path).with_path(path)? {
            let entry = entry.with_path(path)?;
            let file = entry.path();
            match file.extension().and_then(|ext| ext.to_str()) {
                Some("tsv") => read_tsv(&file, resolver, &mut records)?,
                Some("json") => read_json(&file, resolver, &mut records)?,
                _ => {}
            }
        }
        if records.is_empty() {
            read_parent_sample_metadata(path, resolver, &mut records)?;
        }
        if records.is_empty() {
            return Err(Error::Taxo(format!(
                "no feature-level taxonomic records found under {}",
                path.display()
            )));
        }
        Ok(Self {
            default_record: (records.len() == 1)
                .then(|| records.values().next().cloned())
                .flatten(),
            by_feature_id: records,
        })
    }

    /// Returns the taxonomic record for a feature ID.
    pub fn get(&self, feature_id: &str) -> Option<&TaxoRecord> {
        self.by_feature_id.get(feature_id)
    }

    /// Returns the sample-level taxonomic record when exactly one record is loaded.
    pub fn default_record(&self) -> Option<&TaxoRecord> {
        self.default_record.as_ref()
    }

    /// Returns the number of indexed taxonomic records.
    pub fn len(&self) -> usize {
        self.by_feature_id.len()
    }

    /// Returns whether no taxonomic records are indexed.
    pub fn is_empty(&self) -> bool {
        self.by_feature_id.is_empty()
    }
}

/// Reads the parent sample metadata TSV as fallback when OTT output has no matches.
fn read_parent_sample_metadata(
    taxo_output_path: &Path,
    resolver: &NcbiTaxonomyResolver,
    records: &mut HashMap<String, TaxoRecord>,
) -> Result<()> {
    let Some(sample_root) = taxo_output_path.parent() else {
        return Ok(());
    };
    let Some(sample_id) = sample_root.file_name().and_then(|name| name.to_str()) else {
        return Ok(());
    };
    let metadata_path = sample_root.join(format!("{sample_id}_metadata.tsv"));
    if metadata_path.exists() {
        read_tsv(&metadata_path, resolver, records)?;
    }
    Ok(())
}

/// Reads a tab-separated taxonomic metadata file.
fn read_tsv(
    path: &Path,
    resolver: &NcbiTaxonomyResolver,
    records: &mut HashMap<String, TaxoRecord>,
) -> Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)
        .map_err(|source| Error::Csv {
            path: path.to_path_buf(),
            source,
        })?;
    let headers = reader
        .headers()
        .map_err(|source| Error::Csv {
            path: path.to_path_buf(),
            source,
        })?
        .clone();
    for row in reader.records() {
        let row = row.map_err(|source| Error::Csv {
            path: path.to_path_buf(),
            source,
        })?;
        let map = headers
            .iter()
            .zip(row.iter())
            .map(|(key, value)| (key.to_string(), Value::String(value.to_string())))
            .collect::<serde_json::Map<String, Value>>();
        if let Some(record) = record_from_map(&map, resolver)? {
            records.insert(record.feature_id.clone(), record);
        }
    }
    Ok(())
}

/// Reads a JSON taxonomic metadata file.
fn read_json(
    path: &Path,
    resolver: &NcbiTaxonomyResolver,
    records: &mut HashMap<String, TaxoRecord>,
) -> Result<()> {
    let value: Value = serde_json::from_str(&std::fs::read_to_string(path).with_path(path)?)
        .map_err(|source| Error::Json {
            path: path.to_path_buf(),
            source,
        })?;
    visit_json_objects(&value, resolver, records)
}

/// Walks JSON objects and arrays looking for feature-level taxonomic records.
fn visit_json_objects(
    value: &Value,
    resolver: &NcbiTaxonomyResolver,
    records: &mut HashMap<String, TaxoRecord>,
) -> Result<()> {
    match value {
        Value::Object(map) => {
            if let Some(record) = record_from_map(map, resolver)? {
                records.insert(record.feature_id.clone(), record);
            }
            for child in map.values() {
                visit_json_objects(child, resolver, records)?;
            }
        }
        Value::Array(values) => {
            for child in values {
                visit_json_objects(child, resolver, records)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Builds one taxonomic record from a flat JSON object.
fn record_from_map(
    map: &serde_json::Map<String, Value>,
    resolver: &NcbiTaxonomyResolver,
) -> Result<Option<TaxoRecord>> {
    let feature_id = first_string(
        map,
        &[
            "feature_id",
            "FEATURE_ID",
            "featureid",
            "feature",
            "sample_id",
            "sample_filename_pos",
            "sample_filename_neg",
            "cluster_id",
            "spectrum_id",
            "SCANS",
            "scans",
        ],
    );
    let Some(feature_id) = feature_id else {
        return Ok(None);
    };
    let source_taxon_id = first_u64(
        map,
        &[
            "taxon_id",
            "taxid",
            "tax_id",
            "ncbi_taxid",
            "ncbi_taxon_id",
            "NCBI_TAXID",
        ],
    );
    let source_taxon_name = first_string(
        map,
        &[
            "taxon_name",
            "taxon.name",
            "taxon.unique_name",
            "scientific_name",
            "species",
            "species_name",
            "organism_species",
            "organism",
            "name",
        ],
    );
    let source_taxon_id = source_taxon_id.or_else(|| ncbi_id_from_map(map));
    let ncbi = resolver.resolve(source_taxon_id, source_taxon_name.as_deref())?;
    Ok(Some(TaxoRecord {
        feature_id,
        source_taxon_id,
        source_taxon_name,
        ncbi,
    }))
}

/// Extracts an NCBI taxon ID from tax source fields such as `ncbi:110690`.
fn ncbi_id_from_map(map: &serde_json::Map<String, Value>) -> Option<u64> {
    ["taxon.tax_sources", "tax_sources"]
        .iter()
        .find_map(|key| map.get(*key).and_then(ncbi_id_from_value))
}

/// Extracts an NCBI taxon ID from JSON strings or arrays.
fn ncbi_id_from_value(value: &Value) -> Option<u64> {
    match value {
        Value::String(value) => ncbi_id_from_str(value),
        Value::Array(values) => values.iter().find_map(ncbi_id_from_value),
        _ => None,
    }
}

/// Extracts the first `ncbi:<id>` token from a string.
fn ncbi_id_from_str(value: &str) -> Option<u64> {
    value
        .split(|character: char| {
            !(character.is_ascii_alphanumeric() || character == ':' || character == '_')
        })
        .find_map(|token| token.strip_prefix("ncbi:")?.parse::<u64>().ok())
}

/// Returns the first non-empty string-like value for candidate keys.
fn first_string(map: &serde_json::Map<String, Value>, keys: &[&str]) -> Option<String> {
    keys.iter().find_map(|key| {
        map.get(*key).and_then(|value| match value {
            Value::String(value) => non_empty(value),
            Value::Number(value) => Some(value.to_string()),
            _ => None,
        })
    })
}

/// Returns the first unsigned integer value for candidate keys.
fn first_u64(map: &serde_json::Map<String, Value>, keys: &[&str]) -> Option<u64> {
    keys.iter().find_map(|key| {
        map.get(*key).and_then(|value| match value {
            Value::Number(value) => value.as_u64(),
            Value::String(value) => value.trim().parse::<u64>().ok(),
            _ => None,
        })
    })
}

/// Normalizes blank and missing-value marker strings.
fn non_empty(value: &str) -> Option<String> {
    let trimmed = value.trim();
    (!trimmed.is_empty() && trimmed != "NA" && trimmed != "N/A").then(|| trimmed.to_string())
}
