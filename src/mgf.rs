//! MGF enrichment logic built on the GitHub version of `mascot-rs`.

use std::path::{Path, PathBuf};

use mascot_rs::prelude::{MGFVec, SpectrumAlloc};
use serde::Serialize;

use crate::error::IoResultExt;
use crate::taxo::{TaxoDataset, TaxoRecord};
use crate::{Error, Result};

/// Summary statistics for one enriched MGF file.
#[derive(Debug, Clone, Serialize)]
pub struct EnrichmentStats {
    /// Source MGF path.
    pub input_path: PathBuf,
    /// Enriched MGF output path.
    pub output_path: PathBuf,
    /// Number of spectra parsed from the source file.
    pub spectra_total: usize,
    /// Number of spectra that received taxonomic metadata.
    pub spectra_enriched: usize,
    /// Number of spectra with no matching taxonomic row.
    pub spectra_without_taxo: usize,
    /// Number of enriched spectra without NCBI resolution.
    pub spectra_without_ncbi_resolution: usize,
    /// Number of most intense peaks kept per spectrum; zero means all peaks were kept.
    pub top_k_peaks: usize,
}

/// Enriches one MGF file with EMI taxonomy headers and writes a new MGF file.
pub fn enrich_mgf_file(
    input_path: &Path,
    output_dir: &Path,
    taxo: &TaxoDataset,
    top_k_peaks: usize,
    overwrite: bool,
) -> Result<EnrichmentStats> {
    std::fs::create_dir_all(output_dir).map_err(|source| Error::Io {
        path: output_dir.to_path_buf(),
        source,
    })?;
    let output_path = output_path_for(input_path, output_dir);
    if output_path.exists() && !overwrite {
        return Err(Error::OutputExists(output_path));
    }

    let mut spectra = load_mgf_vec(input_path)?;

    let mut stats = EnrichmentStats {
        input_path: input_path.to_path_buf(),
        output_path: output_path.clone(),
        spectra_total: spectra.len(),
        spectra_enriched: 0,
        spectra_without_taxo: 0,
        spectra_without_ncbi_resolution: 0,
        top_k_peaks,
    };

    for spectrum in spectra.iter_mut() {
        let feature_id = spectrum
            .feature_id()
            .or_else(|| spectrum.scans())
            .map(str::to_string);
        let Some(feature_id) = feature_id else {
            stats.spectra_without_taxo += 1;
            continue;
        };
        let Some(record) = taxo.get(&feature_id).or_else(|| taxo.default_record()) else {
            stats.spectra_without_taxo += 1;
            continue;
        };
        insert_emi_metadata(spectrum.metadata_mut(), record).map_err(|source| Error::Mascot {
            path: input_path.to_path_buf(),
            source,
        })?;
        stats.spectra_enriched += 1;
        if record.ncbi.is_none() {
            stats.spectra_without_ncbi_resolution += 1;
        }
    }

    let spectra = apply_top_k_peaks(spectra, top_k_peaks, input_path)?;
    spectra
        .to_path(&output_path)
        .map_err(|source| Error::Mascot {
            path: output_path.clone(),
            source,
        })?;
    Ok(stats)
}

/// Keeps only the top K most intense peaks per spectrum through mascot-rs.
fn apply_top_k_peaks(spectra: MGFVec, top_k_peaks: usize, input_path: &Path) -> Result<MGFVec> {
    if top_k_peaks == 0 {
        return Ok(spectra);
    }
    spectra
        .into_iter()
        .map(|spectrum| spectrum.top_k_peaks(top_k_peaks))
        .collect::<mascot_rs::prelude::Result<Vec<_>>>()
        .map(MGFVec::from)
        .map_err(|source| Error::Mascot {
            path: input_path.to_path_buf(),
            source,
        })
}

/// Loads an MGF file with dataset-specific tolerance for incomplete merged-scan headers.
fn load_mgf_vec(input_path: &Path) -> Result<MGFVec> {
    let content = std::fs::read_to_string(input_path).with_path(input_path)?;
    let filtered = filter_parseable_ms2_blocks(&content);
    MGFVec::try_from_iter(filtered.iter().map(String::as_str)).map_err(|source| Error::Mascot {
        path: input_path.to_path_buf(),
        source,
    })
}

/// Keeps MS2 ion blocks and removes headers that are incomplete for mascot-rs validation.
fn filter_parseable_ms2_blocks(content: &str) -> Vec<String> {
    let mut output = Vec::new();
    let mut block = Vec::new();
    let mut in_block = false;
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed == "BEGIN IONS" {
            block.clear();
            in_block = true;
        }
        if in_block {
            block.push(line.to_string());
        }
        if trimmed == "END IONS" && in_block {
            if block_is_ms2(&block) {
                output.extend(
                    block
                        .iter()
                        .filter(|line| {
                            let line = line.trim_start();
                            !line.starts_with("MERGED_STATS=") && !line.starts_with("MERGED_SCANS=")
                        })
                        .cloned(),
                );
            }
            block.clear();
            in_block = false;
        }
    }
    output
}

/// Returns whether an MGF ion block is an MS2 spectrum.
fn block_is_ms2(block: &[String]) -> bool {
    block
        .iter()
        .any(|line| line.trim().eq_ignore_ascii_case("MSLEVEL=2"))
}

/// Inserts EMI taxonomy metadata through mascot-rs metadata editing APIs.
fn insert_emi_metadata(
    metadata: &mut mascot_rs::prelude::MascotGenericFormatMetadata,
    record: &TaxoRecord,
) -> mascot_rs::prelude::Result<()> {
    metadata.insert_arbitrary_metadata("EMI_FEATURE_ID", &record.feature_id)?;
    if let Some(source_taxon_id) = record.source_taxon_id {
        metadata.insert_arbitrary_metadata("EMI_SOURCE_TAXON_ID", source_taxon_id.to_string())?;
    }
    if let Some(source_taxon_name) = &record.source_taxon_name {
        metadata.insert_arbitrary_metadata("EMI_SOURCE_TAXON_NAME", source_taxon_name)?;
    }
    if let Some(ncbi) = &record.ncbi {
        metadata.insert_arbitrary_metadata("EMI_TAXON_ID", ncbi.taxon_id.to_string())?;
        metadata.insert_arbitrary_metadata("EMI_TAXON_NAME", &ncbi.scientific_name)?;
        metadata.insert_arbitrary_metadata("EMI_TAXON_RANK", &ncbi.rank)?;
        metadata.insert_arbitrary_metadata(
            "EMI_TAXON_LINEAGE_IDS",
            ncbi.lineage_ids
                .iter()
                .map(u64::to_string)
                .collect::<Vec<_>>()
                .join("|"),
        )?;
        metadata
            .insert_arbitrary_metadata("EMI_TAXON_LINEAGE_NAMES", ncbi.lineage_names.join("|"))?;
    }
    Ok(())
}

/// Builds the output path for an enriched MGF file.
fn output_path_for(input_path: &Path, output_dir: &Path) -> PathBuf {
    let stem = input_path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("enriched");
    output_dir.join(format!("{stem}_emi_taxo_enriched.mgf"))
}

#[cfg(test)]
/// Tests for MGF output naming and mascot-rs metadata insertion.
mod tests {
    use super::*;
    use mascot_rs::prelude::Spectrum;

    use crate::taxo::TaxoRecord;
    use crate::taxonomy::ResolvedTaxon;

    #[test]
    /// Verifies enriched MGF file names keep the source stem.
    fn output_path_uses_enriched_suffix() {
        let path = output_path_for(
            Path::new("VGF159_A02_features_ms2_pos.mgf"),
            Path::new("out"),
        );
        assert_eq!(
            path,
            PathBuf::from("out/VGF159_A02_features_ms2_pos_emi_taxo_enriched.mgf")
        );
    }

    #[test]
    /// Verifies EMI headers are inserted through the mascot-rs metadata API.
    fn inserts_emi_metadata_with_mascot_api() {
        let document = concat!(
            "BEGIN IONS\n",
            "FEATURE_ID=feature-1\n",
            "PEPMASS=250.0\n",
            "CHARGE=1\n",
            "MSLEVEL=2\n",
            "100.0 10.0\n",
            "200.0 20.0\n",
            "END IONS\n",
        );
        let mut spectra: MGFVec = document.parse().expect("valid fixture");
        let record = TaxoRecord {
            feature_id: "feature-1".to_string(),
            source_taxon_id: Some(9606),
            source_taxon_name: Some("Homo sapiens".to_string()),
            ncbi: Some(ResolvedTaxon {
                taxon_id: 9606,
                scientific_name: "Homo sapiens".to_string(),
                rank: "species".to_string(),
                lineage_ids: vec![1, 9606],
                lineage_names: vec!["root".to_string(), "Homo sapiens".to_string()],
            }),
        };

        let spectrum = spectra.iter_mut().next().expect("one spectrum");
        insert_emi_metadata(spectrum.metadata_mut(), &record).expect("metadata insert");

        assert_eq!(
            spectra[0]
                .metadata()
                .arbitrary_metadata_value("EMI_TAXON_ID"),
            Some("9606")
        );
        assert_eq!(
            spectra[0]
                .metadata()
                .arbitrary_metadata_value("EMI_TAXON_LINEAGE_NAMES"),
            Some("root|Homo sapiens")
        );
    }

    #[test]
    /// Verifies top-K filtering is delegated to mascot-rs and preserves MGF metadata.
    fn keeps_top_k_peaks_with_mascot_api() {
        let document = concat!(
            "BEGIN IONS\n",
            "FEATURE_ID=feature-1\n",
            "PEPMASS=250.0\n",
            "CHARGE=1\n",
            "MSLEVEL=2\n",
            "100.0 10.0\n",
            "150.0 50.0\n",
            "200.0 20.0\n",
            "END IONS\n",
        );
        let spectra: MGFVec = document.parse().expect("valid fixture");

        let filtered =
            apply_top_k_peaks(spectra, 2, Path::new("fixture.mgf")).expect("top-k filtering");

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].len(), 2);
        assert_eq!(filtered[0].feature_id(), Some("feature-1"));
    }

    #[test]
    /// Verifies incomplete merged-scan statistics do not block dataset parsing.
    fn loads_mgf_with_merged_stats_without_merged_scans() {
        let temp_dir = tempfile::tempdir().expect("temp dir");
        let path = temp_dir.path().join("fixture.mgf");
        std::fs::write(
            &path,
            concat!(
                "BEGIN IONS\n",
                "FEATURE_ID=feature-1\n",
                "PEPMASS=250.0\n",
                "CHARGE=1\n",
                "MSLEVEL=2\n",
                "MERGED_STATS=1 / 2 (0 removed due to low quality, 1 removed due to low cosine).\n",
                "100.0 10.0\n",
                "200.0 20.0\n",
                "END IONS\n",
            ),
        )
        .expect("write fixture");

        let spectra = load_mgf_vec(&path).expect("parse lenient merged stats");

        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].feature_id(), Some("feature-1"));
    }

    #[test]
    /// Verifies incomplete merged-scan lists do not block dataset parsing.
    fn loads_mgf_with_merged_scans_without_merged_stats() {
        let temp_dir = tempfile::tempdir().expect("temp dir");
        let path = temp_dir.path().join("fixture.mgf");
        std::fs::write(
            &path,
            concat!(
                "BEGIN IONS\n",
                "FEATURE_ID=feature-1\n",
                "PEPMASS=250.0\n",
                "CHARGE=1\n",
                "MSLEVEL=2\n",
                "MERGED_SCANS=3362,3340,3402\n",
                "100.0 10.0\n",
                "200.0 20.0\n",
                "END IONS\n",
            ),
        )
        .expect("write fixture");

        let spectra = load_mgf_vec(&path).expect("parse lenient merged scans");

        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].feature_id(), Some("feature-1"));
    }

    #[test]
    /// Verifies correlated MS1 blocks are skipped before mascot-rs parsing.
    fn skips_ms1_blocks_before_parsing() {
        let temp_dir = tempfile::tempdir().expect("temp dir");
        let path = temp_dir.path().join("fixture.mgf");
        std::fs::write(
            &path,
            concat!(
                "BEGIN IONS\n",
                "FEATURE_ID=feature-1\n",
                "PEPMASS=593.2760620117188\n",
                "CHARGE=1\n",
                "MSLEVEL=1\n",
                "593.2760620117188 1.0E8\n",
                "592.2678833007812 1.3E6\n",
                "END IONS\n",
                "BEGIN IONS\n",
                "FEATURE_ID=feature-1\n",
                "PEPMASS=593.2760620117188\n",
                "CHARGE=1\n",
                "MSLEVEL=2\n",
                "100.0 10.0\n",
                "200.0 20.0\n",
                "END IONS\n",
            ),
        )
        .expect("write fixture");

        let spectra = load_mgf_vec(&path).expect("parse ms2 block");

        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].level(), Some(2));
    }
}
