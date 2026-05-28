//! NCBI taxdump parser and taxon resolver.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use serde::Serialize;

use crate::error::{IoResultExt, Result};
use crate::{Error, PipelineConfig};

/// NCBI taxon resolved to a scientific name, rank, and lineage.
#[derive(Debug, Clone, Serialize)]
pub struct ResolvedTaxon {
    /// Current NCBI taxon ID after applying merged-ID mappings.
    pub taxon_id: u64,
    /// Scientific name from `names.dmp`.
    pub scientific_name: String,
    /// Rank from `nodes.dmp`.
    pub rank: String,
    /// Root-to-taxon lineage IDs.
    pub lineage_ids: Vec<u64>,
    /// Root-to-taxon lineage scientific names.
    pub lineage_names: Vec<String>,
}

/// One parsed row from `nodes.dmp`.
#[derive(Debug, Clone)]
struct Node {
    /// Parent taxon ID.
    parent: u64,
    /// NCBI rank string.
    rank: String,
}

/// Map from taxon ID to scientific name.
type ScientificNames = HashMap<u64, String>;
/// Map from normalized taxon name or synonym to candidate taxon IDs.
type NameIndex = HashMap<String, Vec<u64>>;

/// Resolves taxon IDs or names using local NCBI taxdump files.
#[derive(Debug, Clone)]
pub struct NcbiTaxonomyResolver {
    nodes: HashMap<u64, Node>,
    names: ScientificNames,
    name_to_ids: NameIndex,
    merged: HashMap<u64, u64>,
}

impl NcbiTaxonomyResolver {
    /// Creates a resolver from the taxdump directory in the pipeline config.
    pub fn from_config(config: &PipelineConfig) -> Result<Self> {
        Self::from_taxdump_dir(&config.taxdump_dir)
    }

    /// Loads `nodes.dmp`, `names.dmp`, and optional `merged.dmp`.
    pub fn from_taxdump_dir(path: &Path) -> Result<Self> {
        let nodes = parse_nodes(&path.join("nodes.dmp"))?;
        let (names, name_to_ids) = parse_names(&path.join("names.dmp"))?;
        let merged = parse_merged(&path.join("merged.dmp"))?;
        Ok(Self {
            nodes,
            names,
            name_to_ids,
            merged,
        })
    }

    /// Resolves by taxon ID when present, otherwise by unique taxon name.
    pub fn resolve(
        &self,
        taxon_id: Option<u64>,
        name: Option<&str>,
    ) -> Result<Option<ResolvedTaxon>> {
        if let Some(taxon_id) = taxon_id {
            return self.resolve_id(taxon_id).map(Some);
        }
        let Some(name) = name.map(str::trim).filter(|name| !name.is_empty()) else {
            return Ok(None);
        };
        let key = normalize_name(name);
        let Some(ids) = self.name_to_ids.get(&key) else {
            return Ok(None);
        };
        if ids.len() != 1 {
            return Err(Error::Taxonomy(format!(
                "taxon name `{name}` resolves to multiple NCBI taxon IDs: {ids:?}"
            )));
        }
        self.resolve_id(ids[0]).map(Some)
    }

    /// Resolves one NCBI taxon ID and expands its lineage.
    pub fn resolve_id(&self, taxon_id: u64) -> Result<ResolvedTaxon> {
        let taxon_id = self.merged.get(&taxon_id).copied().unwrap_or(taxon_id);
        let node = self.nodes.get(&taxon_id).ok_or_else(|| {
            Error::Taxonomy(format!("NCBI taxon ID {taxon_id} is not in nodes.dmp"))
        })?;
        let scientific_name = self
            .names
            .get(&taxon_id)
            .cloned()
            .unwrap_or_else(|| taxon_id.to_string());
        let mut lineage_ids = Vec::new();
        let mut current = taxon_id;
        let mut seen = HashSet::new();
        loop {
            if !seen.insert(current) {
                return Err(Error::Taxonomy(format!(
                    "cycle detected while resolving NCBI lineage for {taxon_id}"
                )));
            }
            lineage_ids.push(current);
            let Some(node) = self.nodes.get(&current) else {
                break;
            };
            if node.parent == current {
                break;
            }
            current = node.parent;
        }
        lineage_ids.reverse();
        let lineage_names = lineage_ids
            .iter()
            .map(|id| {
                self.names
                    .get(id)
                    .cloned()
                    .unwrap_or_else(|| id.to_string())
            })
            .collect();
        Ok(ResolvedTaxon {
            taxon_id,
            scientific_name,
            rank: node.rank.clone(),
            lineage_ids,
            lineage_names,
        })
    }
}

/// Parses NCBI `nodes.dmp`.
fn parse_nodes(path: &Path) -> Result<HashMap<u64, Node>> {
    let content = std::fs::read_to_string(path).with_path(path)?;
    let mut nodes = HashMap::new();
    for line in content.lines().filter(|line| !line.trim().is_empty()) {
        let parts = taxdump_parts(line);
        if parts.len() < 3 {
            continue;
        }
        let id = parse_u64(parts[0], path)?;
        nodes.insert(
            id,
            Node {
                parent: parse_u64(parts[1], path)?,
                rank: parts[2].to_string(),
            },
        );
    }
    Ok(nodes)
}

/// Parses NCBI `names.dmp`.
fn parse_names(path: &Path) -> Result<(ScientificNames, NameIndex)> {
    let content = std::fs::read_to_string(path).with_path(path)?;
    let mut scientific = HashMap::new();
    let mut name_to_ids: HashMap<String, Vec<u64>> = HashMap::new();
    for line in content.lines().filter(|line| !line.trim().is_empty()) {
        let parts = taxdump_parts(line);
        if parts.len() < 4 {
            continue;
        }
        let id = parse_u64(parts[0], path)?;
        let name = parts[1].to_string();
        name_to_ids
            .entry(normalize_name(&name))
            .or_default()
            .push(id);
        if parts[3] == "scientific name" {
            scientific.insert(id, name);
        }
    }
    Ok((scientific, name_to_ids))
}

/// Parses optional NCBI `merged.dmp`.
fn parse_merged(path: &Path) -> Result<HashMap<u64, u64>> {
    if !path.exists() {
        return Ok(HashMap::new());
    }
    let content = std::fs::read_to_string(path).with_path(path)?;
    let mut merged = HashMap::new();
    for line in content.lines().filter(|line| !line.trim().is_empty()) {
        let parts = taxdump_parts(line);
        if parts.len() >= 2 {
            merged.insert(parse_u64(parts[0], path)?, parse_u64(parts[1], path)?);
        }
    }
    Ok(merged)
}

/// Splits one taxdump line into trimmed fields.
fn taxdump_parts(line: &str) -> Vec<&str> {
    line.split('|').map(str::trim).collect()
}

/// Parses an integer field with taxdump path context.
fn parse_u64(value: &str, path: &Path) -> Result<u64> {
    value.parse::<u64>().map_err(|source| {
        Error::Taxonomy(format!(
            "could not parse `{value}` as integer in {}: {source}",
            PathBuf::from(path).display()
        ))
    })
}

/// Normalizes a taxon name for exact case-insensitive lookup.
fn normalize_name(name: &str) -> String {
    name.trim().to_ascii_lowercase()
}

#[cfg(test)]
/// Tests for NCBI taxonomy resolution.
mod tests {
    use super::*;

    #[test]
    /// Verifies merged taxon IDs resolve to the current ID and lineage.
    fn resolves_merged_taxon_and_lineage() {
        let resolver = NcbiTaxonomyResolver {
            nodes: HashMap::from([
                (
                    1,
                    Node {
                        parent: 1,
                        rank: "no rank".to_string(),
                    },
                ),
                (
                    2,
                    Node {
                        parent: 1,
                        rank: "species".to_string(),
                    },
                ),
            ]),
            names: HashMap::from([(1, "root".to_string()), (2, "Test species".to_string())]),
            name_to_ids: HashMap::from([("test species".to_string(), vec![2])]),
            merged: HashMap::from([(20, 2)]),
        };

        let resolved = resolver.resolve(Some(20), None).expect("resolve").unwrap();

        assert_eq!(resolved.taxon_id, 2);
        assert_eq!(resolved.lineage_ids, vec![1, 2]);
        assert_eq!(resolved.lineage_names, vec!["root", "Test species"]);
    }
}
