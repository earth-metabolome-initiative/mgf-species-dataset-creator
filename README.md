# MGF Species Dataset Creator

Creates taxonomically enriched MGF files for the EMI PF1600 dataset published on
Zenodo record `10827917`.

The tool downloads and caches `pf1600_raw.tar.gz`, extracts the dataset, reads
the `taxo_output` metadata, resolves taxa against a local NCBI taxonomy dump,
parses MGF files with the GitHub version of `mascot-rs`, and writes MGF files
with additional EMI taxonomy headers.

## Requirements

- Rust nightly compatible with edition 2024.
- Enough disk space for the Zenodo archive and extraction. The archive is about
  11.9 GB.
- A local NCBI taxdump directory containing:
  - `nodes.dmp`
  - `names.dmp`
  - `merged.dmp` optional, but recommended

Download the NCBI taxdump from:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir -p taxdump
tar -xzf taxdump.tar.gz -C taxdump
```

## Run

Default run, using Zenodo record `10827917` and processing both positive and
negative MGF files:

```bash
cargo run -- run \
  --taxdump-dir taxdump \
  --output-dir output
```

This writes:

- Enriched MGF files ending in `_emi_taxo_enriched.mgf`, grouped by sample ID
- `output/pipeline_report.json`

The CLI prints each major stage, shows a byte progress bar while downloading
the Zenodo archive, and shows a file progress bar while enriching MGF files.

By default, downloads are cached under `.cache/zenodo` and extracted under
`.cache/extracted`.

## Use Existing Extracted Data

If the Zenodo archive is already extracted, skip download and extraction:

```bash
cargo run -- run \
  --extracted-dataset-dir /path/to/extracted/pf1600 \
  --taxdump-dir taxdump \
  --output-dir output
```

The extracted dataset must contain a `taxo_output` directory and matching MGF
files such as `*_features_ms2_pos.mgf` or `*_features_ms2_neg.mgf`.

## Options

```bash
cargo run -- run --help
```

Common options:

- `--record-id 10827917`: Zenodo record ID.
- `--archive-name pf1600_raw.tar.gz`: file to download from the Zenodo record.
- `--cache-dir .cache/zenodo`: archive cache directory.
- `--extraction-dir .cache/extracted`: extraction directory.
- `--extracted-dataset-dir <path>`: use existing extracted data.
- `--taxdump-dir <path>`: local NCBI taxdump directory.
- `--output-dir output`: output directory.
- `--polarity both`: one of `pos`, `neg`, or `both`.
- `--overwrite`: replace existing enriched MGF files.

## Added MGF Headers

For spectra with matching taxonomic metadata, the output MGF includes EMI fields
inserted through `mascot-rs` metadata APIs:

- `EMI_FEATURE_ID`
- `EMI_SOURCE_TAXON_ID`
- `EMI_SOURCE_TAXON_NAME`
- `EMI_TAXON_ID`
- `EMI_TAXON_NAME`
- `EMI_TAXON_RANK`
- `EMI_TAXON_LINEAGE_IDS`
- `EMI_TAXON_LINEAGE_NAMES`

## Verify

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```
