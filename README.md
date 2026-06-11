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

- `output/emi_taxo_enriched_pos.mgf` when positive-mode files are processed
- `output/emi_taxo_enriched_neg.mgf` when negative-mode files are processed
- `output/pipeline_report.json`

The CLI prints each major stage, shows a byte progress bar while downloading
the Zenodo archive, and shows a file progress bar while enriching MGF files.
Output spectra keep the top 128 most intense peaks by default, using
`mascot-rs` peak editing APIs, and peak intensities are normalized by default.

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

## Plot Input Distributions

To inspect the raw input MGF dataset before enrichment:

```bash
cargo run -- plot-input-stats \
  --dataset-dir .cache/extracted/pf1600_raw \
  --output-dir plots
```

This writes SVG histograms split by positive and negative mode:

- `plots/peak_count_distribution.svg`: number of peaks per spectrum.
- `plots/feature_count_distribution.svg`: number of features per MGF file.
- `plots/parent_mz_distribution.svg`: precursor/parent m/z distribution.

It also writes `plots/input_mgf_file_stats.csv`,
`plots/input_histograms.csv`, and `plots/input_stats_report.json`.

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
- `--top-k-peaks 128`: keep the top K most intense peaks per spectrum; use `0`
  to keep all peaks.
- `--no-normalize-peak-intensities`: keep original peak intensities instead of
  mascot-rs normalized intensities.
- `--overwrite`: replace existing enriched MGF files.

Plotting options are available with:

```bash
cargo run -- plot-input-stats --help
```

## Restore Legacy Sirius Artifacts

If an interrupted repair run leaves the legacy Sirius artifact files in a bad
state, restore only those two files from the original Zenodo archive:

```bash
cargo run -- restore-legacy-sirius \
  --original-archive-path .cache/zenodo/pf1600_raw.tar.gz \
  --work-dataset-dir .cache/corrected/pf1600_raw
```

After these files have been restored in the extracted source dataset, the full
`repair-dataset --overwrite` command preserves them as `*_sirius_pos.mgf`
before downloading the corrected MassIVE feature MS2 files.

When `repair-dataset --publish` is used, the command publishes a new version of
the existing Zenodo record family, defaulting to record `10827917`; it does not
create a separate Zenodo dataset. Use `--record-id <id>` to target a different
existing record.

If the corrected archive already exists and should not be rebuilt, publish only
that archive:

```bash
cargo run -- publish-archive \
  --archive-path .cache/corrected/pf1600_raw_corrected.tar.gz \
  --record-id 10827917
```

The Zenodo token is read from the process environment first, then from `.env`.
Create `.env` locally with:

```bash
ZENODO_TOKEN=your-token
```

Zenodo sandbox is separate from production Zenodo, so production record
`10827917` does not exist there. To test the upload workflow on sandbox, create
a separate sandbox dataset:

```bash
cargo run -- publish-archive \
  --archive-path .cache/corrected/pf1600_raw_corrected.tar.gz \
  --sandbox \
  --new-dataset
```

For the final production release, omit `--sandbox` and `--new-dataset` so the
archive is published as a new version of record `10827917`.

## Added MGF Headers

The output MGF includes selected sample-level metadata from
`<sample_id>_metadata.tsv`:

- `sample_id`
- `sample_type`
- `sample_plate_id`
- `bio_leish_donovani_10ugml_inhibition`
- `bio_leish_donovani_2ugml_inhibition`
- `bio_tryp_brucei_rhodesiense_10ugml_inhibition`
- `bio_tryp_brucei_rhodesiense_2ugml_inhibition`
- `bio_tryp_cruzi_10ugml_inhibition`
- `bio_l6_cytotoxicity_10ugml_inhibition`

Positive-mode spectra also include:

- `sample_filename_pos`
- `pos_injection_date`

Negative-mode spectra also include:

- `sample_filename_neg`
- `neg_injection_date`

For spectra with matching NCBI-resolved taxonomic metadata, the output MGF also
includes the minimal EMI taxonomy field inserted through `mascot-rs` metadata
APIs:

- `EMI_TAXON_ID`

## Verify

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```
