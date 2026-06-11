//! Input dataset statistics and lightweight SVG plotting utilities.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use mascot_rs::prelude::Spectrum;
use serde::Serialize;

use crate::Error;
use crate::archive::find_mgf_files;
use crate::config::Polarity;
use crate::error::{IoResultExt, Result};
use crate::mgf::load_mgf_vec;
use crate::progress::ProgressReporter;

/// Configuration for plotting input MGF dataset distributions.
#[derive(Debug, Clone)]
pub struct PlotConfig {
    /// Root directory containing extracted PF1600 sample directories.
    pub dataset_dir: PathBuf,
    /// Directory where SVG plots and statistics files are written.
    pub output_dir: PathBuf,
    /// MGF polarity selector.
    pub polarity: Polarity,
    /// Number of histogram bins.
    pub bins: usize,
}

impl Default for PlotConfig {
    /// Returns defaults targeting the standard extracted PF1600 dataset.
    fn default() -> Self {
        Self {
            dataset_dir: PathBuf::from(".cache/extracted/pf1600_raw"),
            output_dir: PathBuf::from("plots"),
            polarity: Polarity::Both,
            bins: 60,
        }
    }
}

/// Report produced after plotting input MGF statistics.
#[derive(Debug, Clone, Serialize)]
pub struct PlotReport {
    /// Dataset root used for discovery.
    pub dataset_dir: PathBuf,
    /// Number of MGF files parsed.
    pub mgf_files: usize,
    /// Number of spectra parsed.
    pub spectra: usize,
    /// Number of histogram bins used in each plot.
    pub bins: usize,
    /// Per-polarity summary statistics.
    pub summaries: Vec<PolaritySummary>,
    /// SVG plot paths written by the command.
    pub plots: Vec<PathBuf>,
}

/// Summary statistics for one metric and polarity.
#[derive(Debug, Clone, Serialize)]
pub struct PolaritySummary {
    /// Polarity label, either `pos` or `neg`.
    pub polarity: String,
    /// Metric label.
    pub metric: String,
    /// Number of values summarized.
    pub count: usize,
    /// Minimum value.
    pub min: f64,
    /// First quartile.
    pub q1: f64,
    /// Median value.
    pub median: f64,
    /// Arithmetic mean.
    pub mean: f64,
    /// Third quartile.
    pub q3: f64,
    /// Maximum value.
    pub max: f64,
}

/// Runs input dataset statistics collection and writes plots.
pub fn plot_input_stats(config: PlotConfig, progress: &ProgressReporter) -> Result<PlotReport> {
    if config.bins == 0 {
        return Err(Error::MissingValue("bins must be at least 1".to_string()));
    }
    if !config.dataset_dir.exists() {
        return Err(Error::MissingPath(config.dataset_dir));
    }

    std::fs::create_dir_all(&config.output_dir).with_path(&config.output_dir)?;
    let spinner = progress.spinner(format!(
        "Finding {} MGF files under {}",
        config.polarity,
        config.dataset_dir.display()
    ));
    let mgf_files = find_mgf_files(&config.dataset_dir, config.polarity)?;
    spinner.finish_with_message(format!("Found {} MGF file(s)", mgf_files.len()));

    let bar = progress.bar(mgf_files.len() as u64, "Collecting input MGF statistics");
    let mut stats = PlotStats::default();
    let mut file_rows = Vec::with_capacity(mgf_files.len());

    for path in &mgf_files {
        bar.set_message(format!("Parsing {}", ProgressReporter::file_label(path)));
        let polarity = polarity_for_path(path);
        let spectra = load_mgf_vec(path)?;
        let feature_count = spectra.len();
        let mut total_peaks = 0usize;

        for spectrum in &spectra {
            let peak_count = spectrum.len();
            let parent_mz = spectrum.precursor_mz();
            total_peaks += peak_count;
            stats.push_spectrum(polarity, peak_count, parent_mz);
        }

        stats.push_feature_count(polarity, feature_count);
        file_rows.push(MgfFileStats {
            sample_id: sample_id_for_path(path),
            polarity: polarity.to_string(),
            path: path.clone(),
            features: feature_count,
            total_peaks,
            mean_peaks_per_feature: mean_ratio(total_peaks, feature_count),
        });
        bar.inc(1);
    }
    bar.finish_with_message("Input MGF statistics collected");

    let summaries = stats.summaries();
    let histogram_rows = build_histogram_rows(&stats, config.bins);
    let plots = write_plots(&config.output_dir, &stats, config.bins)?;
    write_file_stats_csv(
        &config.output_dir.join("input_mgf_file_stats.csv"),
        &file_rows,
    )?;
    write_histogram_csv(
        &config.output_dir.join("input_histograms.csv"),
        &histogram_rows,
    )?;

    let report = PlotReport {
        dataset_dir: config.dataset_dir,
        mgf_files: file_rows.len(),
        spectra: stats.spectrum_count(),
        bins: config.bins,
        summaries,
        plots,
    };
    write_report(&config.output_dir.join("input_stats_report.json"), &report)?;
    Ok(report)
}

/// Internal accumulator for distributions split by polarity.
#[derive(Debug, Default)]
struct PlotStats {
    peak_counts: MetricValues,
    feature_counts: MetricValues,
    parent_masses: MetricValues,
}

impl PlotStats {
    /// Adds one spectrum-level observation.
    fn push_spectrum(&mut self, polarity: PlotPolarity, peak_count: usize, parent_mz: f64) {
        self.peak_counts.push(polarity, peak_count as f64);
        self.parent_masses.push(polarity, parent_mz);
    }

    /// Adds one file-level feature count observation.
    fn push_feature_count(&mut self, polarity: PlotPolarity, features: usize) {
        self.feature_counts.push(polarity, features as f64);
    }

    /// Returns the total number of spectrum observations.
    fn spectrum_count(&self) -> usize {
        self.peak_counts.len()
    }

    /// Builds summary statistics for all metrics and polarities.
    fn summaries(&self) -> Vec<PolaritySummary> {
        let mut summaries = Vec::new();
        for (metric, values) in [
            ("peak_count", &self.peak_counts),
            ("feature_count", &self.feature_counts),
            ("parent_mz", &self.parent_masses),
        ] {
            for polarity in [PlotPolarity::Pos, PlotPolarity::Neg] {
                if let Some(summary) = summarize(metric, polarity, values.values(polarity)) {
                    summaries.push(summary);
                }
            }
        }
        summaries
    }
}

/// Values for one metric split into positive and negative polarity.
#[derive(Debug, Default)]
struct MetricValues {
    pos: Vec<f64>,
    neg: Vec<f64>,
}

impl MetricValues {
    /// Adds one value for the selected polarity.
    fn push(&mut self, polarity: PlotPolarity, value: f64) {
        match polarity {
            PlotPolarity::Pos => self.pos.push(value),
            PlotPolarity::Neg => self.neg.push(value),
        }
    }

    /// Returns values for the selected polarity.
    fn values(&self, polarity: PlotPolarity) -> &[f64] {
        match polarity {
            PlotPolarity::Pos => &self.pos,
            PlotPolarity::Neg => &self.neg,
        }
    }

    /// Returns values for both polarities.
    fn all_values(&self) -> impl Iterator<Item = f64> + '_ {
        self.pos.iter().chain(self.neg.iter()).copied()
    }

    /// Returns the number of values across polarities.
    fn len(&self) -> usize {
        self.pos.len() + self.neg.len()
    }
}

/// Binary polarity labels used in plots.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum PlotPolarity {
    /// Positive ionization mode.
    Pos,
    /// Negative ionization mode.
    Neg,
}

impl std::fmt::Display for PlotPolarity {
    /// Formats the polarity for output tables.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Pos => f.write_str("pos"),
            Self::Neg => f.write_str("neg"),
        }
    }
}

/// One row of per-MGF file statistics.
#[derive(Debug)]
struct MgfFileStats {
    sample_id: String,
    polarity: String,
    path: PathBuf,
    features: usize,
    total_peaks: usize,
    mean_peaks_per_feature: f64,
}

/// One histogram bin row.
#[derive(Debug)]
struct HistogramRow {
    metric: String,
    polarity: String,
    bin_start: f64,
    bin_end: f64,
    count: usize,
}

/// Infers positive or negative mode from an MGF path.
fn polarity_for_path(path: &Path) -> PlotPolarity {
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default();
    if name.contains("_ms2_neg") || name.ends_with("_neg.mgf") {
        PlotPolarity::Neg
    } else {
        PlotPolarity::Pos
    }
}

/// Extracts a sample identifier from a dataset MGF path.
fn sample_id_for_path(path: &Path) -> String {
    path.parent()
        .and_then(Path::parent)
        .and_then(Path::file_name)
        .and_then(|name| name.to_str())
        .map(str::to_string)
        .or_else(|| {
            path.file_stem()
                .and_then(|name| name.to_str())
                .map(str::to_string)
        })
        .unwrap_or_else(|| "unknown".to_string())
}

/// Returns `numerator / denominator`, or zero when the denominator is zero.
fn mean_ratio(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

/// Computes summary statistics for one metric/polarity pair.
fn summarize(metric: &str, polarity: PlotPolarity, values: &[f64]) -> Option<PolaritySummary> {
    if values.is_empty() {
        return None;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let sum = sorted.iter().sum::<f64>();
    Some(PolaritySummary {
        polarity: polarity.to_string(),
        metric: metric.to_string(),
        count: sorted.len(),
        min: sorted[0],
        q1: percentile(&sorted, 0.25),
        median: percentile(&sorted, 0.5),
        mean: sum / sorted.len() as f64,
        q3: percentile(&sorted, 0.75),
        max: sorted[sorted.len() - 1],
    })
}

/// Computes a nearest-rank percentile with linear interpolation.
fn percentile(sorted: &[f64], percentile: f64) -> f64 {
    if sorted.len() == 1 {
        return sorted[0];
    }
    let position = percentile * (sorted.len() - 1) as f64;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    if lower == upper {
        sorted[lower]
    } else {
        let fraction = position - lower as f64;
        sorted[lower] * (1.0 - fraction) + sorted[upper] * fraction
    }
}

/// Builds CSV histogram rows for every metric and polarity.
fn build_histogram_rows(stats: &PlotStats, bins: usize) -> Vec<HistogramRow> {
    let mut rows = Vec::new();
    for (metric, values) in [
        ("peak_count", &stats.peak_counts),
        ("feature_count", &stats.feature_counts),
        ("parent_mz", &stats.parent_masses),
    ] {
        let histogram = Histogram::from_metric(values, bins);
        for polarity in [PlotPolarity::Pos, PlotPolarity::Neg] {
            let counts = histogram.counts(polarity);
            for (index, count) in counts.iter().copied().enumerate() {
                rows.push(HistogramRow {
                    metric: metric.to_string(),
                    polarity: polarity.to_string(),
                    bin_start: histogram.bin_start(index),
                    bin_end: histogram.bin_end(index),
                    count,
                });
            }
        }
    }
    rows
}

/// Histogram counts for a metric split by polarity.
#[derive(Debug)]
struct Histogram {
    min: f64,
    max: f64,
    bins: usize,
    counts: BTreeMap<PlotPolarity, Vec<usize>>,
}

impl Histogram {
    /// Creates a histogram over the combined positive and negative metric range.
    fn from_metric(values: &MetricValues, bins: usize) -> Self {
        let (mut min, mut max) = min_max(values.all_values()).unwrap_or((0.0, 1.0));
        if min.total_cmp(&max).is_eq() {
            min -= 0.5;
            max += 0.5;
        }
        let mut counts = BTreeMap::new();
        counts.insert(
            PlotPolarity::Pos,
            bin_values(values.values(PlotPolarity::Pos), min, max, bins),
        );
        counts.insert(
            PlotPolarity::Neg,
            bin_values(values.values(PlotPolarity::Neg), min, max, bins),
        );
        Self {
            min,
            max,
            bins,
            counts,
        }
    }

    /// Returns counts for one polarity.
    fn counts(&self, polarity: PlotPolarity) -> &[usize] {
        self.counts.get(&polarity).map(Vec::as_slice).unwrap_or(&[])
    }

    /// Returns the lower bound for a bin.
    fn bin_start(&self, index: usize) -> f64 {
        self.min + index as f64 * self.width()
    }

    /// Returns the upper bound for a bin.
    fn bin_end(&self, index: usize) -> f64 {
        if index + 1 == self.bins {
            self.max
        } else {
            self.bin_start(index + 1)
        }
    }

    /// Returns the bin width.
    fn width(&self) -> f64 {
        (self.max - self.min) / self.bins as f64
    }

    /// Returns the highest bin count across both polarities.
    fn max_count(&self) -> usize {
        self.counts
            .values()
            .flat_map(|counts| counts.iter())
            .copied()
            .max()
            .unwrap_or(1)
    }
}

/// Computes the minimum and maximum finite values from an iterator.
fn min_max(values: impl Iterator<Item = f64>) -> Option<(f64, f64)> {
    values
        .filter(|value| value.is_finite())
        .fold(None, |acc, value| match acc {
            None => Some((value, value)),
            Some((min, max)) => Some((min.min(value), max.max(value))),
        })
}

/// Computes bin counts for one polarity.
fn bin_values(values: &[f64], min: f64, max: f64, bins: usize) -> Vec<usize> {
    let width = (max - min) / bins as f64;
    let mut counts = vec![0usize; bins];
    for value in values.iter().copied().filter(|value| value.is_finite()) {
        let mut index = ((value - min) / width).floor() as usize;
        if index >= bins {
            index = bins - 1;
        }
        counts[index] += 1;
    }
    counts
}

/// Writes all SVG plots and returns their paths.
fn write_plots(output_dir: &Path, stats: &PlotStats, bins: usize) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    for (metric, title, x_label, values) in [
        (
            "peak_count",
            "Peak Count per Spectrum",
            "Peaks per spectrum",
            &stats.peak_counts,
        ),
        (
            "feature_count",
            "Feature Count per MGF File",
            "Features per MGF file",
            &stats.feature_counts,
        ),
        (
            "parent_mz",
            "Parent Mass Distribution",
            "Precursor m/z",
            &stats.parent_masses,
        ),
    ] {
        let path = output_dir.join(format!("{metric}_distribution.svg"));
        write_histogram_svg(&path, title, x_label, values, bins)?;
        paths.push(path);
    }
    Ok(paths)
}

/// Writes one grouped positive/negative histogram as SVG.
fn write_histogram_svg(
    path: &Path,
    title: &str,
    x_label: &str,
    values: &MetricValues,
    bins: usize,
) -> Result<()> {
    let histogram = Histogram::from_metric(values, bins);
    let width = 960.0;
    let height = 540.0;
    let left = 80.0;
    let right = 30.0;
    let top = 56.0;
    let bottom = 78.0;
    let plot_width = width - left - right;
    let plot_height = height - top - bottom;
    let max_count = histogram.max_count().max(1) as f64;
    let group_width = plot_width / bins as f64;
    let bar_width = (group_width * 0.38).max(1.0);

    let mut svg = String::new();
    svg.push_str(r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    svg.push('\n');
    svg.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"#
    ));
    svg.push_str(r##"<rect width="100%" height="100%" fill="#ffffff"/>"##);
    svg.push_str(&format!(
        r#"<text x="480" y="30" text-anchor="middle" font-family="sans-serif" font-size="22" font-weight="700">{}</text>"#,
        escape_xml(title)
    ));
    let area = PlotArea {
        left,
        top,
        width: plot_width,
        height: plot_height,
    };
    draw_axes(&mut svg, area, max_count, x_label, &histogram);

    for (offset, polarity, color) in [
        (0.08, PlotPolarity::Pos, "#2563eb"),
        (0.54, PlotPolarity::Neg, "#dc2626"),
    ] {
        for (index, count) in histogram.counts(polarity).iter().copied().enumerate() {
            let x = left + index as f64 * group_width + group_width * offset;
            let bar_height = (count as f64 / max_count) * plot_height;
            let y = top + plot_height - bar_height;
            svg.push_str(&format!(
                r#"<rect x="{x:.2}" y="{y:.2}" width="{bar_width:.2}" height="{bar_height:.2}" fill="{color}" opacity="0.78"/>"#
            ));
        }
    }

    svg.push_str(
        r##"<rect x="730" y="58" width="14" height="14" fill="#2563eb" opacity="0.78"/>"##,
    );
    svg.push_str(r#"<text x="752" y="70" font-family="sans-serif" font-size="14">positive</text>"#);
    svg.push_str(
        r##"<rect x="830" y="58" width="14" height="14" fill="#dc2626" opacity="0.78"/>"##,
    );
    svg.push_str(r#"<text x="852" y="70" font-family="sans-serif" font-size="14">negative</text>"#);
    svg.push_str("</svg>\n");

    std::fs::write(path, svg).with_path(path)
}

/// Draws SVG axes, grid lines, and labels.
fn draw_axes(
    svg: &mut String,
    area: PlotArea,
    max_count: f64,
    x_label: &str,
    histogram: &Histogram,
) {
    let left = area.left;
    let top = area.top;
    let plot_width = area.width;
    let plot_height = area.height;
    let axis_color = "#111827";
    let grid_color = "#e5e7eb";
    for tick in 0..=5 {
        let fraction = tick as f64 / 5.0;
        let y = top + plot_height - fraction * plot_height;
        let value = max_count * fraction;
        svg.push_str(&format!(
            r#"<line x1="{left:.2}" y1="{y:.2}" x2="{:.2}" y2="{y:.2}" stroke="{grid_color}" stroke-width="1"/>"#,
            left + plot_width
        ));
        svg.push_str(&format!(
            r#"<text x="{:.2}" y="{:.2}" text-anchor="end" font-family="sans-serif" font-size="12" fill="{axis_color}">{value:.0}</text>"#,
            left - 8.0,
            y + 4.0
        ));
    }
    svg.push_str(&format!(
        r#"<line x1="{left:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{axis_color}" stroke-width="1.5"/>"#,
        top + plot_height,
        left + plot_width,
        top + plot_height
    ));
    svg.push_str(&format!(
        r#"<line x1="{left:.2}" y1="{top:.2}" x2="{left:.2}" y2="{:.2}" stroke="{axis_color}" stroke-width="1.5"/>"#,
        top + plot_height
    ));
    for tick in 0..=4 {
        let fraction = tick as f64 / 4.0;
        let x = left + fraction * plot_width;
        let value = histogram.min + fraction * (histogram.max - histogram.min);
        svg.push_str(&format!(
            r#"<text x="{x:.2}" y="{:.2}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="{axis_color}">{value:.1}</text>"#,
            top + plot_height + 20.0
        ));
    }
    svg.push_str(&format!(
        r#"<text x="{:.2}" y="520" text-anchor="middle" font-family="sans-serif" font-size="15" fill="{axis_color}">{}</text>"#,
        left + plot_width / 2.0,
        escape_xml(x_label)
    ));
    svg.push_str(&format!(
        r#"<text x="22" y="{:.2}" text-anchor="middle" transform="rotate(-90 22 {:.2})" font-family="sans-serif" font-size="15" fill="{axis_color}">Count</text>"#,
        top + plot_height / 2.0,
        top + plot_height / 2.0
    ));
}

/// Plotting region dimensions inside an SVG.
#[derive(Debug, Clone, Copy)]
struct PlotArea {
    /// Left X coordinate.
    left: f64,
    /// Top Y coordinate.
    top: f64,
    /// Plotting width.
    width: f64,
    /// Plotting height.
    height: f64,
}

/// Writes per-MGF file statistics as CSV.
fn write_file_stats_csv(path: &Path, rows: &[MgfFileStats]) -> Result<()> {
    let mut file = std::fs::File::create(path).with_path(path)?;
    writeln!(
        file,
        "sample_id,polarity,path,features,total_peaks,mean_peaks_per_feature"
    )
    .with_path(path)?;
    for row in rows {
        writeln!(
            file,
            "{},{},{},{},{},{:.6}",
            csv_escape(&row.sample_id),
            row.polarity,
            csv_escape(&row.path.display().to_string()),
            row.features,
            row.total_peaks,
            row.mean_peaks_per_feature
        )
        .with_path(path)?;
    }
    Ok(())
}

/// Writes histogram bin counts as CSV.
fn write_histogram_csv(path: &Path, rows: &[HistogramRow]) -> Result<()> {
    let mut file = std::fs::File::create(path).with_path(path)?;
    writeln!(file, "metric,polarity,bin_start,bin_end,count").with_path(path)?;
    for row in rows {
        writeln!(
            file,
            "{},{},{:.8},{:.8},{}",
            row.metric, row.polarity, row.bin_start, row.bin_end, row.count
        )
        .with_path(path)?;
    }
    Ok(())
}

/// Writes the JSON plot report.
fn write_report(path: &Path, report: &PlotReport) -> Result<()> {
    let content = serde_json::to_vec_pretty(report).map_err(|source| Error::Json {
        path: path.to_path_buf(),
        source,
    })?;
    std::fs::write(path, content).with_path(path)
}

/// Escapes a field for simple CSV output.
fn csv_escape(value: &str) -> String {
    if value.contains([',', '"', '\n']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

/// Escapes text for XML output.
fn escape_xml(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
/// Tests for histogram and summary helpers.
mod tests {
    use super::*;

    #[test]
    /// Verifies percentile interpolation on sorted values.
    fn percentile_interpolates() {
        let values = [1.0, 2.0, 3.0, 4.0];

        assert_eq!(percentile(&values, 0.5), 2.5);
    }

    #[test]
    /// Verifies polarity detection from PF1600 MGF file names.
    fn detects_polarity_from_path() {
        assert_eq!(
            polarity_for_path(Path::new("VGF159_A02_features_ms2_neg.mgf")),
            PlotPolarity::Neg
        );
        assert_eq!(
            polarity_for_path(Path::new("VGF159_A02_features_ms2_pos.mgf")),
            PlotPolarity::Pos
        );
    }

    #[test]
    /// Verifies histogram binning clamps the maximum value into the last bin.
    fn bins_include_maximum_value() {
        let counts = bin_values(&[0.0, 0.5, 1.0], 0.0, 1.0, 2);

        assert_eq!(counts, vec![1, 2]);
    }
}
