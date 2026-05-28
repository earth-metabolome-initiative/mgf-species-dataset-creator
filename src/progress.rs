//! Terminal progress reporting helpers used by the CLI pipeline.

use std::path::Path;

use indicatif::{ProgressBar, ProgressStyle};

/// Creates consistently styled spinners and bars for pipeline progress.
#[derive(Debug, Clone)]
pub struct ProgressReporter {
    enabled: bool,
}

impl ProgressReporter {
    /// Creates an enabled progress reporter.
    #[must_use]
    pub const fn enabled() -> Self {
        Self { enabled: true }
    }

    /// Starts a spinner for a named stage.
    #[must_use]
    pub fn spinner(&self, message: impl Into<String>) -> ProgressBar {
        if !self.enabled {
            return ProgressBar::hidden();
        }
        let spinner = ProgressBar::new_spinner();
        spinner.set_style(
            ProgressStyle::with_template("{spinner:.green} {msg}")
                .unwrap_or_else(|_| ProgressStyle::default_spinner()),
        );
        spinner.enable_steady_tick(std::time::Duration::from_millis(120));
        spinner.set_message(message.into());
        spinner
    }

    /// Starts a determinate progress bar.
    #[must_use]
    pub fn bar(&self, len: u64, message: impl Into<String>) -> ProgressBar {
        if !self.enabled {
            return ProgressBar::hidden();
        }
        let bar = ProgressBar::new(len);
        bar.set_style(
            ProgressStyle::with_template(
                "{msg}\n[{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {wide_msg}",
            )
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars("=> "),
        );
        bar.set_message(message.into());
        bar
    }

    /// Starts a byte-count progress bar for downloads.
    #[must_use]
    pub fn bytes_bar(&self, len: u64, message: impl Into<String>) -> ProgressBar {
        if !self.enabled {
            return ProgressBar::hidden();
        }
        let bar = ProgressBar::new(len);
        bar.set_style(
            ProgressStyle::with_template(
                "{msg}\n[{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({bytes_per_sec}, ETA {eta})",
            )
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars("=> "),
        );
        bar.set_message(message.into());
        bar
    }

    /// Formats a file name for compact progress messages.
    #[must_use]
    pub fn file_label(path: &Path) -> String {
        path.file_name()
            .and_then(|name| name.to_str())
            .map(str::to_string)
            .unwrap_or_else(|| path.display().to_string())
    }
}
