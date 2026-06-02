//! `hausdorff.toml` — empirical, machine-tunable parameters.
//!
//! Some kernels carry constants that trade off costs which depend on the target
//! CPU (cache sizes, branch costs, allocation speed), so the right value differs
//! per machine — for example the [`CostModel`] weights that decide when the
//! minimum-distance kernel switches from a brute scan to a spatial index, or the
//! STR-tree [`Index::node_capacity`]. Rather than hard-coding them, they live in
//! the bundled `hausdorff.toml`, loaded via [`cavalieri`] and overridable per
//! host by pointing the [`CONFIG_ENV`] environment variable at an override file.
//!
//! Read the loaded values through [`config()`] (parsed once, lazily).

use std::path::Path;
use std::sync::LazyLock;

use serde::Deserialize;

/// The cavalieri loader contract, re-exported so callers can parse overrides
/// (`Config::from_path`, `Config::from_toml_str`) without depending on cavalieri
/// directly.
pub use cavalieri::{BundledConfig, ConfigLoadError};

/// Environment variable naming an override `hausdorff.toml` (absolute or relative
/// path). Fields it omits fall back to the bundled defaults.
pub const CONFIG_ENV: &str = "HAUSDORFF_CONFIG";

/// All of hausdorff's tunable parameters — the parsed `hausdorff.toml`.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Config {
    /// Cost-model weights for path-selecting (query-planner-style) kernels.
    pub cost_model: CostModel,
    /// Spatial-index defaults.
    pub index: Index,
}

/// Cost-model weights for path-selecting kernels.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CostModel {
    /// The minimum-distance brute-vs-index dispatch.
    pub min_distance: MinDistanceCost,
}

/// Relative cost weights for the minimum-distance brute-vs-index dispatch. The
/// index path runs when `index_build_weight·m + index_query_weight·n·log₂m` is
/// below the brute estimate `n·m` (and the operands' bounding boxes overlap).
#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct MinDistanceCost {
    /// Per-segment weight of building the spatial index.
    pub index_build_weight: f64,
    /// Per-query weight (× `log₂` of the indexed set) of querying it.
    pub index_query_weight: f64,
    /// Segment-pair count (`n·m`) at or below which the brute scan skips its
    /// per-pair bounding-box pruning: for small operands the prune bookkeeping
    /// costs more than the few distance evaluations it saves, so a tight
    /// unpruned `O(n·m)` scan is faster.
    pub prune_pair_threshold: usize,
}

impl Default for MinDistanceCost {
    fn default() -> Self {
        Self { index_build_weight: 6.0, index_query_weight: 6.0, prune_pair_threshold: 256 }
    }
}

/// Spatial-index defaults.
#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Index {
    /// Default items per STR-tree node (fan-out).
    pub node_capacity: usize,
}

impl Default for Index {
    fn default() -> Self {
        Self { node_capacity: 16 }
    }
}

impl BundledConfig for Config {
    const BUNDLED_TOML: &'static str = include_str!("../hausdorff.toml");
    const NAME: &'static str = "hausdorff";
}

impl Config {
    /// An override file named by [`CONFIG_ENV`] if it parses, else the bundled
    /// `hausdorff.toml` defaults (a checked-in invariant, so it must parse).
    fn load() -> Self {
        if let Some(path) = std::env::var_os(CONFIG_ENV)
            && let Ok(cfg) = Config::from_path(Path::new(&path))
        {
            return cfg;
        }
        Config::bundled_or_panic()
    }
}

static CONFIG: LazyLock<Config> = LazyLock::new(Config::load);

/// The process-wide configuration, parsed once on first access from
/// `hausdorff.toml` (or the [`CONFIG_ENV`] override).
#[must_use]
pub fn config() -> &'static Config {
    &CONFIG
}
