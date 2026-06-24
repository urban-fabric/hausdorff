//! Empirical, machine-tunable parameters.
//!
//! Some kernels carry constants that trade off costs which depend on the target
//! CPU (cache sizes, branch costs, allocation speed), so the right value differs
//! per machine — for example the [`CostModel`] weights that decide when the
//! minimum-distance kernel switches from a brute scan to a spatial index, or the
//! STR-tree [`Index::node_capacity`]. They are bundled here as in-code defaults
//! and overridable per host by pointing the [`CONFIG_ENV`] environment variable
//! at an override TOML file.
//!
//! Read the loaded values through [`config()`] (parsed once, lazily).

use std::error::Error;
use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use serde::Deserialize;

/// Error while reading or parsing a bundled config.
#[derive(Debug)]
pub enum ConfigLoadError {
    /// Failed to read an override file.
    Io {
        /// Logical config name.
        name: &'static str,
        /// Path that was read.
        path: PathBuf,
        /// Underlying I/O error.
        source: std::io::Error,
    },
    /// Failed to parse an override file.
    ParseAtPath {
        /// Logical config name.
        name: &'static str,
        /// Path that was parsed.
        path: PathBuf,
        /// Underlying parse error.
        source: ConfigParseError,
    },
    /// Failed to parse a config string.
    Parse {
        /// Logical config name.
        name: &'static str,
        /// Underlying parse error.
        source: ConfigParseError,
    },
}

impl fmt::Display for ConfigLoadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { name, path, source } => {
                write!(f, "read {name} config {}: {source}", path.display())
            }
            Self::ParseAtPath { name, path, source } => {
                write!(f, "parse {name} config {}: {source}", path.display())
            }
            Self::Parse { name, source } => write!(f, "parse {name} config: {source}"),
        }
    }
}

impl Error for ConfigLoadError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::ParseAtPath { source, .. } | Self::Parse { source, .. } => Some(source),
        }
    }
}

/// Error while parsing config TOML.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConfigParseError {
    line: usize,
    message: String,
}

impl ConfigParseError {
    fn new(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for ConfigParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl Error for ConfigParseError {}

/// Loader contract for config types that embed a default TOML body.
pub trait BundledConfig: Sized {
    /// Bundled TOML contents.
    const BUNDLED_TOML: &'static str;

    /// Logical config name used in errors.
    const NAME: &'static str;

    #[doc(hidden)]
    fn parse_toml_str(body: &str) -> Result<Self, ConfigParseError>;

    /// Parse the bundled TOML defaults.
    fn bundled() -> Result<Self, ConfigLoadError> {
        Self::from_toml_str(Self::BUNDLED_TOML)
    }

    /// Parse an override TOML string. Missing fields fall back to defaults.
    fn from_toml_str(body: &str) -> Result<Self, ConfigLoadError> {
        Self::parse_toml_str(body).map_err(|source| ConfigLoadError::Parse {
            name: Self::NAME,
            source,
        })
    }

    /// Parse an override TOML file. Missing fields fall back to defaults.
    fn from_path(path: &Path) -> Result<Self, ConfigLoadError> {
        let body = std::fs::read_to_string(path).map_err(|source| ConfigLoadError::Io {
            name: Self::NAME,
            path: path.to_path_buf(),
            source,
        })?;
        Self::parse_toml_str(&body).map_err(|source| ConfigLoadError::ParseAtPath {
            name: Self::NAME,
            path: path.to_path_buf(),
            source,
        })
    }

    #[allow(
        clippy::panic,
        reason = "bundled config parse failure means a checked-in TOML literal is invalid"
    )]
    /// Parse the bundled TOML defaults, panicking if the checked-in literal is invalid.
    fn bundled_or_panic() -> Self {
        Self::bundled()
            .unwrap_or_else(|err| panic!("bundled {} config must parse: {err}", Self::NAME))
    }
}

/// Environment variable naming an override TOML file (absolute or relative path).
/// Fields it omits fall back to the bundled defaults.
pub const CONFIG_ENV: &str = "HAUSDORFF_CONFIG";

/// All of hausdorff's tunable parameters.
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
        Self {
            index_build_weight: 6.0,
            index_query_weight: 6.0,
            prune_pair_threshold: 256,
        }
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
    const BUNDLED_TOML: &'static str = r#"
[cost_model.min_distance]
index_build_weight = 6.0
index_query_weight = 6.0
prune_pair_threshold = 256

[index]
node_capacity = 16
"#;
    const NAME: &'static str = "hausdorff";

    fn parse_toml_str(body: &str) -> Result<Self, ConfigParseError> {
        parse_config_toml(body)
    }
}

#[derive(Clone, Copy)]
enum Section {
    Root,
    MinDistance,
    Index,
}

#[derive(Default)]
struct Seen {
    index_build_weight: bool,
    index_query_weight: bool,
    prune_pair_threshold: bool,
    node_capacity: bool,
}

fn parse_config_toml(body: &str) -> Result<Config, ConfigParseError> {
    let mut cfg = Config::default();
    let mut section = Section::Root;
    let mut seen = Seen::default();

    for (line_idx, raw) in body.lines().enumerate() {
        let line_no = line_idx + 1;
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with('[') {
            if !line.ends_with(']') {
                return Err(ConfigParseError::new(line_no, "malformed table header"));
            }
            section = match line[1..line.len() - 1].trim() {
                "cost_model.min_distance" => Section::MinDistance,
                "index" => Section::Index,
                "" => {
                    return Err(ConfigParseError::new(line_no, "empty table header"));
                }
                table => {
                    return Err(ConfigParseError::new(
                        line_no,
                        format!("unknown table `{table}`"),
                    ));
                }
            };
            continue;
        }

        let (key, value) = line
            .split_once('=')
            .ok_or_else(|| ConfigParseError::new(line_no, "expected `key = value`"))?;
        let key = key.trim();
        let value = value.trim();
        if key.is_empty() {
            return Err(ConfigParseError::new(line_no, "empty key"));
        }
        if value.is_empty() {
            return Err(ConfigParseError::new(
                line_no,
                format!("empty value for `{key}`"),
            ));
        }

        match section {
            Section::Root => {
                return Err(ConfigParseError::new(
                    line_no,
                    format!("unknown top-level key `{key}`"),
                ));
            }
            Section::MinDistance => match key {
                "index_build_weight" => {
                    mark_unseen(&mut seen.index_build_weight, line_no, key)?;
                    cfg.cost_model.min_distance.index_build_weight =
                        parse_f64(value, line_no, key)?;
                }
                "index_query_weight" => {
                    mark_unseen(&mut seen.index_query_weight, line_no, key)?;
                    cfg.cost_model.min_distance.index_query_weight =
                        parse_f64(value, line_no, key)?;
                }
                "prune_pair_threshold" => {
                    mark_unseen(&mut seen.prune_pair_threshold, line_no, key)?;
                    cfg.cost_model.min_distance.prune_pair_threshold =
                        parse_usize(value, line_no, key)?;
                }
                _ => {
                    return Err(ConfigParseError::new(
                        line_no,
                        format!("unknown key `{key}` in [cost_model.min_distance]"),
                    ));
                }
            },
            Section::Index => match key {
                "node_capacity" => {
                    mark_unseen(&mut seen.node_capacity, line_no, key)?;
                    cfg.index.node_capacity = parse_usize(value, line_no, key)?;
                }
                _ => {
                    return Err(ConfigParseError::new(
                        line_no,
                        format!("unknown key `{key}` in [index]"),
                    ));
                }
            },
        }
    }

    Ok(cfg)
}

fn strip_comment(line: &str) -> &str {
    line.split_once('#').map_or(line, |(head, _)| head)
}

fn mark_unseen(seen: &mut bool, line_no: usize, key: &str) -> Result<(), ConfigParseError> {
    if *seen {
        return Err(ConfigParseError::new(
            line_no,
            format!("duplicate key `{key}`"),
        ));
    }
    *seen = true;
    Ok(())
}

fn parse_f64(value: &str, line_no: usize, key: &str) -> Result<f64, ConfigParseError> {
    value.parse::<f64>().map_err(|_| {
        ConfigParseError::new(line_no, format!("invalid float for `{key}`: `{value}`"))
    })
}

fn parse_usize(value: &str, line_no: usize, key: &str) -> Result<usize, ConfigParseError> {
    value.parse::<usize>().map_err(|_| {
        ConfigParseError::new(line_no, format!("invalid integer for `{key}`: `{value}`"))
    })
}

impl Config {
    /// An override file named by [`CONFIG_ENV`] if it parses, else the bundled
    /// defaults (a checked-in invariant, so they must parse).
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

/// The process-wide configuration, parsed once on first access from the bundled
/// defaults (or the [`CONFIG_ENV`] override).
#[must_use]
pub fn config() -> &'static Config {
    &CONFIG
}
