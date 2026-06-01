//! Axis-aligned geographic bounding box.

use serde::{Deserialize, Serialize};

use super::GeoPoint;

/// An axis-aligned geographic bounding box in WGS84 degrees.
#[derive(Clone, Copy, PartialEq, Debug, Serialize, Deserialize)]
pub struct BBox {
    pub min_lat: f64,
    pub min_lon: f64,
    pub max_lat: f64,
    pub max_lon: f64,
}

impl BBox {
    /// The degenerate box at the origin.
    pub const ZERO: Self = Self {
        min_lat: 0.0,
        min_lon: 0.0,
        max_lat: 0.0,
        max_lon: 0.0,
    };

    /// Construct a bounding box from its latitude/longitude extents.
    pub fn new(min_lat: f64, min_lon: f64, max_lat: f64, max_lon: f64) -> Self {
        Self {
            min_lat,
            min_lon,
            max_lat,
            max_lon,
        }
    }

    /// True iff `point` lies within the box (inclusive); non-finite points are
    /// never contained.
    pub fn contains(&self, point: GeoPoint) -> bool {
        if !point.lat.is_finite() || !point.lon.is_finite() {
            return false;
        }
        point.lat >= self.min_lat
            && point.lat <= self.max_lat
            && point.lon >= self.min_lon
            && point.lon <= self.max_lon
    }

    /// The tight box covering all finite points in `coords`
    /// ([`ZERO`](Self::ZERO) if none are finite).
    pub fn from_coords(coords: &[GeoPoint]) -> Self {
        let mut min_lat = f64::INFINITY;
        let mut min_lon = f64::INFINITY;
        let mut max_lat = f64::NEG_INFINITY;
        let mut max_lon = f64::NEG_INFINITY;
        for c in coords {
            if !c.lat.is_finite() || !c.lon.is_finite() {
                continue;
            }
            if c.lat < min_lat {
                min_lat = c.lat;
            }
            if c.lat > max_lat {
                max_lat = c.lat;
            }
            if c.lon < min_lon {
                min_lon = c.lon;
            }
            if c.lon > max_lon {
                max_lon = c.lon;
            }
        }
        if !min_lat.is_finite() {
            Self::ZERO
        } else {
            Self {
                min_lat,
                min_lon,
                max_lat,
                max_lon,
            }
        }
    }
}

/// Error parsing a [`BBox`] from `"min_lat,min_lon,max_lat,max_lon"`.
#[derive(Debug, Clone, PartialEq)]
pub enum BBoxParseError {
    /// Wrong number of comma-separated components (expected 4).
    WrongArity(usize),
    /// A component failed to parse as `f64`.
    NumberParse {
        index: usize,
        raw: String,
        source: std::num::ParseFloatError,
    },
    /// A component was not finite.
    NotFinite { index: usize, value: f64 },
    /// `min` exceeded `max` on an axis.
    Inverted {
        min_lat: f64,
        max_lat: f64,
        min_lon: f64,
        max_lon: f64,
    },
}

impl std::fmt::Display for BBoxParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::WrongArity(n) => write!(
                f,
                "bbox needs 4 comma-separated values (min_lat,min_lon,max_lat,max_lon); got {n}"
            ),
            Self::NumberParse { index, raw, source } => {
                write!(f, "bbox component {index} ({raw:?}): {source}")
            }
            Self::NotFinite { index, value } => {
                write!(f, "bbox component {index} must be finite (got {value})")
            }
            Self::Inverted {
                min_lat,
                max_lat,
                min_lon,
                max_lon,
            } => write!(
                f,
                "bbox min must be ≤ max in both axes: lat=[{min_lat},{max_lat}] lon=[{min_lon},{max_lon}]"
            ),
        }
    }
}

impl std::error::Error for BBoxParseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::NumberParse { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl std::str::FromStr for BBox {
    type Err = BBoxParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 4 {
            return Err(BBoxParseError::WrongArity(parts.len()));
        }
        let parse = |index: usize| -> Result<f64, BBoxParseError> {
            let raw = parts[index].trim();
            let value = raw
                .parse::<f64>()
                .map_err(|source| BBoxParseError::NumberParse {
                    index,
                    raw: raw.to_string(),
                    source,
                })?;
            if !value.is_finite() {
                return Err(BBoxParseError::NotFinite { index, value });
            }
            Ok(value)
        };
        let min_lat = parse(0)?;
        let min_lon = parse(1)?;
        let max_lat = parse(2)?;
        let max_lon = parse(3)?;
        if min_lat > max_lat || min_lon > max_lon {
            return Err(BBoxParseError::Inverted {
                min_lat,
                max_lat,
                min_lon,
                max_lon,
            });
        }
        Ok(BBox {
            min_lat,
            min_lon,
            max_lat,
            max_lon,
        })
    }
}
