//! WGS84 coordinate and its compact fixed-point storage encoding.

#![allow(unsafe_code)]

use bytemuck::{Pod, Zeroable};
use serde::{Deserialize, Serialize};

/// A WGS84 coordinate in decimal degrees.
#[derive(Clone, Copy, PartialEq, Debug, Pod, Zeroable, Serialize, Deserialize)]
#[repr(C)]
pub struct GeoPoint {
    pub lat: f64,
    pub lon: f64,
}

impl GeoPoint {
    pub fn from_lon_lat(lon: f64, lat: f64) -> Self {
        Self { lat, lon }
    }

    pub fn from_lat_lon(lat: f64, lon: f64) -> Self {
        Self { lat, lon }
    }

    #[inline]
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.lat.is_finite() && self.lon.is_finite()
    }
}

/// Storage encoding of a [`GeoPoint`]: fixed-point degrees scaled by 1e7 — the
/// OSM-native "100 nanodegrees" unit. 8 bytes/point versus 16 for `GeoPoint`,
/// with ~1.1 cm resolution, and `|deg| * 1e7` always fits in `i32`
/// (`180 * 1e7 < 2^31`).
#[derive(Clone, Copy, PartialEq, Eq, Debug, Pod, Zeroable)]
#[repr(C)]
pub struct FixedGeoPoint {
    pub lat_e7: i32,
    pub lon_e7: i32,
}

/// Degrees → fixed-point scale (1e7 ⇒ ~1.1 cm; `±180 * 1e7` fits in `i32`).
const GEO_FIXED_SCALE: f64 = 1e7;

impl FixedGeoPoint {
    #[inline]
    #[must_use]
    pub fn from_geo(p: GeoPoint) -> Self {
        Self {
            lat_e7: (p.lat * GEO_FIXED_SCALE).round() as i32,
            lon_e7: (p.lon * GEO_FIXED_SCALE).round() as i32,
        }
    }

    #[inline]
    #[must_use]
    pub fn to_geo(self) -> GeoPoint {
        GeoPoint {
            lat: f64::from(self.lat_e7) / GEO_FIXED_SCALE,
            lon: f64::from(self.lon_e7) / GEO_FIXED_SCALE,
        }
    }
}

impl From<GeoPoint> for FixedGeoPoint {
    fn from(p: GeoPoint) -> Self {
        Self::from_geo(p)
    }
}

impl From<FixedGeoPoint> for GeoPoint {
    fn from(p: FixedGeoPoint) -> Self {
        p.to_geo()
    }
}
