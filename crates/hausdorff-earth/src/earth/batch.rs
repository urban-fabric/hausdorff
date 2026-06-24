//! `earth::batch` — Structure-of-Arrays batches of WGS84 lon/lat points for
//! allocation-friendly bulk operations.
//!
//! The core point model stays Array-of-Structs ([`GeoPoint`]); a [`LonLatBatch`]
//! is an *optional* columnar view (`lons[]`, `lats[]`) for cache-friendly bulk
//! kernels: distance fans (nearest-point queries), whole-batch H3 encoding, and
//! bounding boxes. The kernels are sequential; a caller wanting thread
//! parallelism splits the columns into chunks and processes those concurrently
//! on its own pool.

use crate::earth::{BBox, GeoPoint, LengthMeters};
use hausdorff_index::index::h3::H3Cell;
use hausdorff_kernels::earth as earth_kernels;

/// A batch of geographic points stored column-wise (`lons`, `lats`, equal length)
/// for bulk kernels. Convert to/from `&[GeoPoint]` at the edges.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LonLatBatch {
    lons: Vec<f64>,
    lats: Vec<f64>,
}

impl LonLatBatch {
    /// Build from parallel `lons`/`lats` columns, or `None` if their lengths
    /// differ.
    #[must_use]
    pub fn from_columns(lons: Vec<f64>, lats: Vec<f64>) -> Option<Self> {
        (lons.len() == lats.len()).then_some(Self { lons, lats })
    }

    /// Build the columnar batch from an Array-of-Structs slice of points.
    #[must_use]
    pub fn from_points(points: &[GeoPoint]) -> Self {
        Self {
            lons: points.iter().map(|p| p.lon).collect(),
            lats: points.iter().map(|p| p.lat).collect(),
        }
    }

    /// Number of points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.lons.len()
    }

    /// Whether the batch has no points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.lons.is_empty()
    }

    /// The longitude column (degrees).
    #[must_use]
    pub fn lons(&self) -> &[f64] {
        &self.lons
    }

    /// The latitude column (degrees).
    #[must_use]
    pub fn lats(&self) -> &[f64] {
        &self.lats
    }

    /// The `i`-th point as a [`GeoPoint`].
    #[must_use]
    pub fn point(&self, i: usize) -> GeoPoint {
        GeoPoint::from_lat_lon(self.lats[i], self.lons[i])
    }

    /// Great-circle distance from every point to a fixed `target` — the
    /// one-to-many "distance fan" behind nearest-point queries.
    #[must_use]
    pub fn haversine_to(&self, target: GeoPoint) -> Vec<LengthMeters> {
        self.haversine_to_m(target)
            .into_iter()
            .map(LengthMeters::new)
            .collect()
    }

    /// Great-circle distance in raw metres from every point to a fixed `target`.
    #[must_use]
    pub fn haversine_to_m(&self, target: GeoPoint) -> Vec<f64> {
        self.lats
            .iter()
            .zip(&self.lons)
            .map(|(&lat, &lon)| earth_kernels::haversine_m(lat, lon, target.lat, target.lon))
            .collect()
    }

    /// Mask for rows that are finite, in-range WGS84 latitude/longitude points.
    #[must_use]
    pub fn valid_lonlat_mask(&self) -> Vec<bool> {
        self.lats
            .iter()
            .zip(&self.lons)
            .map(|(&lat, &lon)| earth_kernels::valid_latlon_degrees(lat, lon))
            .collect()
    }

    /// Mask for rows inside `bbox`, using inclusive latitude/longitude bounds.
    #[must_use]
    pub fn bbox_contains_mask(&self, bbox: BBox) -> Vec<bool> {
        self.lats
            .iter()
            .zip(&self.lons)
            .map(|(&lat, &lon)| {
                earth_kernels::bbox_contains_latlon(
                    bbox.min_lat,
                    bbox.min_lon,
                    bbox.max_lat,
                    bbox.max_lon,
                    lat,
                    lon,
                )
            })
            .collect()
    }

    /// Mask for rows within `radius_m` metres of `target`.
    #[must_use]
    pub fn within_radius_m(&self, target: GeoPoint, radius_m: f64) -> Vec<bool> {
        self.lats
            .iter()
            .zip(&self.lons)
            .map(|(&lat, &lon)| {
                earth_kernels::valid_haversine_within_radius_m(
                    lat, lon, target.lat, target.lon, radius_m,
                )
            })
            .collect()
    }

    /// Mask for rows within `radius` of `target`.
    #[must_use]
    pub fn within_radius(&self, target: GeoPoint, radius: LengthMeters) -> Vec<bool> {
        self.within_radius_m(target, radius.value())
    }

    /// Mask for rows within `radius_m` metres of a target latitude/longitude
    /// pair without constructing a [`GeoPoint`].
    #[must_use]
    pub fn within_radius_latlon_m(
        &self,
        target_lat: f64,
        target_lon: f64,
        radius_m: f64,
    ) -> Vec<bool> {
        self.lats
            .iter()
            .zip(&self.lons)
            .map(|(&lat, &lon)| {
                earth_kernels::valid_haversine_within_radius_m(
                    lat, lon, target_lat, target_lon, radius_m,
                )
            })
            .collect()
    }

    /// The H3 cell containing each point at `resolution` (per-point `None` for an
    /// out-of-range coordinate), encoding the whole batch at once.
    #[must_use]
    pub fn h3_cells(&self, resolution: u8) -> Vec<Option<H3Cell>> {
        (0..self.len())
            .map(|i| H3Cell::for_point(self.lats[i], self.lons[i], resolution))
            .collect()
    }

    /// The axis-aligned lon/lat bounding box of the batch, or `None` if empty.
    /// A single min/max reduction over the lon/lat columns.
    #[must_use]
    pub fn bounding_box(&self) -> Option<BBox> {
        if self.is_empty() {
            return None;
        }
        let (mut min_lat, mut max_lat) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut min_lon, mut max_lon) = (f64::INFINITY, f64::NEG_INFINITY);
        for (&lat, &lon) in self.lats.iter().zip(&self.lons) {
            min_lat = min_lat.min(lat);
            max_lat = max_lat.max(lat);
            min_lon = min_lon.min(lon);
            max_lon = max_lon.max(lon);
        }
        Some(BBox::new(min_lat, min_lon, max_lat, max_lon))
    }
}
