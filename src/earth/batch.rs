//! `earth::batch` — Structure-of-Arrays batches of WGS84 lon/lat points for
//! allocation-friendly, optionally-parallel bulk operations.
//!
//! The core point model stays Array-of-Structs ([`GeoPoint`]); a [`LonLatBatch`]
//! is an *optional* columnar view (`lons[]`, `lats[]`) for cache-friendly,
//! data-parallel bulk kernels: distance fans (nearest-point queries), whole-batch
//! H3 encoding, and bounding boxes. The per-element kernels run on rayon when the
//! `parallel` feature is enabled (otherwise sequentially); results are identical.

use euclid::LengthMeters;

use crate::earth::{BBox, GeoPoint, haversine_m};
use crate::index::h3::H3Cell;
use crate::par::map_indexed;

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
        map_indexed(self.len(), |i| {
            LengthMeters::new(haversine_m(self.point(i), target))
        })
    }

    /// The H3 cell containing each point at `resolution` (per-point `None` for an
    /// out-of-range coordinate), encoding the whole batch at once.
    #[must_use]
    pub fn h3_cells(&self, resolution: u8) -> Vec<Option<H3Cell>> {
        map_indexed(self.len(), |i| {
            H3Cell::for_point(self.lats[i], self.lons[i], resolution)
        })
    }

    /// The axis-aligned lon/lat bounding box of the batch, or `None` if empty.
    /// A sequential min/max reduction (cheap and memory-bound, so not parallel).
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
