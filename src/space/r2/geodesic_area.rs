//! Spherical polygon area via the Chamberlain–Duquette formula — mirror of
//! `geo`'s `ChamberlainDuquetteArea`. Coordinates are interpreted as **lon/lat
//! degrees** (`x = lon`, `y = lat`); the result is in **metres²** on a sphere of
//! the equatorial Earth radius (matching geo, and the JPL Publication 07-03
//! derivation).
//!
//! This is the cheap spherical approximation. The ellipsoidal (Karney) geodesic
//! area is a future addition.

use super::{LineString, MultiPolygon, Polygon};

/// Equatorial Earth radius (WGS84), the radius the Chamberlain–Duquette formula
/// is conventionally evaluated at — kept identical to geo for exact agreement.
pub const EQUATORIAL_EARTH_RADIUS_M: f64 = 6_378_137.0;

/// Spherical (Chamberlain–Duquette) polygon area, signed and unsigned.
pub trait ChamberlainDuquetteArea {
    /// Signed area in m²: sign follows the exterior ring's winding (holes
    /// subtract).
    fn chamberlain_duquette_signed_area(&self) -> f64;

    /// Absolute area in m².
    fn chamberlain_duquette_unsigned_area(&self) -> f64 {
        self.chamberlain_duquette_signed_area().abs()
    }
}

/// Signed spherical area contribution of a single (closed) ring.
fn ring_area(coords: &[LineStringCoord]) -> f64 {
    let n = coords.len();
    if n <= 2 {
        return 0.0;
    }
    let mut total = 0.0;
    for i in 0..n {
        let (lower, middle, upper) = if i == n - 2 {
            (n - 2, n - 1, 0)
        } else if i == n - 1 {
            (n - 1, 0, 1)
        } else {
            (i, i + 1, i + 2)
        };
        let p1 = coords[lower];
        let p2 = coords[middle];
        let p3 = coords[upper];
        total += (p3.0.to_radians() - p1.0.to_radians()) * p2.1.to_radians().sin();
    }
    total * EQUATORIAL_EARTH_RADIUS_M * EQUATORIAL_EARTH_RADIUS_M / -2.0
}

type LineStringCoord = (f64, f64);

fn ring_lonlat(ring: &LineString) -> Vec<LineStringCoord> {
    ring.0.iter().map(|c| (c.x, c.y)).collect()
}

impl ChamberlainDuquetteArea for Polygon {
    fn chamberlain_duquette_signed_area(&self) -> f64 {
        let exterior = ring_area(&ring_lonlat(self.exterior()));
        self.interiors().iter().fold(exterior, |total, hole| {
            total - ring_area(&ring_lonlat(hole))
        })
    }
}

impl ChamberlainDuquetteArea for MultiPolygon {
    fn chamberlain_duquette_signed_area(&self) -> f64 {
        self.iter()
            .map(ChamberlainDuquetteArea::chamberlain_duquette_signed_area)
            .sum()
    }
}
