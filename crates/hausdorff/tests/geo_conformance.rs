//! Conformance for the in-house planar/spherical algorithms, oracle-free.
//! Structural results are checked against independent dependency-free references
//! (shoelace area, min/max bounding box, ray-cast point-in-polygon, RDP
//! post-conditions); spherical/curve results without a simple closed form are
//! anchored by golden values captured from the validated implementation. No
//! third-party geometry crate.

use hausdorff::earth::{GeoPoint, frechet_distance_m, haversine_m};
use hausdorff::{
    Area, BoundingRect, Centroid, ChamberlainDuquetteArea, Contains, ConvexHull, CoordsIter,
    Simplify,
};
use hausdorff::{
    Bearing as _, Destination as _, Distance as _, InterpolatePoint as _, Length as _,
};

const EPS: f64 = 1e-9;

fn h_poly(coords: &[(f64, f64)]) -> hausdorff::Polygon {
    hausdorff::Polygon::new(
        hausdorff::LineString::new(
            coords
                .iter()
                .map(|&(x, y)| hausdorff::Coord::new(x, y))
                .collect(),
        ),
        vec![],
    )
}

fn fixtures() -> Vec<Vec<(f64, f64)>> {
    vec![
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        vec![(0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0)],
        vec![(0.0, 0.0), (6.0, 0.0), (0.0, 3.0)],
        vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 2.0),
            (2.0, 2.0),
            (2.0, 4.0),
            (0.0, 4.0),
        ],
        vec![
            (1.0, -3.1),
            (3.7, 2.7),
            (0.9, 7.6),
            (-4.8, 6.7),
            (-7.5, 0.9),
        ],
    ]
}

/// Independent signed area (shoelace) over the closed ring.
fn shoelace(coords: &[(f64, f64)]) -> f64 {
    let n = coords.len();
    let mut s = 0.0;
    for i in 0..n {
        let (x0, y0) = coords[i];
        let (x1, y1) = coords[(i + 1) % n];
        s += x0 * y1 - x1 * y0;
    }
    s * 0.5
}

/// Independent ray-casting point-in-polygon (interior; boundary is ambiguous).
fn ray_cast(coords: &[(f64, f64)], x: f64, y: f64) -> bool {
    let n = coords.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = coords[i];
        let (xj, yj) = coords[j];
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

#[path = "geo_conformance/basic.rs"]
mod basic;

#[path = "geo_conformance/simplify.rs"]
mod simplify;

#[path = "geo_conformance/geodesic.rs"]
mod geodesic;

#[path = "geo_conformance/fixtures.rs"]
mod fixtures;
