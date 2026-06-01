//! Interchange between the **spherical** WGS84 geometry of [`crate::geo`] and the
//! **planar** Cartesian geometry of [`crate::planar`].
//!
//! The mapping is the GIS convention `x = lon`, `y = lat`: a spherical
//! [`GeoPoint`] becomes a planar [`Coord`]/[`Point`] carrying longitude in `x`
//! and latitude in `y`, and back. This is a *coordinate relabelling*, not a map
//! projection — distances and areas in the planar plane are in degrees, not
//! metres. Project first (or use the spherical algorithms) when you need metric
//! results.

use crate::geo::GeoPoint;
use crate::planar::{Coord, LineString, Point, Polygon};

/// A spherical [`GeoPoint`] as a planar [`Coord`] with `x = lon`, `y = lat`.
#[must_use]
pub fn geo_point_to_coord(p: GeoPoint) -> Coord {
    Coord::new(p.lon, p.lat)
}

/// A planar [`Coord`] (`x = lon`, `y = lat`) as a spherical [`GeoPoint`].
#[must_use]
pub fn coord_to_geo_point(c: Coord) -> GeoPoint {
    GeoPoint::from_lat_lon(c.y, c.x)
}

/// A spherical [`GeoPoint`] as a planar [`Point`] (`x = lon`, `y = lat`).
#[must_use]
pub fn geo_point_to_point(p: GeoPoint) -> Point {
    Point(geo_point_to_coord(p))
}

/// A planar [`Point`] (`x = lon`, `y = lat`) as a spherical [`GeoPoint`].
#[must_use]
pub fn point_to_geo_point(p: Point) -> GeoPoint {
    coord_to_geo_point(p.0)
}

/// A run of spherical points as a planar [`LineString`] (`x = lon`, `y = lat`).
#[must_use]
pub fn geo_points_to_line_string(points: &[GeoPoint]) -> LineString {
    LineString::new(points.iter().copied().map(geo_point_to_coord).collect())
}

/// A planar [`LineString`] (`x = lon`, `y = lat`) as a `Vec` of spherical points.
#[must_use]
pub fn line_string_to_geo_points(ls: &LineString) -> Vec<GeoPoint> {
    ls.coords_copied().map(coord_to_geo_point).collect()
}

/// The exterior ring of a planar [`Polygon`] (`x = lon`, `y = lat`) as spherical
/// points. Holes are ignored — the spherical side is ring-oriented.
#[must_use]
pub fn polygon_exterior_to_geo_points(poly: &Polygon) -> Vec<GeoPoint> {
    line_string_to_geo_points(poly.exterior())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geo_point_round_trips_through_coord() {
        let p = GeoPoint::from_lat_lon(40.0, -74.0);
        let c = geo_point_to_coord(p);
        assert_eq!(c, Coord::new(-74.0, 40.0)); // x = lon, y = lat
        let back = coord_to_geo_point(c);
        assert_eq!(back, p);
    }

    #[test]
    fn line_string_round_trips() {
        let pts = [
            GeoPoint::from_lat_lon(0.0, 0.0),
            GeoPoint::from_lat_lon(1.0, 2.0),
        ];
        let ls = geo_points_to_line_string(&pts);
        assert_eq!(ls.0[1], Coord::new(2.0, 1.0));
        assert_eq!(line_string_to_geo_points(&ls), pts.to_vec());
    }
}
