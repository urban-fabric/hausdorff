//! Cover an r2 [`Geometry`] with S2 cells at a fixed level — the geometry-enum
//! dispatch over [`RegionCoverer`]. The `Series`/WKB plumbing stays in the
//! caller; this is the pure geometry → cells core.

use std::collections::BTreeSet;

use crate::earth::GeoPoint;
use crate::space::r2::{Geometry, LineString, Polygon};

use super::coords::MAX_LEVEL;
use super::{Loop, Polygon as S2Polygon, RegionCoverer, S2CellId, S2Point};

/// S2 cells at `level` covering `geom`:
/// - **points** → the containing cell;
/// - **areal** geometries → the [`RegionCoverer`] cover (every level-`level` cell
///   that intersects the region — boundary-complete, not a centroid sample);
/// - **lines** → `None` (a region cover is undefined for a 1-D geometry);
/// - **collections** recurse.
///
/// Returns `None` for a line geometry (above), when `level` exceeds the leaf
/// level (30), or when a polygon ring has a non-finite or out-of-range
/// coordinate. An out-of-range *point* coordinate is skipped (a point is a single
/// sample, not a malformed ring); a degenerate hole (fewer than three vertices,
/// hence no area) is dropped.
#[must_use]
pub fn cover_r2_geometry_at_level(geom: &Geometry, level: u32) -> Option<Vec<S2CellId>> {
    if level > MAX_LEVEL {
        return None;
    }
    let mut cells = BTreeSet::new();
    collect(geom, level, &mut cells)?;
    Some(cells.into_iter().collect())
}

fn collect(geom: &Geometry, level: u32, out: &mut BTreeSet<S2CellId>) -> Option<()> {
    match geom {
        Geometry::Point(p) => insert_point(p.x(), p.y(), level, out),
        Geometry::MultiPoint(mp) => {
            for p in mp.iter() {
                insert_point(p.x(), p.y(), level, out);
            }
        }
        Geometry::Polygon(poly) => cover_polygon(poly, level, out)?,
        Geometry::MultiPolygon(mp) => {
            for poly in mp.iter() {
                cover_polygon(poly, level, out)?;
            }
        }
        Geometry::Rect(r) => cover_polygon(&r.to_polygon(), level, out)?,
        Geometry::Triangle(t) => cover_polygon(&t.to_polygon(), level, out)?,
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                collect(g, level, out)?;
            }
        }
        Geometry::Line(_) | Geometry::LineString(_) | Geometry::MultiLineString(_) => return None,
    }
    Some(())
}

/// True iff `(lon, lat)` is a finite, in-range WGS84 coordinate.
fn valid_lonlat(lon: f64, lat: f64) -> bool {
    lon.is_finite()
        && lat.is_finite()
        && (-90.0..=90.0).contains(&lat)
        && (-180.0..=180.0).contains(&lon)
}

fn insert_point(lon: f64, lat: f64, level: u32, out: &mut BTreeSet<S2CellId>) {
    if valid_lonlat(lon, lat) {
        let gp = GeoPoint::from_lon_lat(lon, lat);
        out.insert(S2CellId::from_geo_point_at_level(gp, level));
    }
}

fn cover_polygon(poly: &Polygon, level: u32, out: &mut BTreeSet<S2CellId>) -> Option<()> {
    // Reject malformed rings up front: an out-of-range lat/lon would otherwise be
    // silently wrapped by the trig in `S2Point::from_geo_point` rather than
    // flagged. (Unlike a stray point, a bad polygon ring is fatal.)
    for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
        if !ring.coords().all(|c| valid_lonlat(c.x, c.y)) {
            return None;
        }
    }
    let shell = ring_loop(poly.exterior())?;
    // Holes have valid coordinates here; a degenerate one (<3 vertices, no area)
    // simply drops.
    let holes: Vec<Loop> = poly.interiors().iter().filter_map(ring_loop).collect();
    out.extend(RegionCoverer::at_level(level).covering(&S2Polygon::new(shell, holes)));
    Some(())
}

/// A lon/lat ring as a counter-clockwise spherical [`Loop`] (interior to the left
/// — the small side); `None` for a ring with fewer than three vertices.
fn ring_loop(ring: &LineString) -> Option<Loop> {
    let pts: Vec<S2Point> = ring
        .coords()
        .map(|c| S2Point::from_geo_point(GeoPoint::from_lon_lat(c.x, c.y)))
        .collect();
    let l = Loop::new(pts);
    if l.num_vertices() < 3 {
        return None;
    }
    // A GIS ring never encloses more than a hemisphere; if this loop does it is
    // wound interior-to-the-right — reverse it so the interior is the small side.
    if l.area() > std::f64::consts::PI * 2.0 {
        let mut v = l.vertices().to_vec();
        v.reverse();
        Some(Loop::new(v))
    } else {
        Some(l)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::r2::{LineString, Point, Polygon};

    fn at(lon: f64, lat: f64, level: u32) -> S2CellId {
        S2CellId::from_geo_point_at_level(GeoPoint::from_lon_lat(lon, lat), level)
    }

    #[test]
    fn point_polygon_cover_and_line_skip() {
        let pt = Geometry::Point(Point::new(10.0, 50.0));
        assert_eq!(
            cover_r2_geometry_at_level(&pt, 12).unwrap(),
            vec![at(10.0, 50.0, 12)]
        );

        let level = 8;
        let poly = Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (9.5, 49.5),
                (10.5, 49.5),
                (10.5, 50.5),
                (9.5, 50.5),
                (9.5, 49.5),
            ]),
            vec![],
        ));
        let cover = cover_r2_geometry_at_level(&poly, level).unwrap();
        assert!(!cover.is_empty());
        assert!(cover.contains(&at(10.0, 50.0, level)));

        let line = Geometry::LineString(LineString::from(vec![(10.0, 50.0), (11.0, 50.0)]));
        assert!(cover_r2_geometry_at_level(&line, 12).is_none());
    }

    #[test]
    fn level_above_leaf_is_none() {
        let pt = Geometry::Point(Point::new(10.0, 50.0));
        assert!(cover_r2_geometry_at_level(&pt, 31).is_none());
    }

    #[test]
    fn polygon_with_out_of_range_coord_is_none() {
        let poly = Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (200.0, 0.0),
                (200.0, 1.0),
                (0.0, 1.0),
                (0.0, 0.0),
            ]),
            vec![],
        ));
        assert!(cover_r2_geometry_at_level(&poly, 8).is_none());
    }
}
