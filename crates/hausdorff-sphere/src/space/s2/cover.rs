//! Cover an r2 [`Geometry`] with S2 cells — the geometry-enum dispatch over
//! [`RegionCoverer`]. The `Series`/WKB plumbing stays in the caller; this is
//! the pure geometry → cells core. [`cover_r2_geometry_at_level`] emits a
//! fixed-level cover; [`cover_r2_geometry`] emits an adaptive multi-level
//! [`CellUnion`].

use std::collections::BTreeSet;

use hausdorff_base::earth::GeoPoint;
use hausdorff_planar::space::r2::{Geometry, LineString, Polygon};

use super::coords::MAX_LEVEL;
use super::{CellUnion, Loop, Polygon as S2Polygon, Polyline, RegionCoverer, S2CellId, S2Point};

/// S2 cells at `level` covering `geom`:
/// - **points** → the containing cell;
/// - **areal** geometries → the [`RegionCoverer`] cover (every level-`level` cell
///   that intersects the region — boundary-complete, not a centroid sample);
/// - **lines** → the [`Polyline`] boundary cover (every level-`level` cell the
///   chain of great-circle edges passes through);
/// - **collections** recurse.
///
/// Returns `None` when `level` exceeds the leaf level (30), or when a polygon
/// ring or line has a non-finite or out-of-range coordinate. An out-of-range
/// *point* coordinate is skipped (a point is a single sample, not a malformed
/// ring); a degenerate hole (fewer than three vertices, hence no area) is
/// dropped.
#[must_use]
pub fn cover_r2_geometry_at_level(geom: &Geometry, level: u32) -> Option<Vec<S2CellId>> {
    if level > MAX_LEVEL {
        return None;
    }
    let mut cells = BTreeSet::new();
    collect(geom, level, &mut cells)?;
    Some(cells.into_iter().collect())
}

/// An adaptive multi-level cover of `geom` driven by `coverer`
/// (see [`RegionCoverer::adaptive_covering`]): points map to their cell at
/// `coverer.max_level`, areal geometries to the adaptive region cover, lines to
/// the adaptive [`Polyline`] boundary cover, and collections to the normalized
/// union of their parts. Coordinate validation matches
/// [`cover_r2_geometry_at_level`].
#[must_use]
pub fn cover_r2_geometry(geom: &Geometry, coverer: &RegionCoverer) -> Option<CellUnion> {
    let mut cells = Vec::new();
    collect_adaptive(geom, coverer, &mut cells)?;
    Some(CellUnion::from(cells))
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
        Geometry::Line(l) => {
            let coords = [(l.start.x, l.start.y), (l.end.x, l.end.y)];
            out.extend(cover_line_coords(&coords, &RegionCoverer::at_level(level))?);
        }
        Geometry::LineString(ls) => {
            out.extend(cover_line_string(ls, &RegionCoverer::at_level(level))?);
        }
        Geometry::MultiLineString(mls) => {
            for ls in mls.iter() {
                out.extend(cover_line_string(ls, &RegionCoverer::at_level(level))?);
            }
        }
    }
    Some(())
}

fn collect_adaptive(
    geom: &Geometry,
    coverer: &RegionCoverer,
    out: &mut Vec<S2CellId>,
) -> Option<()> {
    let point_level = coverer.max_level.min(MAX_LEVEL);
    match geom {
        Geometry::Point(p) => {
            if valid_lonlat(p.x(), p.y()) {
                let gp = GeoPoint::from_lon_lat(p.x(), p.y());
                out.push(S2CellId::from_geo_point_at_level(gp, point_level));
            }
        }
        Geometry::MultiPoint(mp) => {
            for p in mp.iter() {
                if valid_lonlat(p.x(), p.y()) {
                    let gp = GeoPoint::from_lon_lat(p.x(), p.y());
                    out.push(S2CellId::from_geo_point_at_level(gp, point_level));
                }
            }
        }
        Geometry::Polygon(poly) => out.extend(adaptive_polygon(poly, coverer)?.iter()),
        Geometry::MultiPolygon(mp) => {
            for poly in mp.iter() {
                out.extend(adaptive_polygon(poly, coverer)?.iter());
            }
        }
        Geometry::Rect(r) => out.extend(adaptive_polygon(&r.to_polygon(), coverer)?.iter()),
        Geometry::Triangle(t) => out.extend(adaptive_polygon(&t.to_polygon(), coverer)?.iter()),
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                collect_adaptive(g, coverer, out)?;
            }
        }
        Geometry::Line(l) => {
            let coords = [(l.start.x, l.start.y), (l.end.x, l.end.y)];
            out.extend(cover_line_coords(&coords, coverer)?);
        }
        Geometry::LineString(ls) => out.extend(cover_line_string(ls, coverer)?),
        Geometry::MultiLineString(mls) => {
            for ls in mls.iter() {
                out.extend(cover_line_string(ls, coverer)?);
            }
        }
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
    let s2 = polygon_region(poly)?;
    out.extend(RegionCoverer::at_level(level).covering(&s2));
    Some(())
}

fn adaptive_polygon(poly: &Polygon, coverer: &RegionCoverer) -> Option<CellUnion> {
    Some(coverer.adaptive_covering(&polygon_region(poly)?))
}

/// The lon/lat polygon as a spherical [`S2Polygon`]; `None` when a ring has an
/// invalid coordinate.
fn polygon_region(poly: &Polygon) -> Option<S2Polygon> {
    for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
        if !ring.coords().all(|c| valid_lonlat(c.x, c.y)) {
            return None;
        }
    }
    let shell = ring_loop(poly.exterior())?;
    let holes: Vec<Loop> = poly.interiors().iter().filter_map(ring_loop).collect();
    Some(S2Polygon::new(shell, holes))
}

/// The cells covering a chain of lon/lat coordinates: empty for no
/// coordinates, the containing cell for one, and the [`Polyline`] cover
/// otherwise; `None` when any coordinate is invalid.
fn cover_line_coords(coords: &[(f64, f64)], coverer: &RegionCoverer) -> Option<Vec<S2CellId>> {
    if !coords.iter().all(|&(lon, lat)| valid_lonlat(lon, lat)) {
        return None;
    }
    let level = coverer.max_level.min(MAX_LEVEL);
    match coords {
        [] => Some(Vec::new()),
        [(lon, lat)] => {
            let gp = GeoPoint::from_lon_lat(*lon, *lat);
            Some(vec![S2CellId::from_geo_point_at_level(gp, level)])
        }
        _ => {
            let pts: Vec<S2Point> = coords
                .iter()
                .map(|&(lon, lat)| S2Point::from_geo_point(GeoPoint::from_lon_lat(lon, lat)))
                .collect();
            let line = Polyline::new(pts);
            if coverer.min_level == coverer.max_level && coverer.max_cells == usize::MAX {
                Some(coverer.covering(&line))
            } else {
                Some(coverer.adaptive_covering(&line).cell_ids().to_vec())
            }
        }
    }
}

fn cover_line_string(ls: &LineString, coverer: &RegionCoverer) -> Option<Vec<S2CellId>> {
    let coords: Vec<(f64, f64)> = ls.coords().map(|c| (c.x, c.y)).collect();
    cover_line_coords(&coords, coverer)
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
    use hausdorff_planar::space::r2::{Line, LineString, MultiLineString, Point, Polygon};

    fn at(lon: f64, lat: f64, level: u32) -> S2CellId {
        S2CellId::from_geo_point_at_level(GeoPoint::from_lon_lat(lon, lat), level)
    }

    #[test]
    fn point_polygon_and_line_covers() {
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
        let line_cover = cover_r2_geometry_at_level(&line, 10).unwrap();
        assert!(!line_cover.is_empty());
        for c in &line_cover {
            assert_eq!(c.level(), 10);
        }
    }

    #[test]
    fn line_cover_contains_sampled_point_ancestors() {
        let level = 10;
        let line = Geometry::LineString(LineString::from(vec![
            (10.0, 50.0),
            (10.6, 50.3),
            (11.0, 49.9),
        ]));
        let cover = cover_r2_geometry_at_level(&line, level).unwrap();
        for seg in [((10.0, 50.0), (10.6, 50.3)), ((10.6, 50.3), (11.0, 49.9))] {
            for k in 0..=300 {
                let t = f64::from(k) / 300.0;
                let lon = seg.0.0 + (seg.1.0 - seg.0.0) * t;
                let lat = seg.0.1 + (seg.1.1 - seg.0.1) * t;
                let ancestor = at(lon, lat, level);
                assert!(
                    cover.contains(&ancestor),
                    "sampled point ({lon},{lat}) not covered at level {level}"
                );
            }
        }
    }

    #[test]
    fn line_and_multilinestring_geometries_are_covered() {
        let seg = Geometry::Line(Line::new((10.0, 50.0), (10.2, 50.1)));
        let cover = cover_r2_geometry_at_level(&seg, 12).unwrap();
        assert!(cover.contains(&at(10.0, 50.0, 12)));
        assert!(cover.contains(&at(10.2, 50.1, 12)));

        let mls = Geometry::MultiLineString(MultiLineString::new(vec![
            LineString::from(vec![(0.0, 0.0), (0.5, 0.0)]),
            LineString::from(vec![(20.0, 10.0), (20.0, 10.5)]),
        ]));
        let cover = cover_r2_geometry_at_level(&mls, 9).unwrap();
        assert!(cover.contains(&at(0.2, 0.0, 9)));
        assert!(cover.contains(&at(20.0, 10.2, 9)));

        let degenerate = Geometry::LineString(LineString::from(Vec::<(f64, f64)>::new()));
        assert_eq!(cover_r2_geometry_at_level(&degenerate, 9).unwrap(), vec![]);
        let single = Geometry::LineString(LineString::from(vec![(10.0, 50.0)]));
        assert_eq!(
            cover_r2_geometry_at_level(&single, 9).unwrap(),
            vec![at(10.0, 50.0, 9)]
        );
    }

    #[test]
    fn adaptive_cover_of_collection() {
        let geom = Geometry::GeometryCollection(vec![
            Geometry::Polygon(Polygon::new(
                LineString::from(vec![
                    (9.5, 49.5),
                    (10.5, 49.5),
                    (10.5, 50.5),
                    (9.5, 50.5),
                    (9.5, 49.5),
                ]),
                vec![],
            )),
            Geometry::LineString(LineString::from(vec![(12.0, 50.0), (13.0, 50.5)])),
            Geometry::Point(Point::new(8.0, 48.0)),
        ]);
        let coverer = RegionCoverer {
            min_level: 0,
            max_level: 12,
            max_cells: 24,
        };
        let union = cover_r2_geometry(&geom, &coverer).unwrap();
        assert!(!union.is_empty());
        let p = |lon: f64, lat: f64| S2Point::from_geo_point(GeoPoint::from_lon_lat(lon, lat));
        assert!(union.contains_point(p(10.0, 50.0)), "polygon interior");
        assert!(union.contains_point(p(12.5, 50.25)), "line midpoint");
        assert!(union.contains_point(p(8.0, 48.0)), "point member");
        assert!(!union.contains_point(p(-60.0, -30.0)));
    }

    #[test]
    fn level_above_leaf_is_none() {
        let pt = Geometry::Point(Point::new(10.0, 50.0));
        assert!(cover_r2_geometry_at_level(&pt, 31).is_none());
    }

    #[test]
    fn out_of_range_coords_are_rejected() {
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

        let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (200.0, 0.0)]));
        assert!(cover_r2_geometry_at_level(&line, 8).is_none());
        assert!(cover_r2_geometry(&line, &RegionCoverer::default()).is_none());
    }
}
