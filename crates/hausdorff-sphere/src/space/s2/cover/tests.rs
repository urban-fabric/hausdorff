//! S2 covering and bounding tests.

use super::*;

use super::*;
use hausdorff_planar::space::r2::{Line, LineString, MultiLineString, Point, Polygon};

fn at(lon: f64, lat: f64, level: u32) -> S2CellId {
    S2CellId::from_geo_point_at_level(GeoPoint::from_lon_lat(lon, lat), level)
}

#[test]
fn bound_matches_planar_rect_away_from_antimeridian() {
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
    let bound = bound_r2_geometry(&poly).unwrap();
    assert!(!bound.is_lng_inverted());
    assert!((bound.lo().lon - 9.5).abs() < 1e-9 && (bound.hi().lon - 10.5).abs() < 1e-9);
    assert!((bound.lo().lat - 49.5).abs() < 1e-9 && (bound.hi().lat - 50.5).abs() < 1e-9);
    assert!(bound.contains_lat_lng(GeoPoint::from_lon_lat(10.0, 50.0)));
}

#[test]
fn bound_wraps_antimeridian_polygon() {
    // A polygon spanning 178° → -178° (4° wide across ±180°). A planar min/max
    // box would report lng [-178, 178] (356° wide); the wrap-aware bound is 4°.
    let poly = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (178.0, 10.0),
            (-178.0, 10.0),
            (-178.0, 20.0),
            (178.0, 20.0),
            (178.0, 10.0),
        ]),
        vec![],
    ));
    let bound = bound_r2_geometry(&poly).unwrap();
    assert!(bound.is_lng_inverted(), "antimeridian polygon must wrap");
    assert!(bound.contains_lat_lng(GeoPoint::from_lon_lat(180.0, 15.0)));
    assert!(bound.contains_lat_lng(GeoPoint::from_lon_lat(-179.0, 15.0)));
    assert!(!bound.contains_lat_lng(GeoPoint::from_lon_lat(0.0, 15.0)));
}

#[test]
fn bound_none_for_empty_or_invalid() {
    let empty = Geometry::MultiPoint(hausdorff_planar::space::r2::MultiPoint::new(vec![]));
    assert!(bound_r2_geometry(&empty).is_none());
    let bad = Geometry::Point(Point::new(200.0, 100.0)); // out of range
    assert!(bound_r2_geometry(&bad).is_none());
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
