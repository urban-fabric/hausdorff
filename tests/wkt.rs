//! Native WKT codec: round-trips, standard strings, and interop with the
//! vendored `wkt` crate (so it can be deprecated in favour of the native one).

use hausdorff::format::wkt::{from_wkt, to_wkt};
use hausdorff::{Coord, Geometry, LineString, MultiPolygon, Point, Polygon};

mod common;

fn poly(rings: Vec<Vec<(f64, f64)>>) -> Polygon {
    let mut it = rings.into_iter();
    let ext = LineString::from(it.next().unwrap());
    Polygon::new(ext, it.map(LineString::from).collect())
}

fn fixtures() -> Vec<Geometry> {
    vec![
        Geometry::Point(Point::new(1.0, 2.5)),
        Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])),
        Geometry::Polygon(poly(vec![
            vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
            vec![(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)],
        ])),
        Geometry::MultiPolygon(MultiPolygon::new(vec![
            poly(vec![vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]]),
            poly(vec![vec![(5.0, 5.0), (6.0, 5.0), (6.0, 6.0)]]),
        ])),
        Geometry::GeometryCollection(vec![
            Geometry::Point(Point::new(-3.0, 7.0)),
            Geometry::LineString(LineString::from(vec![(0.0, 0.0), (10.0, 10.0)])),
        ]),
    ]
}

#[test]
fn native_round_trips() {
    for g in fixtures() {
        let s = to_wkt(&g);
        assert_eq!(from_wkt(&s).unwrap(), g, "round trip via `{s}`");
    }
}

#[test]
fn parses_standard_strings() {
    let g = from_wkt("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))").unwrap();
    let Geometry::Polygon(p) = g else { panic!() };
    assert_eq!(p.exterior().coords_count(), 5);
    assert!(matches!(
        from_wkt("POINT (3 4)").unwrap(),
        Geometry::Point(_)
    ));
    assert!(from_wkt("LINESTRING EMPTY").is_ok());
    assert!(from_wkt("MULTIPOINT (1 2, 3 4)").is_ok()); // bare form
    assert!(from_wkt("NONSENSE (1 2)").is_err());
    let _ = Coord::new(0.0, 0.0);
}

#[test]
fn matches_wkt_reference_crate() {
    use wkt::{ToWkt, TryFromWkt};
    for g in fixtures() {
        // reference WKT (the `wkt` crate, via geo_types) parses through native to g
        let reference = common::to_geo_geometry(&g).wkt_string();
        assert_eq!(
            from_wkt(&reference).unwrap(),
            g,
            "native parse of reference WKT"
        );
        // native WKT parses through the `wkt` crate back to g
        let parsed = geo_types::Geometry::<f64>::try_from_wkt_str(&to_wkt(&g)).unwrap();
        assert_eq!(
            common::from_geo_geometry(parsed),
            g,
            "reference parse of native WKT"
        );
    }
}
