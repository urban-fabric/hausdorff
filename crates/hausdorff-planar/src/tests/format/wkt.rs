//! Native WKT codec: round-trips, standard strings, and conformance against
//! golden standard-form WKT (captured from the `wkt` crate at the time the
//! native codec replaced it; the live dependency is no longer needed).

use crate::format::wkt::{from_wkt, to_wkt};
use crate::space::r2::{Coord, Geometry, LineString, MultiPolygon, Point, Polygon};

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
    assert!(from_wkt("MULTIPOINT (1 2, 3 4)").is_ok());
    assert!(from_wkt("NONSENSE (1 2)").is_err());
    let _ = Coord::new(0.0, 0.0);
}

#[test]
fn parses_golden_standard_wkt() {
    // Standard-form WKT for each fixture (golden captured from the `wkt` crate).
    const GOLDEN: [&str; 5] = [
        "POINT(1 2.5)",
        "LINESTRING(0 0,1 1,2 0)",
        "POLYGON((0 0,4 0,4 4,0 4,0 0),(1 1,2 1,2 2,1 2,1 1))",
        "MULTIPOLYGON(((0 0,1 0,1 1,0 0)),((5 5,6 5,6 6,5 5)))",
        "GEOMETRYCOLLECTION(POINT(-3 7),LINESTRING(0 0,10 10))",
    ];
    for (g, golden) in fixtures().into_iter().zip(GOLDEN) {
        assert_eq!(from_wkt(golden).unwrap(), g, "native parse of `{golden}`");
        // And the native encoder still round-trips its own output.
        assert_eq!(from_wkt(&to_wkt(&g)).unwrap(), g);
    }
}
