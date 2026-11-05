//! Native GeoJSON codec: round-trips, standard strings, and conformance against
//! golden GeoJSON (captured from the `geojson` crate at the time the native codec
//! replaced it; the live dependency is no longer needed).

use crate::format::geojson::{from_geojson, to_geojson};
use crate::space::r2::{Geometry, LineString, MultiPolygon, Point, Polygon};

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
        Geometry::MultiPolygon(MultiPolygon::new(vec![poly(vec![vec![
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
        ]])])),
        Geometry::GeometryCollection(vec![Geometry::Point(Point::new(-3.0, 7.0))]),
    ]
}

#[test]
fn native_round_trips() {
    for g in fixtures() {
        assert_eq!(from_geojson(&to_geojson(&g)).unwrap(), g);
    }
}

#[test]
fn parses_standard_strings() {
    let g = from_geojson(r#"{"type":"Point","coordinates":[3,4]}"#).unwrap();
    assert_eq!(g, Geometry::Point(Point::new(3.0, 4.0)));
    let g = from_geojson(r#"{"type":"Point","coordinates":[3,4,99]}"#).unwrap();
    assert_eq!(g, Geometry::Point(Point::new(3.0, 4.0)));
    assert!(from_geojson(r#"{"type":"Nope","coordinates":[]}"#).is_err());
    assert!(from_geojson("not json").is_err());
}

#[test]
fn parses_golden_reference_geojson() {
    // Standard GeoJSON for each fixture (golden captured from the `geojson` crate).
    const GOLDEN: [&str; 5] = [
        r#"{"type":"Point","coordinates":[1.0,2.5]}"#,
        r#"{"type":"LineString","coordinates":[[0.0,0.0],[1.0,1.0],[2.0,0.0]]}"#,
        r#"{"type":"Polygon","coordinates":[[[0.0,0.0],[4.0,0.0],[4.0,4.0],[0.0,4.0],[0.0,0.0]],[[1.0,1.0],[2.0,1.0],[2.0,2.0],[1.0,2.0],[1.0,1.0]]]}"#,
        r#"{"type":"MultiPolygon","coordinates":[[[[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,0.0]]]]}"#,
        r#"{"type":"GeometryCollection","geometries":[{"type":"Point","coordinates":[-3.0,7.0]}]}"#,
    ];
    for (g, golden) in fixtures().into_iter().zip(GOLDEN) {
        assert_eq!(
            from_geojson(golden).unwrap(),
            g,
            "native parse of `{golden}`"
        );
        // And the native encoder still round-trips its own output.
        assert_eq!(from_geojson(&to_geojson(&g)).unwrap(), g);
    }
}
