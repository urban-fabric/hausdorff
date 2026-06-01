//! Native GeoJSON codec: round-trips, standard strings, and interop with the
//! vendored `geojson` crate (so it can be deprecated).

use hausdorff::r2_format::geojson::{from_geojson, to_geojson};
use hausdorff::{Geometry, LineString, MultiPolygon, Point, Polygon};

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
        ])),
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
    // elevation (3rd ordinate) is ignored
    let g = from_geojson(r#"{"type":"Point","coordinates":[3,4,99]}"#).unwrap();
    assert_eq!(g, Geometry::Point(Point::new(3.0, 4.0)));
    assert!(from_geojson(r#"{"type":"Nope","coordinates":[]}"#).is_err());
    assert!(from_geojson("not json").is_err());
}

#[test]
fn interops_with_geojson_crate() {
    for g in fixtures() {
        // native output parses via the geojson crate (through compat) back to g
        assert_eq!(hausdorff::compat::geojson_to_geometry(&to_geojson(&g)).unwrap(), g);
        // the geojson crate's output parses via native back to g
        assert_eq!(from_geojson(&hausdorff::compat::geometry_to_geojson(&g)).unwrap(), g);
    }
}
