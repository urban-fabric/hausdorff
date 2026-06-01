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
fn matches_geojson_reference_crate() {
    for g in fixtures() {
        let geo_g = hausdorff::compat::to_geo_geometry(&g);
        // reference GeoJSON (the `geojson` crate) parses through native to g
        let reference = geojson::GeometryValue::from(&geo_g).to_string();
        assert_eq!(from_geojson(&reference).unwrap(), g, "native parse of reference GeoJSON");
        // native GeoJSON parses through the `geojson` crate back to g
        let parsed: geojson::Geometry = to_geojson(&g).parse().unwrap();
        let geo_parsed = geo_types::Geometry::<f64>::try_from(parsed.value).unwrap();
        assert_eq!(
            hausdorff::compat::from_geo_geometry(geo_parsed),
            g,
            "reference parse of native GeoJSON"
        );
    }
}
