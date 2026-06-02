//! Native WKB codec: round-trips and interop with the vendored `wkb` crate (so
//! it can be deprecated in favour of the native one).

use hausdorff::format::wkb::{from_wkb, to_wkb};
use hausdorff::{Geometry, LineString, MultiPolygon, Point, Polygon};

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
        assert_eq!(from_wkb(&to_wkb(&g)).unwrap(), g);
    }
}

#[test]
fn matches_wkb_reference_crate() {
    use geo_traits::to_geo::ToGeoGeometry;
    use wkb::writer::{WriteOptions, write_geometry};
    for g in fixtures() {
        // reference bytes (the `wkb` crate, via geo_types) decode through native to g
        let mut reference = Vec::new();
        write_geometry(
            &mut reference,
            &common::to_geo_geometry(&g),
            &WriteOptions::default(),
        )
        .unwrap();
        assert_eq!(from_wkb(&reference).unwrap(), g, "native decode of reference WKB");
        // native bytes decode through the `wkb` crate back to g
        let decoded = wkb::reader::read_wkb(&to_wkb(&g)).unwrap().to_geometry();
        assert_eq!(
            common::from_geo_geometry(decoded),
            g,
            "reference decode of native WKB"
        );
    }
}

#[test]
fn reads_big_endian() {
    // a big-endian (XDR) POINT(1 2): order=0, type=1, then two big-endian f64s
    let mut bytes = vec![0u8, 0, 0, 0, 1];
    bytes.extend_from_slice(&1.0f64.to_be_bytes());
    bytes.extend_from_slice(&2.0f64.to_be_bytes());
    assert_eq!(from_wkb(&bytes).unwrap(), Geometry::Point(Point::new(1.0, 2.0)));
}

#[test]
fn from_wkb_rejects_trailing_bytes() {
    use hausdorff::format::wkb::from_wkb_prefix;
    let g = Geometry::Point(Point::new(1.0, 2.0));
    let mut bytes = to_wkb(&g);
    let exact_len = bytes.len();
    // exact bytes decode
    assert_eq!(from_wkb(&bytes).unwrap(), g);
    // a trailing byte is rejected by the strict decoder
    bytes.push(0xFF);
    assert!(from_wkb(&bytes).is_err(), "trailing byte must be rejected");
    // but the prefix decoder parses the geometry and reports bytes consumed
    let (decoded, consumed) = from_wkb_prefix(&bytes).unwrap();
    assert_eq!(decoded, g);
    assert_eq!(consumed, exact_len);
}
