//! Native WKB codec: round-trips and interop with the vendored `wkb` crate (so
//! it can be deprecated in favour of the native one).

use hausdorff::r2_format::wkb::{from_wkb, to_wkb};
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
fn interops_with_wkb_crate() {
    for g in fixtures() {
        // native bytes decode via the wkb crate (through compat) back to the same geometry
        assert_eq!(hausdorff::compat::decode_wkb(&to_wkb(&g)).unwrap(), g);
        // the wkb crate's bytes decode via the native reader back to the same geometry
        assert_eq!(from_wkb(&hausdorff::compat::encode_wkb(&g).unwrap()).unwrap(), g);
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
