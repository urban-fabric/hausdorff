//! Native WKB codec: round-trips and conformance against golden WKB bytes
//! (captured from the `wkb` crate at the time the native codec replaced it; the
//! live dependency is no longer needed).

use crate::format::wkb::{from_wkb, to_wkb};
use crate::space::r2::{Geometry, LineString, MultiPolygon, Point, Polygon};

fn hex(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
        .collect()
}

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
fn decodes_golden_reference_wkb() {
    // Little-endian WKB for each fixture (golden captured from the `wkb` crate).
    const GOLDEN: [&str; 5] = [
        "0101000000000000000000f03f0000000000000440",
        "01020000000300000000000000000000000000000000000000000000000000f03f000000000000f03f00000000000000400000000000000000",
        "01030000000200000005000000000000000000000000000000000000000000000000001040000000000000000000000000000010400000000000001040000000000000000000000000000010400000000000000000000000000000000005000000000000000000f03f000000000000f03f0000000000000040000000000000f03f00000000000000400000000000000040000000000000f03f0000000000000040000000000000f03f000000000000f03f",
        "0106000000020000000103000000010000000400000000000000000000000000000000000000000000000000f03f0000000000000000000000000000f03f000000000000f03f000000000000000000000000000000000103000000010000000400000000000000000014400000000000001440000000000000184000000000000014400000000000001840000000000000184000000000000014400000000000001440",
        "010700000002000000010100000000000000000008c00000000000001c400102000000020000000000000000000000000000000000000000000000000024400000000000002440",
    ];
    for (g, golden) in fixtures().into_iter().zip(GOLDEN) {
        assert_eq!(
            from_wkb(&hex(golden)).unwrap(),
            g,
            "native decode of golden WKB"
        );
        // And the native encoder still round-trips its own output.
        assert_eq!(from_wkb(&to_wkb(&g)).unwrap(), g);
    }
}

#[test]
fn reads_big_endian() {
    let mut bytes = vec![0u8, 0, 0, 0, 1];
    bytes.extend_from_slice(&1.0f64.to_be_bytes());
    bytes.extend_from_slice(&2.0f64.to_be_bytes());
    assert_eq!(
        from_wkb(&bytes).unwrap(),
        Geometry::Point(Point::new(1.0, 2.0))
    );
}

#[test]
fn from_wkb_rejects_trailing_bytes() {
    use crate::format::wkb::from_wkb_prefix;
    let g = Geometry::Point(Point::new(1.0, 2.0));
    let mut bytes = to_wkb(&g);
    let exact_len = bytes.len();
    assert_eq!(from_wkb(&bytes).unwrap(), g);
    bytes.push(0xFF);
    assert!(from_wkb(&bytes).is_err(), "trailing byte must be rejected");
    let (decoded, consumed) = from_wkb_prefix(&bytes).unwrap();
    assert_eq!(decoded, g);
    assert_eq!(consumed, exact_len);
}
