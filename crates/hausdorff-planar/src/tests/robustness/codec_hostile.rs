//! Codec behaviour on **hostile** input: the codecs are the only place in this
//! workspace that reads bytes it did not produce, so they are the only place where
//! "what does this do on garbage" is a correctness question rather than a
//! hypothetical.
//!
//! The bar is narrow and absolute: **every input must produce `Ok` or `Err`.** Not
//! a panic, and in particular not an abort. That last distinction is the point of
//! this file. A panic unwinds and a caller can catch it; the two failure modes
//! found here could not be caught by anything:
//!
//! - **Stack overflow from unbounded recursion.** WKB and WKT both nest through
//!   their collection types, and both decoders recursed to match with no limit.
//!   Nine bytes bought one level of WKB nesting, so ~180 KB of input killed the
//!   process with `fatal runtime error: stack overflow`.
//! - **Allocation abort from an unvalidated count.** Both WKB and the shapefile
//!   reader passed a count read straight off the wire to `Vec::with_capacity`.
//!   Nine bytes claiming `u32::MAX` members asked the allocator for 38–240 GB;
//!   when that is refused the process aborts.
//!
//! Both are now bounded, and the tests below are what keep them bounded. They are
//! written to *fail by crashing* if a guard is removed, which is unusual and worth
//! knowing: a regression here does not show up as a red assertion, it shows up as
//! the test binary dying and taking the rest of the suite with it.
//!
//! The `Vec::with_capacity` risk is subtler than it looks, which is why the counts
//! are checked at every container type separately. `coords` builds its vector with
//! `(0..n).map(..).collect::<Result<Vec<_>, _>>()`, and that was *already* safe
//! before the fix — the `Result` shunt iterator reports a `size_hint` lower bound
//! of 0, so `collect` reserves nothing. The four explicit `with_capacity(n)` calls
//! in the multi-geometry and collection arms were not. Two code paths that look
//! identical in a skim differed in whether they could abort the process, so each
//! arm gets its own case rather than one representative.

use crate::format::geojson::from_geojson;
use crate::format::shapefile::read_geometries;
use crate::format::wkb::{from_wkb, from_wkb_prefix, to_wkb};
use crate::format::wkt::{from_wkt, to_wkt};
use crate::space::r2::{Coord, Geometry, LineString, Point, Polygon};

/// A little-endian WKB header followed by a 32-bit count and nothing else — the
/// nine-byte shape of the allocation bomb.
fn header_with_count(type_code: u32, count: u32) -> Vec<u8> {
    let mut b = vec![0x01u8];
    b.extend_from_slice(&type_code.to_le_bytes());
    b.extend_from_slice(&count.to_le_bytes());
    b
}

#[test]
fn wkb_absurd_element_counts_are_rejected_without_allocating() {
    // 2 = LineString (coordinates), 4/5/6 = Multi*, 7 = GeometryCollection,
    // 3 = Polygon (ring count). Every one of these reaches a sizing decision made
    // from a number the input chose.
    for type_code in [2u32, 3, 4, 5, 6, 7] {
        for count in [u32::MAX, u32::MAX / 2, 1 << 30, 1 << 20] {
            let bytes = header_with_count(type_code, count);
            let got = from_wkb(&bytes);
            assert!(
                got.is_err(),
                "type {type_code} claiming {count} elements in {} bytes must be rejected, got {got:?}",
                bytes.len()
            );
        }
    }
}

#[test]
fn wkb_nesting_is_bounded() {
    // Each level is a GeometryCollection header plus a member count of 1.
    let level = header_with_count(7, 1);
    for depth in [100usize, 10_000, 200_000] {
        let mut bytes = level.repeat(depth);
        bytes.push(0x01);
        bytes.extend_from_slice(&1u32.to_le_bytes());
        bytes.extend_from_slice(&0f64.to_le_bytes());
        bytes.extend_from_slice(&0f64.to_le_bytes());
        let got = from_wkb(&bytes);
        assert!(
            got.is_err(),
            "{depth} levels of nesting must be refused, not followed"
        );
    }
}

#[test]
fn wkb_shallow_nesting_still_parses() {
    // The depth guard must not have cost us legitimate nesting. Two levels is
    // already more than any real producer emits.
    let inner = to_wkb(&Geometry::Point(Point(Coord::new(1.0, 2.0))));
    let mut lvl1 = header_with_count(7, 1);
    lvl1.extend_from_slice(&inner);
    let mut lvl2 = header_with_count(7, 1);
    lvl2.extend_from_slice(&lvl1);
    match from_wkb(&lvl2) {
        Ok(Geometry::GeometryCollection(gs)) => assert_eq!(gs.len(), 1),
        other => panic!("two-level collection should round-trip, got {other:?}"),
    }
}

#[test]
fn wkb_truncated_at_every_prefix_is_rejected_not_panicked() {
    // A valid encoding cut at every possible point. Each prefix is either a
    // complete shorter geometry or an error; none may panic.
    let geoms = vec![
        Geometry::Point(Point(Coord::new(1.0, 2.0))),
        Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.5)])),
        Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (4.0, 0.0),
                (4.0, 4.0),
                (0.0, 4.0),
                (0.0, 0.0),
            ]),
            vec![LineString::from(vec![
                (1.0, 1.0),
                (2.0, 1.0),
                (2.0, 2.0),
                (1.0, 1.0),
            ])],
        )),
    ];
    for g in &geoms {
        let full = to_wkb(g);
        for cut in 0..full.len() {
            // Must not panic. Truncation should be an error, though a prefix that
            // happens to be a whole geometry is legitimately Ok for the prefix API.
            let _ = from_wkb(&full[..cut]);
            let _ = from_wkb_prefix(&full[..cut]);
        }
        assert!(
            from_wkb(&full).is_ok(),
            "the untruncated encoding must parse"
        );
    }
}

#[test]
fn wkb_bad_byte_order_and_unknown_types_are_rejected() {
    for order in [2u8, 0x7f, 0xff] {
        let mut b = vec![order];
        b.extend_from_slice(&1u32.to_le_bytes());
        b.extend_from_slice(&0f64.to_le_bytes());
        b.extend_from_slice(&0f64.to_le_bytes());
        assert!(from_wkb(&b).is_err(), "byte order {order} is not 0 or 1");
    }
    // 0 and 8..=999 are not OGC geometry types; the decoder takes `raw % 1000` as
    // the base type, so 8 and 1008 must both be refused rather than aliasing.
    for type_code in [0u32, 8, 99, 1008, 0x0000_ffff] {
        let mut b = vec![0x01u8];
        b.extend_from_slice(&type_code.to_le_bytes());
        b.extend_from_slice(&0f64.to_le_bytes());
        b.extend_from_slice(&0f64.to_le_bytes());
        assert!(
            from_wkb(&b).is_err(),
            "geometry type {type_code} should be unsupported"
        );
    }
}

#[test]
fn wkb_empty_and_tiny_buffers_are_rejected() {
    for len in 0..5usize {
        let b = vec![0x01u8; len];
        assert!(from_wkb(&b).is_err(), "a {len}-byte buffer cannot be WKB");
    }
}

#[test]
fn wkb_non_finite_coordinates_round_trip_bit_exactly() {
    // NaN and the infinities are representable in WKB and the codec is a codec, not
    // a validator: it must carry them through rather than silently repairing them.
    // Pinned as bits because NaN != NaN.
    for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.0] {
        let g = Geometry::Point(Point(Coord::new(v, v)));
        match from_wkb(&to_wkb(&g)) {
            Ok(Geometry::Point(p)) => {
                assert_eq!(
                    p.0.x.to_bits(),
                    v.to_bits(),
                    "x should survive bit-exactly for {v}"
                );
                assert_eq!(p.0.y.to_bits(), v.to_bits());
            }
            other => panic!("expected a point, got {other:?}"),
        }
    }
}

#[test]
fn wkt_nesting_is_bounded() {
    for depth in [100usize, 10_000, 200_000] {
        let text = format!(
            "{}POINT(0 0){}",
            "GEOMETRYCOLLECTION(".repeat(depth),
            ")".repeat(depth)
        );
        assert!(
            from_wkt(&text).is_err(),
            "{depth} levels of WKT nesting must be refused, not followed"
        );
    }
    // Still accepts nesting a real producer might emit.
    assert!(from_wkt("GEOMETRYCOLLECTION(GEOMETRYCOLLECTION(POINT(1 2)))").is_ok());
}

#[test]
fn wkt_malformed_text_is_rejected_not_panicked() {
    let cases = [
        "",
        " ",
        "\0",
        "POINT",
        "POINT(",
        "POINT()",
        "POINT(1)",
        "POINT(1 2",
        "POINT 1 2)",
        "POINT(1 2))",
        "POINT(1 2) POINT(3 4)",
        "POINT(a b)",
        "POINT(1,2)",
        "NOTAGEOMETRY(1 2)",
        "LINESTRING",
        "LINESTRING(",
        "LINESTRING(1 2",
        "LINESTRING(1 2,)",
        "LINESTRING(,)",
        "POLYGON",
        "POLYGON(())",
        "POLYGON((1 2))",
        "POLYGON((1 2,3 4),)",
        "GEOMETRYCOLLECTION",
        "GEOMETRYCOLLECTION(",
        "GEOMETRYCOLLECTION()",
        "GEOMETRYCOLLECTION(POINT(1 2),)",
        "MULTIPOINT",
        "MULTIPOINT((1 2),)",
        // numeric edges the lexer has to survive
        "POINT(1e 2)",
        "POINT(1e999999 2)",
        "POINT(- 2)",
        "POINT(. .)",
        "POINT(1.2.3 4)",
        "POINT(--1 2)",
        "POINT(0x10 2)",
        // non-ASCII, which the keyword scanner must not slice through
        "PÖINT(1 2)",
        "POINT(1 2)\u{0}",
        "POINT（1 2）",
    ];
    for text in cases {
        // Only requirement: it returns. A panic here fails the test.
        let _ = from_wkt(text);
    }
    // And a sanity anchor: the well-formed neighbours of those cases do parse.
    assert!(from_wkt("POINT(1 2)").is_ok());
    assert!(from_wkt("LINESTRING(1 2,3 4)").is_ok());
    assert!(from_wkt("POLYGON((0 0,1 0,1 1,0 0))").is_ok());
}

#[test]
fn wkt_round_trips_through_its_own_writer() {
    // `to_wkt` output must be accepted by `from_wkt` — a writer and parser that
    // disagree is how a format quietly becomes a private dialect.
    let geoms = vec![
        Geometry::Point(Point(Coord::new(1.5, -2.5))),
        Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0)])),
        Geometry::Polygon(Polygon::new(
            LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 0.0)]),
            vec![],
        )),
        Geometry::GeometryCollection(vec![Geometry::Point(Point(Coord::new(7.0, 8.0)))]),
    ];
    for g in &geoms {
        let text = to_wkt(g);
        assert!(
            from_wkt(&text).is_ok(),
            "our own writer produced WKT our parser rejects: {text}"
        );
    }
}

#[test]
fn geojson_malformed_text_is_rejected_not_panicked() {
    let cases = [
        "",
        "{}",
        "null",
        "[]",
        "1",
        r#"{"type":"Point"}"#,
        r#"{"type":"Point","coordinates":null}"#,
        r#"{"type":"Point","coordinates":[]}"#,
        r#"{"type":"Point","coordinates":[1]}"#,
        r#"{"type":"Point","coordinates":["a","b"]}"#,
        r#"{"type":"Point","coordinates":{"x":1,"y":2}}"#,
        r#"{"type":"Nonsense","coordinates":[1,2]}"#,
        r#"{"type":"LineString","coordinates":[[1,2]"#,
        r#"{"type":"Polygon","coordinates":[[[1,2]]]}"#,
        r#"{"type":"GeometryCollection"}"#,
        r#"{"type":"GeometryCollection","geometries":"no"}"#,
        r#"{"coordinates":[1,2]}"#,
        "{\"type\":\"Point\",\"coordinates\":[1e999999,2]}",
    ];
    for text in cases {
        let _ = from_geojson(text);
    }
    assert!(from_geojson(r#"{"type":"Point","coordinates":[1,2]}"#).is_ok());
}

#[test]
fn geojson_nesting_is_bounded() {
    // serde_json enforces its own recursion limit, so this has always been safe —
    // pinned so that a future hand-rolled JSON layer cannot regress it silently.
    let head = r#"{"type":"GeometryCollection","geometries":["#;
    let text = format!(
        "{}{}{}",
        head.repeat(50_000),
        r#"{"type":"Point","coordinates":[0,0]}"#,
        "]}".repeat(50_000)
    );
    assert!(from_geojson(&text).is_err(), "deep nesting must be refused");
}

#[test]
fn shapefile_hostile_bytes_are_rejected_not_panicked() {
    // Truncated and malformed headers.
    for len in 0..120usize {
        let _ = read_geometries(&vec![0u8; len][..]);
    }
    // A plausible 100-byte header whose record then claims an absurd point count.
    // This is the shapefile shape of the same allocation bomb as WKB's.
    let mut buf = vec![0u8; 100];
    buf[0..4].copy_from_slice(&9994i32.to_be_bytes()); // file code
    buf[28..32].copy_from_slice(&1000i32.to_le_bytes()); // version
    buf[32..36].copy_from_slice(&5i32.to_le_bytes()); // polygon
    // record header: number 1, content length in 16-bit words
    buf.extend_from_slice(&1i32.to_be_bytes());
    buf.extend_from_slice(&1000i32.to_be_bytes());
    buf.extend_from_slice(&5i32.to_le_bytes()); // shape type: polygon
    buf.extend_from_slice(&[0u8; 32]); // bbox
    buf.extend_from_slice(&i32::MAX.to_le_bytes()); // num_parts
    buf.extend_from_slice(&i32::MAX.to_le_bytes()); // num_points
    let got = read_geometries(&buf[..]);
    assert!(
        got.is_err(),
        "a record claiming i32::MAX parts and points must be rejected, got {got:?}"
    );
}
