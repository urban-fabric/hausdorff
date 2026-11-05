//! CRS interchange-format conformance, split by the format under test.
//!
//! Each of the three formats encodes the same CRS model, so each file asks the
//! same two questions of its own syntax: does parsing a string produced by an
//! independent implementation (pyproj) yield a CRS that projects identically to
//! the registry's, and are malformed inputs rejected rather than silently
//! half-parsed.
//!
//! | file | format |
//! |---|---|
//! | `proj_string` | PROJ-strings (`+proj=… +datum=…`) |
//! | `wkt` | WKT2:2019 (ISO 19162) |
//! | `projjson` | PROJJSON, plus datum-shift preservation |
//!
//! The cross-format round trip stays here rather than in any one file, because
//! it emits *every* registry CRS to PROJ-string and WKT2 and re-parses both: it
//! is the test that would catch a serializer and its parser drifting together
//! into a private dialect that round-trips but matches nobody.

use super::*;

mod proj_string;
mod projjson;
mod wkt;

/// Compare a parsed CRS to the registry CRS by projecting a test point on each
/// CRS's *own* datum (isolating projection + ellipsoid), returning the worst
/// component difference. Panics if the two CRSes are different kinds.
pub(super) fn projection_delta(from_code: Crs, from_text: Crs, lat: f64, lon: f64) -> f64 {
    let geo = |c: Crs| Crs::Geographic {
        datum: c.datum(),
        epsg: None,
        dim: GeoDim::D2,
    };
    let p = CrsCoord::geographic(GeoPoint::from_lat_lon(lat, lon));
    let a = transform(geo(from_code), from_code, p);
    let b = transform(geo(from_text), from_text, p);
    match (a, b) {
        (CrsCoord::Projected(a), CrsCoord::Projected(b)) => {
            (a.x - b.x).abs().max((a.y - b.y).abs())
        }
        (CrsCoord::Geocentric(a), CrsCoord::Geocentric(b)) => (a.x - b.x)
            .abs()
            .max((a.y - b.y).abs())
            .max((a.z - b.z).abs()),
        (CrsCoord::Geographic { point: a, .. }, CrsCoord::Geographic { point: b, .. }) => {
            (a.lat - b.lat).abs().max((a.lon - b.lon).abs())
        }
        _ => panic!("parsed CRS is a different kind than the registry CRS"),
    }
}

#[test]
fn serializers_round_trip_through_the_parsers() {
    // Each registry CRS emitted to PROJ-string and WKT2, then re-parsed, must
    // project identically. (Spans geographic 2-D/3-D, geocentric, and every
    // projection family.)
    let cases: &[(u32, f64, f64)] = &[
        (4326, 51.5, -0.1),
        (4979, 51.5, -0.1),
        (4978, 45.0, 45.0),
        (4269, 40.0, -100.0),
        (3857, 51.5, -0.1),
        (32_630, 51.5, -0.1),
        (27_700, 52.0, -1.0),
        (2154, 46.5, 2.5),
        (5070, 40.0, -100.0),
        (3031, -75.0, 30.0),
    ];
    for &(code, lat, lon) in cases {
        let crs = Crs::from_epsg(code).unwrap();
        let ps = crs.to_proj_string();
        let from_ps = Crs::from_proj_string(&ps)
            .unwrap_or_else(|| panic!("EPSG:{code}: re-parse of proj-string failed: {ps}"));
        assert!(
            projection_delta(crs, from_ps, lat, lon) < 1e-6,
            "EPSG:{code}: proj-string round-trip off — {ps}"
        );
        let w = crs.to_wkt2();
        let from_w = Crs::from_wkt2(&w)
            .unwrap_or_else(|| panic!("EPSG:{code}: re-parse of WKT2 failed: {w}"));
        assert!(
            projection_delta(crs, from_w, lat, lon) < 1e-6,
            "EPSG:{code}: WKT2 round-trip off — {w}"
        );
        // WKT2 carries the ID, so the code survives the round-trip.
        assert_eq!(
            from_w.epsg(),
            Some(code),
            "EPSG:{code}: WKT2 keeps the code"
        );
    }
}
