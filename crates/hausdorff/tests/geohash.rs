//! Native geohash codec: round-trips and conformance against golden reference
//! strings. The reference values below were captured from the `geohash` crate
//! (and cross-validated by the native encoder at capture time); the live
//! dependency is no longer needed.

use hausdorff::index::geohash::{decode, encode, try_encode};

#[test]
fn matches_golden_reference() {
    let samples = [
        (13.405, 52.52),
        (-0.1278, 51.5074),
        (151.2093, -33.8688),
        (0.0, 0.0),
        (-122.4194, 37.7749),
    ];
    // Reference geohashes at precisions [1, 5, 9, 12] (geohash-crate golden).
    const GOLDEN: [[&str; 4]; 5] = [
        ["u", "u33dc", "u33dc0cpp", "u33dc0cppjs7"],
        ["g", "gcpvj", "gcpvj0duq", "gcpvj0duq533"],
        ["r", "r3gx2", "r3gx2f77b", "r3gx2f77bn44"],
        ["s", "s0000", "s00000000", "s00000000000"],
        ["9", "9q8yy", "9q8yyk8yt", "9q8yyk8ytpxr"],
    ];
    for (&(lon, lat), golden) in samples.iter().zip(GOLDEN.iter()) {
        for (&p, &want) in [1usize, 5, 9, 12].iter().zip(golden.iter()) {
            let got = encode(lon, lat, p);
            assert_eq!(got, want, "geohash mismatch at ({lon},{lat}) prec {p}");
        }
    }
}

#[test]
fn decode_returns_cell_center_near_point() {
    let (lon, lat) = (13.405, 52.52);
    let h = encode(lon, lat, 9);
    let (dlon, dlat) = decode(&h).unwrap();
    assert!((dlon - lon).abs() < 1e-3 && (dlat - lat).abs() < 1e-3);
    assert!(decode("a!b").is_none());
    assert!(decode("").is_none());
}

#[test]
fn try_encode_validates_input() {
    assert_eq!(try_encode(13.405, 52.52, 9), Some(encode(13.405, 52.52, 9)));
    assert!(try_encode(180.001, 0.0, 9).is_none());
    assert!(try_encode(0.0, 90.001, 9).is_none());
    assert!(try_encode(f64::NAN, 0.0, 9).is_none());
    assert!(try_encode(0.0, f64::INFINITY, 9).is_none());
    assert!(try_encode(0.0, 0.0, 0).is_none());
    assert!(try_encode(-180.0, -90.0, 1).is_some());
}
