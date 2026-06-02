//! Native geohash codec: round-trips and conformance against the `geohash` crate
//! reference.

use hausdorff::index::geohash::{decode, encode, try_encode};

#[test]
fn matches_geohash_crate() {
    let samples = [
        (13.405, 52.52),   // Berlin
        (-0.1278, 51.5074),
        (151.2093, -33.8688),
        (0.0, 0.0),
        (-122.4194, 37.7749),
    ];
    for &(lon, lat) in &samples {
        for prec in [1usize, 5, 9, 12] {
            let ours = encode(lon, lat, prec);
            let theirs = geohash::encode(geohash::Coord { x: lon, y: lat }, prec).unwrap();
            assert_eq!(ours, theirs, "geohash mismatch at ({lon},{lat}) prec {prec}");
        }
    }
}

#[test]
fn decode_returns_cell_center_near_point() {
    let (lon, lat) = (13.405, 52.52);
    let h = encode(lon, lat, 9);
    let (dlon, dlat) = decode(&h).unwrap();
    // a 9-char geohash is ~2-5 m; center is well within 0.001 deg of the point
    assert!((dlon - lon).abs() < 1e-3 && (dlat - lat).abs() < 1e-3);
    assert!(decode("a!b").is_none()); // invalid char
    assert!(decode("").is_none()); // empty geohash is invalid
}

#[test]
fn try_encode_validates_input() {
    // valid input agrees with the unchecked form
    assert_eq!(try_encode(13.405, 52.52, 9), Some(encode(13.405, 52.52, 9)));
    // out-of-range, non-finite, or zero-precision → None
    assert!(try_encode(180.001, 0.0, 9).is_none());
    assert!(try_encode(0.0, 90.001, 9).is_none());
    assert!(try_encode(f64::NAN, 0.0, 9).is_none());
    assert!(try_encode(0.0, f64::INFINITY, 9).is_none());
    assert!(try_encode(0.0, 0.0, 0).is_none());
    // inclusive bounds are accepted
    assert!(try_encode(-180.0, -90.0, 1).is_some());
}
