//! Native geohash codec: round-trips and conformance against the `geohash` crate
//! reference.

use hausdorff::geohash::{decode, encode};

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
}
