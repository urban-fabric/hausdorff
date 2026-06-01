//! Map-projection interchange: WGS84 geographic ↔ projected planar (Web
//! Mercator, UTM). Validated by analytic anchors, round-trips, and monotonicity.

use std::f64::consts::PI;

use hausdorff::earth::GeoPoint;
use hausdorff::projection::{Ellipsoid, Projection, TransverseMercator, Utm, WebMercator};

fn g(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

#[test]
fn web_mercator_anchors_and_round_trip() {
    let wm = WebMercator::EPSG_3857;
    // origin maps to (0, 0)
    let o = wm.project(g(0.0, 0.0));
    assert!(o.x.abs() < 1e-6 && o.y.abs() < 1e-6);
    // the antimeridian is half the equatorial circumference east
    let e = wm.project(g(0.0, 180.0));
    assert!((e.x - wm.radius * PI).abs() < 1e-3, "x = {}", e.x);
    // round trip
    for &(lat, lon) in &[(0.0, 0.0), (51.5, -0.12), (-33.87, 151.21), (60.0, 24.9)] {
        let back = wm.inverse(wm.project(g(lat, lon)));
        assert!((back.lat - lat).abs() < 1e-9 && (back.lon - lon).abs() < 1e-9);
    }
}

#[test]
fn utm_central_meridian_anchor() {
    // On a zone's central meridian at the equator: easting = false easting
    // (500 km), northing = 0, in the expected zone.
    let u = Utm::from_geo(g(0.0, 3.0)); // zone 31, CM 3°E
    assert_eq!(u.zone, 31);
    assert!(u.north);
    assert!(
        (u.easting - 500_000.0).abs() < 1e-3,
        "easting {}",
        u.easting
    );
    assert!(u.northing.abs() < 1e-3, "northing {}", u.northing);
}

#[test]
fn utm_round_trips_to_mm() {
    for &(lat, lon) in &[
        (40.7128, -74.0060),  // New York (zone 18N)
        (48.8566, 2.3522),    // Paris (zone 31N)
        (-33.8688, 151.2093), // Sydney (zone 56S)
        (35.6762, 139.6503),  // Tokyo (zone 54N)
    ] {
        let back = Utm::from_geo(g(lat, lon)).to_geo();
        assert!((back.lat - lat).abs() < 1e-8, "lat {lat} -> {}", back.lat);
        assert!((back.lon - lon).abs() < 1e-8, "lon {lon} -> {}", back.lon);
    }
}

#[test]
fn utm_easting_increases_eastward() {
    let cm = Utm::from_geo(g(45.0, 3.0)).easting; // on the central meridian
    let east = Utm::from_geo(g(45.0, 5.0)).easting; // east of it
    assert!(east > cm && (cm - 500_000.0).abs() < 1.0);
}

#[test]
fn projection_is_ellipsoid_flexible() {
    // The same TransverseMercator math runs on any ellipsoid — a different datum
    // (or planet) yields different planar coordinates for the same lat/lon.
    let wgs84 = TransverseMercator {
        ellipsoid: Ellipsoid::WGS84,
        central_meridian: 0.0,
        k0: 0.9996,
        false_easting: 500_000.0,
        false_northing: 0.0,
    };
    let sphere = TransverseMercator {
        ellipsoid: Ellipsoid::SPHERE,
        ..wgs84
    };
    let p = g(45.0, 2.0);
    assert!((wgs84.project(p).y - sphere.project(p).y).abs() > 1.0);
    // a spherical-ellipsoid projection still round-trips
    let back = sphere.inverse(sphere.project(p));
    assert!((back.lat - 45.0).abs() < 1e-8 && (back.lon - 2.0).abs() < 1e-8);
}
