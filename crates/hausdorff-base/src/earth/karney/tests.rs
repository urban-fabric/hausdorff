use std::f64::consts::PI;

use crate::earth::ellipsoid::Ellipsoid;

use super::math::sincosd;
use super::solver::{Karney, wgs84};

#[test]
fn sincosd_is_exact_at_quadrants() {
    assert_eq!(sincosd(0.0), (0.0, 1.0));
    assert_eq!(sincosd(90.0), (1.0, 0.0));
    assert_eq!(sincosd(180.0), (0.0, -1.0));
    assert_eq!(sincosd(-90.0), (-1.0, 0.0));
    assert_eq!(sincosd(720.0), (0.0, 1.0));
}

#[test]
fn equator_arc_is_exact() {
    // Along the equator the geodesic is the equator itself: one degree of
    // longitude is exactly a/√(1) scaled by the radian measure.
    let k = Karney::new(Ellipsoid::WGS84);
    let (s12, azi1, azi2) = k.inverse(0.0, 0.0, 0.0, 1.0);
    assert!((azi1 - 90.0).abs() < 1e-13 && (azi2 - 90.0).abs() < 1e-13);
    let expected = Ellipsoid::WGS84.a * 1f64.to_radians();
    assert!((s12 - expected).abs() < 1e-6, "{s12} vs {expected}");
}

#[test]
fn sphere_matches_haversine() {
    // With f = 0 the Karney machinery degenerates to great circles.
    let k = Karney::new(Ellipsoid::SPHERE);
    let a = crate::earth::GeoPoint::from_lat_lon(48.85, 2.35);
    let b = crate::earth::GeoPoint::from_lat_lon(-33.87, 151.21);
    let (s12, _, _) = k.inverse(a.lat, a.lon, b.lat, b.lon);
    let hav = crate::earth::geodesic::haversine_on_m(a, b, Ellipsoid::SPHERE.a);
    assert!((s12 - hav).abs() < 1e-6, "{s12} vs {hav}");
}

#[test]
fn direct_inverts_inverse() {
    let k = wgs84();
    let (lat1, lon1, lat2, lon2) = (40.6398, -73.7787, 51.4775, -0.4614); // JFK → LHR
    let (s12, azi1, _) = k.inverse(lat1, lon1, lat2, lon2);
    let (dlat, dlon, _) = k.direct(lat1, lon1, azi1, s12);
    assert!((dlat - lat2).abs() < 1e-12, "lat {dlat}");
    assert!((dlon - lon2).abs() < 1e-12, "lon {dlon}");
}

#[test]
fn ring_area_matches_reference_value() {
    // The documented geographiclib-rs PolygonArea example: a 1°×1° square
    // on the equator, counter-clockwise.
    let k = wgs84();
    let ring = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)];
    let (perimeter, area) = k.ring_perimeter_area(&ring, true, true);
    assert!(
        (perimeter - 443_770.917_248_302).abs() < 1e-6,
        "{perimeter}"
    );
    assert!((area - 12_308_778_361.469_452).abs() < 1e-3, "{area}");

    // Reversed winding: same magnitude, negative sign when signed.
    let rev = [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)];
    let (_, neg) = k.ring_perimeter_area(&rev, true, true);
    assert!((neg + area).abs() < 1e-3, "{neg}");

    // Unsigned: the reversed ring is read as enclosing the rest of the
    // planet (the large-polygon example: 510 trillion m²).
    let (_, large) = k.ring_perimeter_area(&rev, true, false);
    assert!((large - 510_053_312_945_726.94).abs() < 1.0, "{large}");
}

#[test]
fn nearly_antipodal_converges() {
    // The classic Vincenty-failure regime; Karney must converge and the
    // distance must be slightly under half the ellipsoid circumference.
    let k = wgs84();
    let (s12, azi1, azi2) = k.inverse(0.5, 0.0, -0.4999, 179.9999);
    assert!(s12.is_finite() && azi1.is_finite() && azi2.is_finite());
    let half_meridian = PI * Ellipsoid::WGS84.b(); // lower bound scale
    assert!(s12 > half_meridian, "{s12}");
    assert!(s12 < PI * Ellipsoid::WGS84.a, "{s12}");
}
