//! True spherical geometry: robust triangle area / orientation on unit vectors
//! (`s2` measures + predicates), validated against analytic values — no
//! projection, no oracle needed.

use std::f64::consts::PI;

use hausdorff::geo::{EARTH_RADIUS_M, GeoPoint};
use hausdorff::s2::{S2Point, girard_area, point_area, robust_sign, signed_area};

fn p(lat: f64, lon: f64) -> S2Point {
    S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
}

#[test]
fn octant_area_is_half_pi() {
    // (1,0,0), (0,1,0), (0,0,1) span one octant of the sphere = π/2 steradians.
    let (a, b, c) = (p(0.0, 0.0), p(0.0, 90.0), p(90.0, 0.0));
    assert!((point_area(a, b, c) - PI / 2.0).abs() < 1e-12, "l'Huilier: {}", point_area(a, b, c));
    assert!((girard_area(a, b, c) - PI / 2.0).abs() < 1e-12, "Girard: {}", girard_area(a, b, c));
}

#[test]
fn signed_area_and_sign_follow_orientation() {
    let (a, b, c) = (p(0.0, 0.0), p(0.0, 90.0), p(90.0, 0.0));
    assert_eq!(robust_sign(a, b, c), 1);
    assert_eq!(robust_sign(a, c, b), -1);
    assert!((signed_area(a, b, c) - PI / 2.0).abs() < 1e-12);
    assert!((signed_area(a, c, b) + PI / 2.0).abs() < 1e-12);
}

#[test]
fn small_triangle_matches_planar_area() {
    // At this scale a spherical triangle's area (× R²) equals its planar lon/lat
    // area in m² to within rounding — the l'Huilier path stays accurate for tiny
    // triangles where Girard's formula would lose precision.
    let (a, b, c) = (p(0.0, 0.0), p(0.0, 0.001), p(0.001, 0.0));
    let spherical_m2 = point_area(a, b, c) * EARTH_RADIUS_M * EARTH_RADIUS_M;
    let leg = (0.001_f64).to_radians() * EARTH_RADIUS_M;
    let planar_m2 = 0.5 * leg * leg;
    let rel = (spherical_m2 - planar_m2).abs() / planar_m2;
    assert!(rel < 1e-4, "spherical {spherical_m2} vs planar {planar_m2} (rel {rel})");
}

#[test]
fn degenerate_triangle_has_zero_signed_area() {
    let a = p(10.0, 20.0);
    assert_eq!(robust_sign(a, a, a), 0);
    assert_eq!(signed_area(a, a, a), 0.0);
}

// ---- s1::Angle (the angle type, brought in from euclid) ---------------------

#[test]
fn angle_radians_degrees_round_trip() {
    use hausdorff::s1::{Angle, AngleDegrees};
    let a = Angle::from_degrees(90.0);
    assert!((a.radians() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    assert!((a.degrees() - 90.0).abs() < 1e-12);
    assert_eq!(a.to_degrees(), AngleDegrees::new(90.0));
    assert!((AngleDegrees::new(180.0).to_radians().radians() - PI).abs() < 1e-12);
}

#[test]
fn angle_constants_and_ops() {
    use hausdorff::s1::Angle;
    assert_eq!(Angle::PI.radians(), PI);
    assert_eq!((Angle::from_radians(1.0) + Angle::from_radians(2.0)).radians(), 3.0);
    assert_eq!((Angle::PI * 0.5).radians(), std::f64::consts::FRAC_PI_2);
    assert_eq!(Angle::from_radians(-1.5).abs().radians(), 1.5);
}
