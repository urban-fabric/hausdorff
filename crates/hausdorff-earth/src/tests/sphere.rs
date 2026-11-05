//! The planet-flexible metric layer: `Sphere` turning the radius-free spherical
//! geometry below it into physical lengths and areas.
//!
//! This is the one test in the spherical family that needs a radius to state its
//! expectation, which is why it lives here and not with the `s2` suites in
//! `hausdorff-sphere` — that crate is deliberately radius-free, and `Sphere` is
//! the bridge out of it.

use std::f64::consts::PI;

use crate::earth::EARTH_RADIUS_M;
use crate::earth::sphere::Sphere;
use hausdorff_base::space::s1::Angle;

#[test]
fn sphere_is_planet_flexible() {
    assert!((Sphere::EARTH.length(Angle::from_radians(1.0)).value() - EARTH_RADIUS_M).abs() < 1e-6);
    let full = Sphere::EARTH.area_m2(4.0 * PI);
    assert!((full - 4.0 * PI * EARTH_RADIUS_M * EARTH_RADIUS_M).abs() / full < 1e-12);
    assert_eq!(Sphere::EARTH.area(4.0 * PI).value(), full);
    let a = Angle::from_degrees(1.0);
    assert!(Sphere::MARS.radius().value() < Sphere::EARTH.radius().value());
    assert!(Sphere::MARS.length(a).value() < Sphere::EARTH.length(a).value());
}
