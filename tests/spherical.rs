//! True spherical geometry: robust triangle area / orientation on unit vectors
//! (`s2` measures + predicates), validated against analytic values — no
//! projection, no oracle needed.

use std::f64::consts::PI;

use hausdorff::earth::{EARTH_RADIUS_M, GeoPoint};
use hausdorff::space::s2::{S2Point, girard_area, point_area, robust_sign, signed_area};

fn p(lat: f64, lon: f64) -> S2Point {
    S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
}

#[test]
fn octant_area_is_half_pi() {
    // (1,0,0), (0,1,0), (0,0,1) span one octant of the sphere = π/2 steradians.
    let (a, b, c) = (p(0.0, 0.0), p(0.0, 90.0), p(90.0, 0.0));
    assert!(
        (point_area(a, b, c) - PI / 2.0).abs() < 1e-12,
        "l'Huilier: {}",
        point_area(a, b, c)
    );
    assert!(
        (girard_area(a, b, c) - PI / 2.0).abs() < 1e-12,
        "Girard: {}",
        girard_area(a, b, c)
    );
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
    assert!(
        rel < 1e-4,
        "spherical {spherical_m2} vs planar {planar_m2} (rel {rel})"
    );
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
    use hausdorff::space::s1::{Angle, AngleDegrees};
    let a = Angle::from_degrees(90.0);
    assert!((a.radians() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    assert!((a.degrees() - 90.0).abs() < 1e-12);
    assert_eq!(a.to_degrees(), AngleDegrees::new(90.0));
    assert!((AngleDegrees::new(180.0).to_radians().radians() - PI).abs() < 1e-12);
}

#[test]
fn angle_constants_and_ops() {
    use hausdorff::space::s1::Angle;
    assert_eq!(Angle::PI.radians(), PI);
    assert_eq!(
        (Angle::from_radians(1.0) + Angle::from_radians(2.0)).radians(),
        3.0
    );
    assert_eq!((Angle::PI * 0.5).radians(), std::f64::consts::FRAC_PI_2);
    assert_eq!(Angle::from_radians(-1.5).abs().radians(), 1.5);
}

// ---- s2::Loop (spherical polygon, unit sphere) -----------------------------

#[test]
fn loop_area_octant_and_hemisphere() {
    use hausdorff::space::s2::Loop;
    let octant = Loop::new(vec![p(0.0, 0.0), p(0.0, 90.0), p(90.0, 0.0)]);
    assert!(
        (octant.area() - PI / 2.0).abs() < 1e-9,
        "octant area {}",
        octant.area()
    );
    // four equator points, CCW as seen from the north pole = northern hemisphere
    let hemi = Loop::new(vec![
        p(0.0, 0.0),
        p(0.0, 90.0),
        p(0.0, 180.0),
        p(0.0, -90.0),
    ]);
    assert!(
        (hemi.area() - 2.0 * PI).abs() < 1e-9,
        "hemisphere area {}",
        hemi.area()
    );
}

#[test]
fn loop_contains_and_perimeter() {
    use hausdorff::space::s2::Loop;
    let hemi = Loop::new(vec![
        p(0.0, 0.0),
        p(0.0, 90.0),
        p(0.0, 180.0),
        p(0.0, -90.0),
    ]);
    assert!(
        hemi.contains_point(p(89.0, 12.0)),
        "north is inside northern hemisphere"
    );
    assert!(!hemi.contains_point(p(-89.0, 12.0)), "south is outside");
    // perimeter is the full equator great circle = 2π
    assert!((hemi.perimeter().radians() - 2.0 * PI).abs() < 1e-9);
}

// ---- planet-flexible metric (Sphere) ----------------------------------------

#[test]
fn sphere_is_planet_flexible() {
    use hausdorff::earth::sphere::Sphere;
    use hausdorff::space::s1::Angle;

    // arc length of one radian equals the radius
    assert!(
        (Sphere::EARTH.length(Angle::from_radians(1.0)).as_meters() - EARTH_RADIUS_M).abs() < 1e-6
    );
    // a full sphere is 4π steradians
    let full = Sphere::EARTH.area_m2(4.0 * PI);
    assert!((full - 4.0 * PI * EARTH_RADIUS_M * EARTH_RADIUS_M).abs() / full < 1e-12);
    // the typed form (euclid AreaSquareMeters, via the LengthMeters × LengthMeters
    // algebraic connection) carries the unit and agrees with the raw m² kernel
    assert_eq!(Sphere::EARTH.area(4.0 * PI).as_square_meters(), full);
    // the same angle is a shorter distance on smaller Mars than on Earth
    let a = Angle::from_degrees(1.0);
    assert!(Sphere::MARS.radius().as_meters() < Sphere::EARTH.radius().as_meters());
    assert!(Sphere::MARS.length(a).as_meters() < Sphere::EARTH.length(a).as_meters());
}

// ---- s2::Cap (spherical disc) ----------------------------------------------

#[test]
fn cap_hemisphere() {
    use hausdorff::space::s1::Angle;
    use hausdorff::space::s2::Cap;
    // a cap of angular radius π/2 about the north pole = the northern hemisphere
    let north = p(90.0, 0.0);
    let cap = Cap::from_center_angle(north, Angle::from_radians(PI / 2.0));
    assert!(
        (cap.area() - 2.0 * PI).abs() < 1e-12,
        "hemisphere cap area {}",
        cap.area()
    );
    assert!(cap.contains_point(north));
    assert!(cap.contains_point(p(45.0, 30.0))); // mid-northern latitude
    assert!(!cap.contains_point(p(-45.0, 30.0))); // southern
}

#[test]
fn cap_empty_full_and_complement() {
    use hausdorff::space::s1::Angle;
    use hausdorff::space::s2::Cap;
    let c = p(10.0, 20.0);
    assert!(Cap::empty(c).is_empty() && Cap::empty(c).area() == 0.0);
    assert!(Cap::full(c).is_full() && (Cap::full(c).area() - 4.0 * PI).abs() < 1e-12);
    // a point outside a cap is inside its complement
    let cap = Cap::from_center_angle(p(90.0, 0.0), Angle::from_radians(PI / 4.0));
    let outside = p(0.0, 0.0);
    assert!(!cap.contains_point(outside));
    assert!(cap.complement().contains_point(outside));
}

// ---- s2::Polygon (shell + holes) -------------------------------------------

#[test]
fn polygon_area_subtracts_holes_and_contains() {
    use hausdorff::space::s2::{Loop, Polygon};
    // big shell triangle (an octant), with a small triangular hole inside it
    let shell = Loop::new(vec![p(1.0, 1.0), p(1.0, 80.0), p(80.0, 1.0)]);
    let hole = Loop::new(vec![p(20.0, 20.0), p(20.0, 30.0), p(30.0, 20.0)]);
    let poly = Polygon::new(shell.clone(), vec![hole.clone()]);

    assert!((poly.area() - (shell.area() - hole.area())).abs() < 1e-12);
    assert!(poly.area() < shell.area());

    // a point in the shell but outside the hole is contained; a point in the
    // hole is not; a point outside the shell is not.
    assert!(poly.contains_point(p(10.0, 10.0)));
    assert!(!poly.contains_point(p(24.0, 24.0))); // inside the hole
    assert!(!poly.contains_point(p(-10.0, -10.0))); // outside the shell
}

// ---- s2::Region (bounding caps) --------------------------------------------

#[test]
fn region_bounding_cap_contains_geometry() {
    use hausdorff::space::s1::Angle;
    use hausdorff::space::s2::{Cap, Loop, Region};
    let lp = Loop::new(vec![p(0.0, 0.0), p(0.0, 30.0), p(20.0, 15.0)]);
    let cap = lp.bounding_cap();
    for v in lp.vertices() {
        assert!(
            cap.contains_point(*v),
            "bounding cap must contain every vertex"
        );
    }
    // a Cap bounds itself
    let c = Cap::from_center_angle(p(45.0, 45.0), Angle::from_radians(0.3));
    assert_eq!(c.bounding_cap(), c);
}

// ---- s2 edge crossings ------------------------------------------------------

#[test]
fn edge_crossing_and_intersection() {
    use hausdorff::space::s2::{Crossing, crossing_sign, intersection};
    // AB: the equator from lon 0 to lon 90.  CD: a meridian arc at lon 45.
    let (a, b) = (p(0.0, 0.0), p(0.0, 90.0));
    let (c, d) = (p(-10.0, 45.0), p(10.0, 45.0));
    assert_eq!(crossing_sign(a, b, c, d), Crossing::Cross);

    let x = intersection(a, b, c, d).to_geo_point();
    assert!(
        (x.lat).abs() < 1e-6 && (x.lon - 45.0).abs() < 1e-6,
        "crossed at ({}, {})",
        x.lat,
        x.lon
    );

    // CD entirely north of AB -> no crossing
    let (c2, d2) = (p(10.0, 10.0), p(10.0, 80.0));
    assert_eq!(crossing_sign(a, b, c2, d2), Crossing::DoNotCross);

    // shared vertex -> Maybe (boundary case)
    assert_eq!(crossing_sign(a, b, a, d), Crossing::Maybe);
}

#[test]
fn loop_intersects() {
    use hausdorff::space::s2::Loop;
    let a = Loop::new(vec![p(0.0, 0.0), p(0.0, 20.0), p(20.0, 10.0)]);
    let overlapping = Loop::new(vec![p(5.0, 5.0), p(5.0, 25.0), p(25.0, 15.0)]);
    let disjoint = Loop::new(vec![p(0.0, 60.0), p(0.0, 80.0), p(20.0, 70.0)]);
    let nested = Loop::new(vec![p(2.0, 5.0), p(2.0, 10.0), p(6.0, 7.0)]); // inside a

    assert!(a.intersects(&overlapping));
    assert!(!a.intersects(&disjoint));
    assert!(a.intersects(&nested));
    assert!(nested.intersects(&a)); // symmetric
}
