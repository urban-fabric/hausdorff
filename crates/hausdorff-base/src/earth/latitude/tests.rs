use super::*;

fn close(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

/// On a sphere (`f = 0`) the three latitudes coincide and both radii are `a`.
#[test]
fn sphere_collapses_the_latitude_system() {
    let s = Ellipsoid::SPHERE;
    for &lat in &[0.0, 12.3, 45.0, -60.0, 89.9] {
        assert!(close(s.reduced_latitude(lat), lat, 1e-12));
        assert!(close(s.geocentric_latitude(lat), lat, 1e-12));
        assert!(close(s.prime_vertical_radius(lat), s.a, 1e-6));
        assert!(close(s.meridian_radius(lat), s.a, 1e-6));
    }
}

/// Each conversion inverts cleanly (WGS84).
#[test]
fn conversions_round_trip() {
    let e = Ellipsoid::WGS84;
    for &lat in &[0.0, 12.3, 45.0, -60.0, 80.0] {
        assert!(close(
            e.geodetic_from_reduced(e.reduced_latitude(lat)),
            lat,
            1e-9
        ));
        assert!(close(
            e.geodetic_from_geocentric(e.geocentric_latitude(lat)),
            lat,
            1e-9
        ));
    }
}

/// For `0 < φ < 90` the spheroid squashes the latitudes toward the equator,
/// `|ψ| < |β| < |φ|`, with the standard WGS84 values at 45°.
#[test]
fn latitudes_ordered_and_known() {
    let e = Ellipsoid::WGS84;
    let beta = e.reduced_latitude(45.0);
    let psi = e.geocentric_latitude(45.0);
    assert!(psi < beta && beta < 45.0);
    assert!(close(beta, 44.9038, 1e-2));
    assert!(close(psi, 44.8076, 1e-2));
}

/// The geodetic normal is a unit vector whose elevation (z-component) is
/// `sin φ` — i.e. its angle above the equatorial plane is the geodetic
/// latitude, the defining property.
#[test]
fn normal_is_unit_with_latitude_elevation() {
    for &(lat, lon) in &[(0.0, 0.0), (45.0, 30.0), (-60.0, 170.0), (89.0, -120.0)] {
        let n = geodetic_normal(lat, lon);
        let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(close(norm, 1.0, 1e-12));
        assert!(close(n[2], lat.to_radians().sin(), 1e-12));
    }
}

/// `N ≥ a` (equality at the equator); `M` brackets `a` (smaller at the
/// equator, larger at the pole).
#[test]
fn radii_bracket_the_semi_major_axis() {
    let e = Ellipsoid::WGS84;
    assert!(close(e.prime_vertical_radius(0.0), e.a, 1e-6));
    assert!(e.prime_vertical_radius(45.0) > e.a);
    assert!(e.meridian_radius(0.0) < e.a);
    assert!(e.meridian_radius(90.0) > e.a);
}

/// The generic [`SurfaceOfRevolution::total_area`] over the ellipsoid (the
/// trait's numeric integrator) reproduces the closed-form `total_area_m2`, and
/// the sphere instance reproduces `4πR²`.
#[test]
fn surface_of_revolution_total_area() {
    use crate::space::surface::SurfaceOfRevolution;
    let e = Ellipsoid::WGS84;
    let rel = (e.total_area(200_000) - e.total_area_m2()).abs() / e.total_area_m2();
    assert!(rel < 1e-7, "WGS84 SoR area rel_err = {rel}");
    // F = 0 (orthogonal net); the trait's generic area element (√(EG)) agrees
    // with the inherent closed form M·N·cosφ.
    let (_, f, _) = e.first_fundamental_form(45.0);
    assert_eq!(f, 0.0);
    let generic = SurfaceOfRevolution::area_element(&e, 45.0);
    assert!((generic - e.area_element(45.0)).abs() / generic < 1e-12);
    // Sphere: total area 4πR².
    let s = Ellipsoid::SPHERE;
    let exact = 4.0 * core::f64::consts::PI * s.a * s.a;
    assert!((s.total_area(200_000) - exact).abs() / exact < 1e-7);
}

/// The local area element `M·N·cosφ`, integrated over the whole spheroid,
/// reproduces the closed-form `total_area_m2` — two independent area
/// computations (the first-fundamental-form integrand vs. the authalic
/// radius) agree.
#[test]
fn area_element_integrates_to_total_area() {
    use core::f64::consts::{FRAC_PI_2, PI};
    let e = Ellipsoid::WGS84;
    let n = 100_000;
    let dphi = PI / n as f64; // φ from -π/2 to π/2 (radians)
    let mut sum = 0.0;
    for i in 0..n {
        let phi = -FRAC_PI_2 + (i as f64 + 0.5) * dphi; // midpoint
        sum += e.area_element(phi.to_degrees()) * dphi;
    }
    let area = 2.0 * PI * sum; // ∫ dλ = 2π
    let rel_err = (area - e.total_area_m2()).abs() / e.total_area_m2();
    assert!(rel_err < 1e-7, "rel_err = {rel_err}");
}

/// On a sphere (`f = 0`) the conformal / authalic / rectifying latitudes also
/// collapse to `φ` (the `e → 0` limit).
#[test]
fn sphere_collapses_auxiliary_latitudes() {
    let s = Ellipsoid::SPHERE;
    for &lat in &[0.0, 12.3, 45.0, -60.0, 89.9] {
        assert!(close(s.conformal_latitude(lat), lat, 1e-12));
        assert!(close(s.authalic_latitude(lat), lat, 1e-12));
        assert!(close(s.rectifying_latitude(lat), lat, 1e-12));
    }
}

/// Each auxiliary latitude inverts cleanly (WGS84), both hemispheres. The
/// conformal/authalic inverses iterate to convergence (≤ 1e-9°); the
/// rectifying inverse is Snyder's footpoint *series* (the TM's), exact only to
/// the series order (~sub-mm, ≤ 1e-7°).
#[test]
fn auxiliary_latitudes_round_trip() {
    let e = Ellipsoid::WGS84;
    for &lat in &[0.0, 12.3, 45.0, -60.0, 80.0, -88.0] {
        assert!(close(
            e.geodetic_from_conformal(e.conformal_latitude(lat)),
            lat,
            1e-9
        ));
        assert!(close(
            e.geodetic_from_authalic(e.authalic_latitude(lat)),
            lat,
            1e-9
        ));
        assert!(close(
            e.geodetic_from_rectifying(e.rectifying_latitude(lat)),
            lat,
            1e-7
        ));
    }
}

/// WGS84 auxiliary latitudes vs an **independent numerical integration** of the
/// area / meridian-arc elements (2e6-step midpoint — not the Snyder series the
/// methods use). At 45° the conformal squashes most, the authalic least, so
/// `χ < μ < ξ < φ`.
#[test]
fn auxiliary_latitudes_match_independent_integration() {
    let e = Ellipsoid::WGS84;
    for &(lat, chi, xi, mu) in &[
        (30.0, 29.833_682, 29.888_997, 29.875_148),
        (45.0, 44.807_684, 44.871_703, 44.855_682),
        (60.0, 59.833_216, 59.888_786, 59.874_886),
    ] {
        assert!(close(e.conformal_latitude(lat), chi, 1e-4), "χ({lat})");
        assert!(close(e.authalic_latitude(lat), xi, 1e-4), "ξ({lat})");
        assert!(close(e.rectifying_latitude(lat), mu, 1e-4), "μ({lat})");
    }
    // The squash ordering at 45°: χ < μ < ξ < φ.
    assert!(e.conformal_latitude(45.0) < e.rectifying_latitude(45.0));
    assert!(e.rectifying_latitude(45.0) < e.authalic_latitude(45.0));
    assert!(e.authalic_latitude(45.0) < 45.0);
    // Quarter-meridian (equator → pole), WGS84 ≈ 10 001 965.73 m.
    assert!(close(e.meridian_arc_m(90.0), 10_001_965.73, 0.1));
    assert!(close(e.meridian_arc_m(-90.0), -10_001_965.73, 0.1));
}
