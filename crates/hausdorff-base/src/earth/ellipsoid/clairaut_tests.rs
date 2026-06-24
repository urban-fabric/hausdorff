use super::*;

#[test]
fn constant_is_invariant_along_a_geodesic() {
    // Clairaut's relation: ρ·sin α is conserved along the whole geodesic.
    // Verify it numerically against the Karney solver — the runtime invariant.
    let e = Ellipsoid::WGS84;
    let start = GeoPoint::from_lat_lon(40.0, -75.0);
    let az = Bearing::from_degrees(63.0);
    let line = e.line(start, az);
    let c0 = e.clairaut_constant(start, az);
    for km in (500..=9000).step_by(500) {
        let gd = line.position_m(f64::from(km) * 1000.0);
        let c = e.clairaut_constant(gd.point, gd.final_bearing);
        // ρ ~ 6.4e6 m; the constant holds to well under a millimetre.
        assert!((c - c0).abs() < 1e-3, "at {km} km: clairaut {c} vs {c0}");
    }
}

#[test]
fn auxiliary_sphere_form_matches() {
    // cos β · sin α == (ρ sin α) / a — the auxiliary-sphere reduction.
    let e = Ellipsoid::WGS84;
    let p = GeoPoint::from_lat_lon(35.0, 12.0);
    let az = Bearing::from_degrees(48.0);
    assert!((e.clairaut_sin_alpha0(p, az) - e.clairaut_constant(p, az) / e.a).abs() < 1e-12);
    // On a sphere β = φ, so the two reductions coincide there.
    let s = Ellipsoid::SPHERE;
    assert!(
        (s.clairaut_sin_alpha0(p, az) - p.lat.to_radians().cos() * az.radians().sin()).abs()
            < 1e-12
    );
}

#[test]
fn meridian_and_equator_are_degenerate_geodesics() {
    let e = Ellipsoid::WGS84;
    // Due north (α = 0): the meridian is a geodesic with Clairaut constant 0.
    let north = e.clairaut_constant(
        GeoPoint::from_lat_lon(20.0, 5.0),
        Bearing::from_degrees(0.0),
    );
    assert!(north.abs() < 1e-9);
    // Due east at the equator (α = 90°, φ = 0): the equator is a geodesic,
    // Clairaut constant = a (ρ = N = a there), and it stays a downrange.
    let eq = GeoPoint::from_lat_lon(0.0, 0.0);
    let east = Bearing::from_degrees(90.0);
    assert!((e.clairaut_constant(eq, east) - e.a).abs() < 1e-6);
    let gd = e.line(eq, east).position_m(3_000_000.0);
    assert!((e.clairaut_constant(gd.point, gd.final_bearing) - e.a).abs() < 1e-3);
}

#[test]
fn vertex_latitude_matches_clairaut_and_sampling() {
    let e = Ellipsoid::WGS84;
    let start = GeoPoint::from_lat_lon(0.0, 0.0);
    let az = Bearing::from_degrees(45.0);
    let phi_max = e.geodesic_vertex_latitude(start, az);
    // Geodetic exceeds reduced, so a 45° equatorial launch tops out just above 45°.
    assert!(phi_max > 45.0 && phi_max < 45.3, "phi_max = {phi_max}");
    // Exact Clairaut check: ρ·sin α at the vertex (φ_max, due east) equals the
    // launch constant.
    let c0 = e.clairaut_constant(start, az);
    let c_vertex = e.clairaut_constant(
        GeoPoint::from_lat_lon(phi_max, 0.0),
        Bearing::from_degrees(90.0),
    );
    assert!((c_vertex - c0).abs() < 1e-6, "{c_vertex} vs {c0}");
    // Sampling: the geodesic never exceeds the vertex, and a fine walk reaches it.
    let line = e.line(start, az);
    let mut sampled_max = f64::MIN;
    for km in (0..=12_000).step_by(50) {
        let lat = line.point_m(f64::from(km) * 1000.0).lat;
        sampled_max = sampled_max.max(lat);
        assert!(lat <= phi_max + 1e-6, "lat {lat} exceeds vertex {phi_max}");
    }
    assert!(
        (phi_max - sampled_max).abs() < 0.02,
        "sampled {sampled_max} vs {phi_max}"
    );
}

#[test]
fn meridian_reaches_the_pole_equator_stays() {
    let e = Ellipsoid::WGS84;
    // Due north → a meridian → reaches the pole (vertex latitude 90°).
    let pole = e.geodesic_vertex_latitude(
        GeoPoint::from_lat_lon(10.0, 20.0),
        Bearing::from_degrees(0.0),
    );
    assert!((pole - 90.0).abs() < 1e-6, "{pole}");
    // Due east at the equator → the equator → vertex latitude 0.
    let eq = e.geodesic_vertex_latitude(
        GeoPoint::from_lat_lon(0.0, 0.0),
        Bearing::from_degrees(90.0),
    );
    assert!(eq.abs() < 1e-6, "{eq}");
}
