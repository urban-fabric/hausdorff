use super::*;
use std::f64::consts::PI;

/// Clarke 1866, the ellipsoid of Snyder's worked numerical examples.
fn clarke1866() -> Ellipsoid {
    Ellipsoid {
        a: 6_378_206.4,
        f: 1.0 / 294.978_698_2,
    }
}

fn approx(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

#[test]
fn lcc_matches_snyder_worked_example() {
    // Snyder PP1395 p. 295: Clarke 1866, standard parallels 33°/45°,
    // origin (23°, −96°), point (35°, −75°) → (1 894 410.9, 1 564 649.5) m.
    let lcc = LambertConformalConic {
        ellipsoid: clarke1866(),
        lat1: 33.0,
        lat2: 45.0,
        lat0: 23.0,
        lon0: -96.0,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    let c = lcc.project(GeoPoint::from_lat_lon(35.0, -75.0));
    assert!(approx(c.x, 1_894_410.9, 0.2), "x = {}", c.x);
    assert!(approx(c.y, 1_564_649.5, 0.2), "y = {}", c.y);
    // On the central meridian, x = false easting (θ = 0).
    let cm = lcc.project(GeoPoint::from_lat_lon(40.0, -96.0));
    assert!(approx(cm.x, 0.0, 1e-6), "cm x = {}", cm.x);
    // Round-trip.
    let back = lcc.inverse(c);
    assert!(approx(back.lat, 35.0, 1e-9) && approx(back.lon, -75.0, 1e-9));
}

#[test]
fn albers_matches_snyder_worked_example() {
    // Snyder PP1395 p. 292: Clarke 1866, standard parallels 29.5°/45.5°,
    // origin (23°, −96°), point (35°, −75°) → (1 885 472.7, 1 535 925.0) m.
    let albers = AlbersEqualArea {
        ellipsoid: clarke1866(),
        lat1: 29.5,
        lat2: 45.5,
        lat0: 23.0,
        lon0: -96.0,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    let c = albers.project(GeoPoint::from_lat_lon(35.0, -75.0));
    assert!(approx(c.x, 1_885_472.7, 0.2), "x = {}", c.x);
    assert!(approx(c.y, 1_535_925.0, 0.2), "y = {}", c.y);
    let back = albers.inverse(c);
    assert!(approx(back.lat, 35.0, 1e-9) && approx(back.lon, -75.0, 1e-9));
}

#[test]
fn conics_round_trip_both_sides_and_hemispheres() {
    // Regression for the inverse θ sign: points west of the central meridian
    // (θ < 0) and a southern-hemisphere conic (cone constant n < 0).
    let north = LambertConformalConic {
        ellipsoid: Ellipsoid::GRS80,
        lat1: 33.0,
        lat2: 45.0,
        lat0: 39.0,
        lon0: -96.0,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    let south = AlbersEqualArea {
        ellipsoid: Ellipsoid::GRS80,
        lat1: -33.0,
        lat2: -45.0,
        lat0: -39.0,
        lon0: 145.0,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    // West and east of the central meridian, both hemispheres.
    for &(lat, lon) in &[(40.0, -110.0), (40.0, -82.0), (38.5, -96.0)] {
        let b = north.inverse(north.project(GeoPoint::from_lat_lon(lat, lon)));
        assert!(
            approx(b.lat, lat, 1e-8) && approx(b.lon, lon, 1e-8),
            "LCC {lat},{lon}"
        );
    }
    for &(lat, lon) in &[(-40.0, 130.0), (-40.0, 160.0), (-38.5, 145.0)] {
        let b = south.inverse(south.project(GeoPoint::from_lat_lon(lat, lon)));
        assert!(
            approx(b.lat, lat, 1e-8) && approx(b.lon, lon, 1e-8),
            "Albers S {lat},{lon}"
        );
    }
}

#[test]
fn polar_stereographic_spherical_limit_and_round_trip() {
    // On a sphere (e = 0) the ellipsoidal formula collapses to the closed
    // form ρ = a·cos(φc)·tan(π/4 − φ/2)/tan(π/4 − φc/2).
    let sphere = PolarStereographic {
        ellipsoid: Ellipsoid::SPHERE,
        lat_ts: 60.0,
        lon0: 0.0,
        north: true,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    let a = Ellipsoid::SPHERE.a;
    let c = sphere.project(GeoPoint::from_lat_lon(80.0, 0.0));
    let phic = 60.0_f64.to_radians();
    let phi = 80.0_f64.to_radians();
    let expect_rho = a * phic.cos() * (PI / 4.0 - phi / 2.0).tan() / (PI / 4.0 - phic / 2.0).tan();
    // At lon = lon0 the point is straight "down" from the pole: x = 0, y = −ρ.
    assert!(approx(c.x, 0.0, 1e-6), "x = {}", c.x);
    assert!(
        approx(c.y, -expect_rho, 1.0),
        "y = {} vs {}",
        c.y,
        -expect_rho
    );

    // Ellipsoidal round-trips (WGS84), both hemispheres.
    let north = PolarStereographic {
        ellipsoid: Ellipsoid::WGS84,
        lat_ts: 71.0,
        lon0: 0.0,
        north: true,
        false_easting: 2_000_000.0,
        false_northing: 2_000_000.0,
    };
    for &(lat, lon) in &[(85.0, 30.0), (72.0, -120.0), (89.0, 175.0)] {
        let back = north.inverse(north.project(GeoPoint::from_lat_lon(lat, lon)));
        assert!(
            approx(back.lat, lat, 1e-7) && approx(back.lon, lon, 1e-6),
            "N {lat},{lon}"
        );
    }
    let south = PolarStereographic {
        ellipsoid: Ellipsoid::WGS84,
        lat_ts: -71.0,
        lon0: 0.0,
        north: false,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    for &(lat, lon) in &[(-85.0, 30.0), (-72.0, -120.0), (-89.0, 175.0)] {
        let back = south.inverse(south.project(GeoPoint::from_lat_lon(lat, lon)));
        assert!(
            approx(back.lat, lat, 1e-7) && approx(back.lon, lon, 1e-6),
            "S {lat},{lon}"
        );
    }
}

#[test]
fn local_tangent_plane_matches_proj_topocentric() {
    // PROJ `+proj=topocentric +lat_0=52 +lon_0=10 +ellps=WGS84`, East/North.
    let ltp = LocalTangentPlane::new(Ellipsoid::WGS84, GeoPoint::from_lat_lon(52.0, 10.0));
    let refs = [
        (52.0, 10.0, 0.0, 0.0),
        (52.01, 10.0, 0.0, 1_112.674_477),
        (52.0, 10.01, 686.780_157, 0.047_228),
        (52.05, 10.07, 4_802.102_024, 5_565.702_286),
    ];
    for &(lat, lon, e, n) in &refs {
        let c = ltp.project(GeoPoint::from_lat_lon(lat, lon));
        assert!(approx(c.x, e, 1e-3), "E {} vs {e}", c.x);
        assert!(approx(c.y, n, 1e-3), "N {} vs {n}", c.y);
        // Round-trip back to the ellipsoid surface (the inverse solves for the
        // surface-landing `up`), recovering the original lat/lon.
        let back = ltp.inverse(c);
        assert!(
            approx(back.lat, lat, 1e-9) && approx(back.lon, lon, 1e-9),
            "rt {lat},{lon}"
        );
    }
    // Planar round-trip project(inverse(c)) = c is exact for any local coord.
    for &(e, n) in &[(0.0, 0.0), (1234.0, -5678.0), (9000.0, 7000.0)] {
        let c = ltp.project(ltp.inverse(Coord::new(e, n)));
        assert!(approx(c.x, e, 1e-6) && approx(c.y, n, 1e-6));
    }
    // The origin maps to the local origin (0, 0).
    let o = ltp.project(GeoPoint::from_lat_lon(52.0, 10.0));
    assert!(o.x.abs() < 1e-9 && o.y.abs() < 1e-9);
}

/// A conformal map sends the (perpendicular) meridian and parallel directions
/// to perpendicular projected directions — a metric-free conformality check.
fn assert_conformal<P: ConformalProjection>(p: &P, lat: f64, lon: f64) {
    let d = 1e-4;
    let g = |la: f64, lo: f64| p.project(GeoPoint::from_lat_lon(la, lo));
    let a = g(lat + d, lon);
    let b = g(lat - d, lon);
    let dm = [a.x - b.x, a.y - b.y]; // projected meridian direction
    let a = g(lat, lon + d);
    let b = g(lat, lon - d);
    let dp = [a.x - b.x, a.y - b.y]; // projected parallel direction
    let cos = (dm[0] * dp[0] + dm[1] * dp[1]) / (dm[0].hypot(dm[1]) * dp[0].hypot(dp[1]));
    assert!(
        cos.abs() < 1e-6,
        "not conformal at {lat},{lon}: cos∠ = {cos}"
    );
}

#[test]
fn utm_wrapper_wires_kernel_parameters() {
    // Locks the 8-argument delegation to the Transverse Mercator kernel: a
    // wrong central meridian / false origin / k0 would break these.
    // On the central meridian at the equator the projection lands on the false
    // origin (500000, 0) for the northern zone.
    let zone31 = TransverseMercator::utm(31, true); // central meridian 3°E
    let o = zone31.project(GeoPoint::from_lat_lon(0.0, 3.0));
    assert!(approx(o.x, 500_000.0, 1e-6), "easting = {}", o.x);
    assert!(approx(o.y, 0.0, 1e-6), "northing = {}", o.y);
    // Southern zone carries the 10 000 000 m false northing at the equator.
    let south = TransverseMercator::utm(31, false);
    let os = south.project(GeoPoint::from_lat_lon(0.0, 3.0));
    assert!(approx(os.y, 10_000_000.0, 1e-6), "S northing = {}", os.y);
    // Round-trip through the `Utm` convenience across hemispheres/zones.
    for &(lat, lon) in &[(48.8566, 2.3522), (-33.8688, 151.2093), (64.0, -147.0)] {
        let back = Utm::from_geo(GeoPoint::from_lat_lon(lat, lon)).to_geo();
        assert!(
            approx(back.lat, lat, 1e-7) && approx(back.lon, lon, 1e-7),
            "rt {lat},{lon}"
        );
    }
}

#[test]
fn conformal_projections_preserve_angles_and_scale() {
    // Web Mercator: the point scale is exactly sec φ (proof-backed), and the
    // map is conformal everywhere.
    let wm = WebMercator::EPSG_3857;
    for &lat in &[0.0, 20.0, 45.0, 70.0] {
        assert!((wm.point_scale(lat) - 1.0 / lat.to_radians().cos()).abs() < 1e-12);
        assert_conformal(&wm, lat, 10.0);
    }

    // Lambert Conformal Conic: scale is unity at both standard parallels and
    // the map is conformal.
    let lcc = LambertConformalConic {
        ellipsoid: Ellipsoid::GRS80,
        lat1: 44.0,
        lat2: 49.0,
        lat0: 46.5,
        lon0: 3.0,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    assert!((lcc.point_scale(44.0) - 1.0).abs() < 1e-9);
    assert!((lcc.point_scale(49.0) - 1.0).abs() < 1e-9);
    assert!(lcc.point_scale(46.5) < 1.0); // between the parallels, k < 1
    assert_conformal(&lcc, 46.5, 3.0);

    // Polar Stereographic: scale is unity at the latitude of true scale.
    let ps = PolarStereographic {
        ellipsoid: Ellipsoid::WGS84,
        lat_ts: 71.0,
        lon0: 0.0,
        north: true,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    assert!((ps.point_scale(71.0) - 1.0).abs() < 1e-9);
    assert_conformal(&ps, 80.0, 30.0);
}

#[test]
fn albers_preserves_area() {
    // Equal-area: a graticule cell projects to the *same* area as it has on
    // the ellipsoid (`Ellipsoid::area_element` = the true area per radian²).
    let albers = AlbersEqualArea {
        ellipsoid: Ellipsoid::GRS80,
        lat1: 29.5,
        lat2: 45.5,
        lat0: 23.0,
        lon0: -96.0,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    assert!((albers.area_scale() - 1.0).abs() < 1e-12);
    let h = 0.01; // half-cell, degrees
    for &(lat, lon) in &[(35.0, -100.0), (45.0, -90.0), (30.0, -110.0)] {
        // Project the four cell corners and take the shoelace area.
        let c = |dla: f64, dlo: f64| albers.project(GeoPoint::from_lat_lon(lat + dla, lon + dlo));
        let p = [c(-h, -h), c(-h, h), c(h, h), c(h, -h)];
        let mut shoelace = 0.0;
        for i in 0..4 {
            let a = p[i];
            let b = p[(i + 1) % 4];
            shoelace += a.x * b.y - b.x * a.y;
        }
        let projected_area = shoelace.abs() / 2.0;
        // True area: area element (m²/rad²) × Δφ × Δλ (radians).
        let cell = (2.0 * h).to_radians();
        let true_area = albers.ellipsoid.area_element(lat) * cell * cell;
        let rel = (projected_area - true_area).abs() / true_area;
        assert!(rel < 1e-4, "area not preserved at {lat},{lon}: rel = {rel}");
    }
}
