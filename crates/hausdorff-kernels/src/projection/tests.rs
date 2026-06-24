use super::*;

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual} != expected {expected} within {tolerance}"
    );
}

const WGS84_A: f64 = 6_378_137.0;
const WGS84_F: f64 = 1.0 / 298.257_223_563;

/// The Transverse Mercator parameters of UTM zone `zone` (north).
fn utm(zone: u8) -> (f64, KrugerSeries) {
    let central_meridian = f64::from(i32::from(zone) * 6 - 183);
    (
        central_meridian,
        KrugerSeries::for_ellipsoid(WGS84_A, WGS84_F),
    )
}

#[test]
fn transverse_mercator_central_meridian_equator_is_false_origin() {
    // On the central meridian at the equator (UTM lat_origin = 0) the
    // projection lands exactly on the false origin: (500000, 0).
    let (cm, ks) = utm(31); // zone 31N → central meridian 3°E
    let (e, n) = transverse_mercator_forward(0.0_f64, cm, cm, 0.0, 0.9996, 500_000.0, 0.0, &ks);
    assert_close(e, 500_000.0, 1e-6);
    assert_close(n, 0.0, 1e-6);
}

#[test]
fn transverse_mercator_round_trips() {
    let (cm, ks) = utm(31);
    for &(lat, lon) in &[(0.0, 3.0), (48.5, 2.0), (-33.9, 4.5), (60.0, 1.0)] {
        let (e, n) = transverse_mercator_forward(lat, lon, cm, 0.0, 0.9996, 500_000.0, 0.0, &ks);
        let (blat, blon) = transverse_mercator_inverse(e, n, cm, 0.0, 0.9996, 500_000.0, 0.0, &ks);
        assert_close(blat, lat, 1e-9);
        assert_close(blon, lon, 1e-9);
    }
}

#[test]
fn web_mercator_round_trips_and_matches_closed_form() {
    let r = WGS84_A;
    // y at the equator is 0; x is R·λ.
    let (x, y) = web_mercator_forward(0.0_f64, 45.0, r);
    assert_close(x, r * 45.0_f64.to_radians(), 1e-6);
    assert_close(y, 0.0, 1e-9);
    for &(lat, lon) in &[(0.0, 45.0), (40.0, -75.0), (-12.5, 130.0)] {
        let (x, y) = web_mercator_forward(lat, lon, r);
        let (blat, blon) = web_mercator_inverse(x, y, r);
        assert_close(blat, lat, 1e-9);
        assert_close(blon, lon, 1e-9);
    }
}

#[test]
fn projections_run_for_f32() {
    // Genericity over `Real`: the same kernels run at `f32` (hence on
    // any lazy backend — `f64`/`f32`/tensor).
    let ks32 = KrugerSeries::for_ellipsoid(WGS84_A, WGS84_F);
    let (e, n) =
        transverse_mercator_forward(48.5_f32, 2.0_f32, 3.0, 0.0, 0.9996, 500_000.0, 0.0, &ks32);
    assert!(e.is_finite() && n.is_finite());
    let (x, y) = web_mercator_forward(40.0_f32, -75.0_f32, WGS84_A as f32);
    assert!(x.is_finite() && y.is_finite());
}

#[test]
fn lat_lon_entry_matches_scalar() {
    let r = WGS84_A;
    let p = (40.0_f64, -75.0_f64); // (lat, lon) implements LatLon
    let viatuple = web_mercator_project_at(&p, r);
    let direct = web_mercator_forward(40.0_f64, -75.0_f64, r);
    assert_close(viatuple.0, direct.0, 1e-12);
    assert_close(viatuple.1, direct.1, 1e-12);
}

// Clarke 1866 — the ellipsoid of Snyder's worked numerical examples.
const CLARKE_A: f64 = 6_378_206.4;
const CLARKE_F: f64 = 1.0 / 294.978_698_2;

fn ecc(f: f64) -> f64 {
    (f * (2.0 - f)).sqrt()
}

#[test]
fn lcc_matches_snyder_and_round_trips() {
    // Snyder PP1395 p. 295: parallels 33°/45°, origin (23°, −96°),
    // point (35°, −75°) → (1 894 410.9, 1 564 649.5) m.
    let c = LccConstants::new(CLARKE_A, ecc(CLARKE_F), 33.0, 45.0, 23.0, -96.0, 0.0, 0.0);
    let (x, y) = lcc_forward(35.0_f64, -75.0, &c);
    assert_close(x, 1_894_410.9, 0.2);
    assert_close(y, 1_564_649.5, 0.2);
    let (lat, lon) = lcc_inverse(x, y, &c);
    assert_close(lat, 35.0, 1e-9);
    assert_close(lon, -75.0, 1e-9);
}

#[test]
fn albers_matches_snyder_and_round_trips() {
    // Snyder PP1395 p. 292: parallels 29.5°/45.5°, origin (23°, −96°),
    // point (35°, −75°) → (1 885 472.7, 1 535 925.0) m.
    let c = AlbersConstants::new(CLARKE_A, ecc(CLARKE_F), 29.5, 45.5, 23.0, -96.0, 0.0, 0.0);
    let (x, y) = albers_forward(35.0_f64, -75.0, &c);
    assert_close(x, 1_885_472.7, 0.2);
    assert_close(y, 1_535_925.0, 0.2);
    let (lat, lon) = albers_inverse(x, y, &c);
    assert_close(lat, 35.0, 1e-9);
    assert_close(lon, -75.0, 1e-9);
}

#[test]
fn polar_stereographic_round_trips() {
    let c = PolarStereoConstants::new(
        WGS84_A,
        ecc(WGS84_F),
        71.0,
        0.0,
        true,
        2_000_000.0,
        2_000_000.0,
    );
    for &(lat, lon) in &[(85.0_f64, 30.0), (72.0, -120.0), (89.0, 175.0)] {
        let (x, y) = polar_stereographic_forward(lat, lon, &c);
        let (blat, blon) = polar_stereographic_inverse(x, y, &c);
        assert_close(blat, lat, 1e-7);
        assert_close(blon, lon, 1e-6);
    }
}

#[test]
fn auxiliary_latitudes_round_trip() {
    let e = ecc(WGS84_F);
    for &lat in &[0.0_f64, 12.3, 45.0, -60.0, 80.0] {
        let phi = lat.to_radians();
        assert_close(
            lat_from_conformal_t(conformal_t(phi, e), e).to_degrees(),
            lat,
            1e-9,
        );
        assert_close(
            lat_from_authalic_q(authalic_q(phi, e), e).to_degrees(),
            lat,
            1e-9,
        );
    }
}

#[test]
fn families_run_for_f32() {
    // Genericity over `Real`: the family kernels run at `f32` (constants
    // stay `f64`, the coordinate is the lazy scalar) — so they lower to any
    // lazy backend, tensors included.
    let e = ecc(WGS84_F);
    let lcc = LccConstants::new(WGS84_A, e, 33.0, 45.0, 39.0, -96.0, 0.0, 0.0);
    let (x, y) = lcc_forward(40.0_f32, -100.0_f32, &lcc);
    assert!(x.is_finite() && y.is_finite());
    let alb = AlbersConstants::new(WGS84_A, e, 29.5, 45.5, 23.0, -96.0, 0.0, 0.0);
    let (ax, ay) = albers_forward(35.0_f32, -100.0_f32, &alb);
    assert!(ax.is_finite() && ay.is_finite());
    let ps = PolarStereoConstants::new(WGS84_A, e, 71.0, 0.0, true, 0.0, 0.0);
    let (px, py) = polar_stereographic_forward(85.0_f32, 30.0_f32, &ps);
    assert!(px.is_finite() && py.is_finite());
    let e32 = e as f32;
    assert!(m_factor(0.7_f32, e32).is_finite());
    assert!(conformal_t(0.7_f32, e32).is_finite());
    assert!(authalic_q(0.7_f32, e32).is_finite());
}

#[test]
fn typed_projection_constants_match_convenience_entries() {
    let web = WebMercatorConstants::from_radius(WGS84_A);
    let web_direct = web_mercator_forward(40.0_f64, -75.0, WGS84_A);
    let web_typed = web_mercator_forward_typed(40.0_f64, -75.0, &web);
    assert_close(web_typed.0, web_direct.0, 1e-12);
    assert_close(web_typed.1, web_direct.1, 1e-12);
    let web_p = (40.0_f64, -75.0_f64);
    let web_at = web_mercator_project_at_typed(&web_p, &web);
    assert_close(web_at.0, web_direct.0, 1e-12);
    assert_close(web_at.1, web_direct.1, 1e-12);
    let web_inv_direct = web_mercator_inverse(web_direct.0, web_direct.1, WGS84_A);
    let web_inv_typed = web_mercator_inverse_typed(web_typed.0, web_typed.1, &web);
    assert_close(web_inv_typed.0, web_inv_direct.0, 1e-12);
    assert_close(web_inv_typed.1, web_inv_direct.1, 1e-12);

    let (cm, ks) = utm(31);
    let tm = ks.lift::<f64>(cm, 0.0, 0.9996, 500_000.0, 0.0);
    let direct = transverse_mercator_forward(48.5_f64, 2.0, cm, 0.0, 0.9996, 500_000.0, 0.0, &ks);
    let typed = transverse_mercator_forward_typed(48.5_f64, 2.0, &tm);
    assert_close(typed.0, direct.0, 1e-12);
    assert_close(typed.1, direct.1, 1e-12);
    let p = (48.5_f64, 2.0_f64);
    let typed_at = transverse_mercator_project_at_typed(&p, &tm);
    assert_close(typed_at.0, direct.0, 1e-12);
    assert_close(typed_at.1, direct.1, 1e-12);
    let inv_direct =
        transverse_mercator_inverse(direct.0, direct.1, cm, 0.0, 0.9996, 500_000.0, 0.0, &ks);
    let inv_typed = transverse_mercator_inverse_typed(typed.0, typed.1, &tm);
    assert_close(inv_typed.0, inv_direct.0, 1e-12);
    assert_close(inv_typed.1, inv_direct.1, 1e-12);

    let e = ecc(WGS84_F);
    let lcc = LccConstants::new(WGS84_A, e, 33.0, 45.0, 39.0, -96.0, 0.0, 0.0);
    let lcc_t = lcc.lift::<f64>();
    let lcc_direct = lcc_forward(40.0_f64, -100.0, &lcc);
    let lcc_typed = lcc_forward_typed(40.0_f64, -100.0, &lcc_t);
    assert_close(lcc_typed.0, lcc_direct.0, 1e-12);
    assert_close(lcc_typed.1, lcc_direct.1, 1e-12);
    let lcc_inv_direct = lcc_inverse(lcc_direct.0, lcc_direct.1, &lcc);
    let lcc_inv_typed = lcc_inverse_typed(lcc_typed.0, lcc_typed.1, &lcc_t);
    assert_close(lcc_inv_typed.0, lcc_inv_direct.0, 1e-12);
    assert_close(lcc_inv_typed.1, lcc_inv_direct.1, 1e-12);

    let alb = AlbersConstants::new(WGS84_A, e, 29.5, 45.5, 23.0, -96.0, 0.0, 0.0);
    let alb_t = alb.lift::<f64>();
    let alb_direct = albers_forward(35.0_f64, -100.0, &alb);
    let alb_typed = albers_forward_typed(35.0_f64, -100.0, &alb_t);
    assert_close(alb_typed.0, alb_direct.0, 1e-12);
    assert_close(alb_typed.1, alb_direct.1, 1e-12);
    let alb_inv_direct = albers_inverse(alb_direct.0, alb_direct.1, &alb);
    let alb_inv_typed = albers_inverse_typed(alb_typed.0, alb_typed.1, &alb_t);
    assert_close(alb_inv_typed.0, alb_inv_direct.0, 1e-12);
    assert_close(alb_inv_typed.1, alb_inv_direct.1, 1e-12);

    let ps = PolarStereoConstants::new(WGS84_A, e, 71.0, 0.0, true, 0.0, 0.0);
    let ps_t = ps.lift::<f64>();
    let ps_direct = polar_stereographic_forward(85.0_f64, 30.0, &ps);
    let ps_typed = polar_stereographic_forward_typed(85.0_f64, 30.0, &ps_t);
    assert_close(ps_typed.0, ps_direct.0, 1e-12);
    assert_close(ps_typed.1, ps_direct.1, 1e-12);
    let ps_inv_direct = polar_stereographic_inverse(ps_direct.0, ps_direct.1, &ps);
    let ps_inv_typed = polar_stereographic_inverse_typed(ps_typed.0, ps_typed.1, &ps_t);
    assert_close(ps_inv_typed.0, ps_inv_direct.0, 1e-12);
    assert_close(ps_inv_typed.1, ps_inv_direct.1, 1e-12);
}
