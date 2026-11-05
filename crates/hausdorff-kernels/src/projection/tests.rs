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

// ---- sphere degeneracy in the auxiliary latitudes ----
//
// A sphere (`e = 0`) is a legitimate ellipsoid, and `Ellipsoid::SPHERE` is a public
// constant, so every auxiliary latitude has to have a defined value there. The
// authalic pair did not: Snyder 3-12 and 3-16 carry a literal `1/(2e)`, which is a
// pole at `e = 0`, and the result was `inf · ln(1)` — NaN for every equal-area
// quantity on a sphere, including the whole Albers projection. These pin the limits.

#[test]
fn auxiliary_latitudes_are_finite_on_a_sphere() {
    // Every auxiliary latitude must be finite at e = 0, not just the ones whose
    // formula happens to avoid a division.
    let mut phi = -1.5_f64;
    while phi <= 1.5 {
        for f in [
            authalic_q(phi, 0.0),
            lat_from_authalic_q(phi.sin() * 2.0, 0.0),
            conformal_t(phi, 0.0),
            lat_from_conformal_t(phi.tan(), 0.0),
            m_factor(phi, 0.0),
            taupf(phi.tan(), 0.0),
            tauf(phi.tan(), 0.0),
        ] {
            assert!(
                f.is_finite(),
                "non-finite auxiliary latitude at phi = {phi}: {f}"
            );
        }
        phi += 0.1;
    }
}

#[test]
fn authalic_q_on_a_sphere_is_twice_the_sine() {
    // The e -> 0 limit of Snyder 3-12. Expanding atanh(e·s)/e as a power series
    // gives s + e²s³/3 + ..., and s/(1 - e²s²) gives s + e²s³ + ..., so q -> 2s.
    // With q(pi/2) = 2 this makes sin(xi) = sin(phi): the authalic latitude
    // collapses to the geodetic one, which is what a sphere means.
    assert_close(authalic_q(std::f64::consts::FRAC_PI_2, 0.0), 2.0, 1e-15);
    let mut phi = -1.5_f64;
    while phi <= 1.5 {
        assert_close(authalic_q(phi, 0.0), 2.0 * phi.sin(), 1e-15);
        // And the inverse recovers phi exactly.
        assert_close(lat_from_authalic_q(2.0 * phi.sin(), 0.0), phi, 1e-12);
        phi += 0.1;
    }
}

#[test]
fn authalic_q_is_continuous_as_eccentricity_approaches_zero() {
    // The old form did not merely fail *at* e = 0 — it lost about eps/e relative
    // precision *near* it, because it took the log of a ratio close to 1 and then
    // divided by a small number. Approaching the limit from above must converge
    // smoothly to 2·sin(phi) rather than degrading.
    let phi = 0.7_f64;
    let limit = 2.0 * phi.sin();
    // The leading correction is second order, not first: expanding both terms gives
    // q = 2s + e²(4s³/3 - 2s) + O(e⁴). So the bound is O(e²) with that exact
    // coefficient, which also means this fails if the convergence *rate* degrades
    // rather than only if the limit is wrong.
    let s = phi.sin();
    let coefficient = (4.0 * s * s * s / 3.0 - 2.0 * s).abs();
    for e in [1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-14, 0.0] {
        let q = authalic_q(phi, e);
        let bound = 1.01 * coefficient * e * e + 1e-15;
        assert!(
            (q - limit).abs() <= bound,
            "authalic_q({phi}, {e}) = {q} deviates from {limit} by {}, above the \
             second-order bound {bound}",
            (q - limit).abs()
        );
    }
}

#[test]
fn authalic_q_still_matches_snyder_on_a_real_ellipsoid() {
    // The pole removal rewrote the formula, so the non-degenerate case needs an
    // anchor: the authalic latitude of 40 deg on WGS84, and the monotonicity and
    // antisymmetry that hold for any eccentricity.
    let e = (WGS84_F * (2.0 - WGS84_F)).sqrt();
    // At the pole s = 1, so Snyder 3-12 collapses to the closed form
    // q(pi/2) = 1 + (1 - e²)·atanh(e)/e. Written out independently here — a
    // different expression from the one the kernel evaluates — and cross-checked
    // against its decimal value, so a rewrite of the kernel cannot drag both sides
    // along with it.
    let q_pole = authalic_q(std::f64::consts::FRAC_PI_2, e);
    let closed_form = 1.0 + (1.0 - e * e) * e.atanh() / e;
    assert_close(q_pole, closed_form, 1e-15);
    assert_close(q_pole, 1.995_531_087_502_837, 1e-14);
    // Antisymmetry: q(-phi) = -q(phi), exactly.
    let mut phi = 0.05_f64;
    while phi < 1.5 {
        assert_close(authalic_q(-phi, e), -authalic_q(phi, e), 1e-15);
        phi += 0.1;
    }
    // Strict monotonicity, and the round trip.
    let mut prev = f64::NEG_INFINITY;
    let mut i = -89;
    while i <= 89 {
        let p = f64::from(i).to_radians();
        let q = authalic_q(p, e);
        assert!(q > prev, "q must increase with latitude, broke at {i} deg");
        prev = q;
        assert_close(lat_from_authalic_q(q, e), p, 1e-11);
        i += 1;
    }
}

// ---- conic injectivity across the antimeridian ----
//
// A conic scales the longitude difference by a cone constant and feeds it to
// sin/cos, so the cone angle is `n·dlon`. With `dlon` unwrapped it reaches +-360,
// and for a high-latitude conic (n near 1) the cone angle then passes +-180: the
// forward wraps and the inverse's atan2 recovers a *different meridian*. EPSG:3978
// (Canada Atlas Lambert, lon0 = -95) mapped longitude 180 and recovered 140.35 --
// finite, plausible, and 40 degrees wrong. Both halves are now normalized: the
// difference on the way in, the sum on the way out.

#[test]
fn central_lon_delta_is_bit_exact_inside_half_a_turn() {
    // The reason the wrap is selected rather than applied unconditionally. The
    // algebraic identity `(x + 180) % 360 - 180` is not bit-exact for small x --
    // at 1e-10 it returns 1.0000000827e-10 -- so applying it to every coordinate
    // would shift results this crate pins against PROJ.
    for raw in [
        0.0_f64,
        1e-12,
        1e-10,
        1e-7,
        0.5,
        1.0,
        45.0,
        90.0,
        179.0,
        179.999_999_999,
        180.0,
    ] {
        for sign in [1.0_f64, -1.0] {
            let d = raw * sign;
            assert_eq!(
                central_lon_delta(d, 0.0).to_bits(),
                d.to_bits(),
                "a difference of {d} is already in range and must pass through untouched"
            );
        }
    }
}

#[test]
fn central_lon_delta_wraps_only_what_is_out_of_range() {
    for (lon, lon0, expected) in [
        (180.0, -95.0, -85.0), // the EPSG:3978 case: 275 -> -85
        (179.0, -95.0, -86.0),
        (-180.0, 95.0, 85.0),   // -275 -> 85
        (170.0, -170.0, -20.0), // 340 -> -20
        (-170.0, 170.0, 20.0),
        (90.0, -179.0, -91.0), // 269 -> -91
    ] {
        let got = central_lon_delta(lon, lon0);
        assert_close(got, expected, 1e-9);
        assert!(
            got.abs() <= 180.0,
            "a wrapped difference must land in [-180, 180], got {got}"
        );
        // And it is the same meridian: the two differ by a whole turn.
        let raw = lon - lon0;
        assert_close((raw - got).abs(), 360.0, 1e-9);
    }
}

#[test]
fn conics_round_trip_from_any_longitude() {
    // The property the normalization buys: injectivity over the whole sphere
    // rather than only within 180/n of the central meridian. Swept over central
    // meridians and cone constants, since the defect only bites when n is large
    // and lon0 is far from the point.
    let e = (WGS84_F * (2.0 - WGS84_F)).sqrt();
    for lon0 in [-179.0, -95.0, 0.0, 95.0, 179.0] {
        // Standard parallels near the pole give n close to 1, the worst case.
        let lcc = LccConstants::new(WGS84_A, e, 77.0, 83.0, 80.0, lon0, 0.0, 0.0);
        let albers = AlbersConstants::new(WGS84_A, e, 77.0, 83.0, 80.0, lon0, 0.0, 0.0);
        let mut lon = -180.0_f64;
        while lon <= 180.0 {
            for lat in [60.0_f64, 75.0, 80.0, 85.0] {
                let (x, y) = lcc_forward(lat, lon, &lcc);
                assert!(
                    x.is_finite() && y.is_finite(),
                    "LCC lon0={lon0} at ({lat}, {lon}) produced ({x}, {y})"
                );
                let (rlat, rlon) = lcc_inverse(x, y, &lcc);
                assert_close(rlat, lat, 1e-7);
                // Compare as meridians: +-180 are the same place.
                let d = (rlon - lon).rem_euclid(360.0);
                let d = if d > 180.0 { 360.0 - d } else { d };
                assert!(
                    d < 1e-7,
                    "LCC lon0={lon0}: ({lat}, {lon}) recovered longitude {rlon}, \
                     off by {d} degrees"
                );
                assert!(
                    (-180.0..=180.0).contains(&rlon),
                    "LCC lon0={lon0}: recovered longitude {rlon} is outside [-180, 180]"
                );

                let (ax, ay) = albers_forward(lat, lon, &albers);
                let (alat, alon) = albers_inverse(ax, ay, &albers);
                assert_close(alat, lat, 1e-7);
                let da = (alon - lon).rem_euclid(360.0);
                let da = if da > 180.0 { 360.0 - da } else { da };
                assert!(
                    da < 1e-7,
                    "Albers lon0={lon0}: ({lat}, {lon}) recovered longitude {alon}, \
                     off by {da} degrees"
                );
                assert!(
                    (-180.0..=180.0).contains(&alon),
                    "Albers lon0={lon0}: recovered longitude {alon} is outside [-180, 180]"
                );
            }
            lon += 7.0;
        }
    }
}

#[test]
fn conic_inverses_return_canonical_longitudes() {
    // The inverse recovers a *difference* and adds lon0 back, which can leave the
    // sum outside (-180, 180] even when the meridian is right. A caller comparing
    // against its input would reject -181 for 179.
    let e = (WGS84_F * (2.0 - WGS84_F)).sqrt();
    let lcc = LccConstants::new(WGS84_A, e, 29.5, 45.5, 23.0, -96.0, 0.0, 0.0);
    for lon in [179.0_f64, 180.0, -179.0, -180.0, 150.0, -150.0] {
        let (x, y) = lcc_forward(40.0, lon, &lcc);
        let (_, rlon) = lcc_inverse(x, y, &lcc);
        assert!(
            (-180.0..=180.0).contains(&rlon),
            "longitude {lon} recovered as {rlon}, outside the canonical range"
        );
    }
}
