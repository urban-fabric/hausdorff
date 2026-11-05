//! The crate's central claim, tested: these formulas are written once and run on
//! any `euclid-numeric` carrier.
//!
//! Every other test in this crate instantiates the kernels at `f64`, which means
//! the genericity was *claimed and not backed*. A function can be generic in its
//! signature and still be wrong or uncompilable at another carrier — a stray
//! `as f64`, an `f64::consts` reference, a `0.5` literal where `T::half()` was
//! needed. None of that shows up until something instantiates it.
//!
//! So this suite exists to make the claim falsifiable. It instantiates the public
//! surface at **`f32`** and checks each result against the `f64` answer. `f32` is
//! the right second carrier for this: it is a genuinely different type with
//! genuinely different precision, so it catches the mistakes above, while still
//! being a host scalar whose answers can be compared numerically. A tensor or
//! expression carrier would test compilation but have no number to compare.
//!
//! **Tolerances are quoted in `f32` ulp**, not as bare numbers. `f32` has about
//! 7 decimal digits, so a relative agreement of ~1e-7 is *exact* at this
//! precision and a bare `1e-7` would read as suspiciously loose. Quoting ulp makes
//! the question legible: how many trailing bits does single precision lose here?
//! A composition of several transcendentals loses a few, which is why the bounds
//! below are tens of ulp rather than one.
//!
//! What this suite deliberately does **not** claim is that the whole workspace is
//! carrier-generic. It is not, and cannot be. The planar, spherical, and index
//! layers branch on comparisons, sort, and use ordered containers, none of which a
//! carrier without a host `bool` can do. That boundary is the reason this crate
//! exists as a separate layer, and it is documented in the crate header.

use crate::earth::{haversine_m, initial_bearing_degrees};
use crate::planar::{distance_squared2d, orient2d, point_segment_distance_squared2d};
use crate::projection::{
    KrugerSeries, authalic_q, conformal_t, lat_from_authalic_q, m_factor, tauf, taupf,
    transverse_mercator_forward, transverse_mercator_inverse, web_mercator_forward,
    web_mercator_inverse,
};
use crate::spherical1::normalize_radians;
use euclid_numeric::interval;

/// One `f32` ulp, relative.
const ULP: f64 = 1.192_092_895_507_812_5e-7;

/// Assert an `f32` result agrees with its `f64` peer to within `ulps` of single
/// precision, under the mixed absolute/relative criterion that degrades
/// gracefully through zero.
#[track_caller]
fn agrees(got: f32, want: f64, ulps: f64, what: &str) {
    let err = (f64::from(got) - want).abs();
    let bound = ulps * ULP * want.abs().max(1.0);
    assert!(
        err <= bound,
        "{what}: f32 {got:e} vs f64 {want:e} (error {err:.3e} > {bound:.3e}, {ulps} ulp)"
    );
}

/// Places spread across the hemispheres and the seams, as (lat, lon) degrees.
const PLACES: &[(f32, f32)] = &[
    (51.5074, -0.1278),
    (48.8566, 2.3522),
    (-33.8688, 151.2093),
    (35.6762, 139.6503),
    (0.0, 0.0),
    (10.0, 179.5),
    (-10.0, -179.5),
    (60.0, -45.0),
    (-60.0, 45.0),
];

#[test]
fn haversine_and_bearing_agree_at_f32() {
    for &(lat1, lon1) in PLACES {
        for &(lat2, lon2) in PLACES {
            let d32 = haversine_m(lat1, lon1, lat2, lon2);
            let d64 = haversine_m(
                f64::from(lat1),
                f64::from(lon1),
                f64::from(lat2),
                f64::from(lon2),
            );
            // A haversine composes sin, cos, sqrt, and atan2, so a handful of
            // trailing bits go. Tens of ulp on a 10 000 km distance is a few
            // metres, which is what single precision buys on a planet.
            agrees(d32, d64, 64.0, "haversine");

            if (lat1, lon1) != (lat2, lon2) {
                let b32 = initial_bearing_degrees(lat1, lon1, lat2, lon2);
                let b64 = initial_bearing_degrees(
                    f64::from(lat1),
                    f64::from(lon1),
                    f64::from(lat2),
                    f64::from(lon2),
                );
                // A bearing is an atan2 of two differences, so near-antipodal or
                // near-coincident pairs are genuinely ill-conditioned; the bound
                // is on the angle, where 1e-3 degrees is 4 metres at 1 km.
                let delta = (f64::from(b32) - b64)
                    .abs()
                    .min(360.0 - (f64::from(b32) - b64).abs());
                assert!(
                    delta < 1e-2,
                    "bearing: f32 {b32} vs f64 {b64} differ by {delta} degrees"
                );
            }
        }
    }
}

#[test]
fn web_mercator_agrees_and_round_trips_at_f32() {
    let radius32 = 6_378_137.0_f32;
    let radius64 = 6_378_137.0_f64;
    for &(lat, lon) in PLACES {
        // Web Mercator is undefined at the poles and our sample stays well inside.
        let (x32, y32) = web_mercator_forward(lat, lon, radius32);
        let (x64, y64) = web_mercator_forward(f64::from(lat), f64::from(lon), radius64);
        agrees(x32, x64, 32.0, "web mercator x");
        agrees(y32, y64, 32.0, "web mercator y");

        // And the inverse must undo the forward *within the same carrier*, which
        // is a different claim from matching f64: it says the pair is consistent
        // at this precision rather than that it tracks a wider one.
        let (rlat, rlon) = web_mercator_inverse(x32, y32, radius32);
        assert!(
            (rlat - lat).abs() < 1e-3 && (rlon - lon).abs() < 1e-3,
            "web mercator round trip at f32: ({lat}, {lon}) -> ({rlat}, {rlon})"
        );
    }
}

#[test]
fn transverse_mercator_agrees_and_round_trips_at_f32() {
    // The Krüger series constants are computed once in f64 by design — they are
    // ellipsoid parameters, not coordinates — and lifted into the carrier. This is
    // the case most likely to have an `as f64` hiding in it.
    let series = KrugerSeries::for_ellipsoid(6_378_137.0, 1.0 / 298.257_223_563);
    for &(lat, lon) in &PLACES[..4] {
        let (x32, y32) = transverse_mercator_forward(
            lat, lon, -3.0, 0.0, 0.9996, 400_000.0, -100_000.0, &series,
        );
        let (x64, y64) = transverse_mercator_forward(
            f64::from(lat),
            f64::from(lon),
            -3.0,
            0.0,
            0.9996,
            400_000.0,
            -100_000.0,
            &series,
        );
        agrees(x32, x64, 64.0, "transverse mercator easting");
        agrees(y32, y64, 64.0, "transverse mercator northing");

        let (rlat, rlon) = transverse_mercator_inverse(
            x32, y32, -3.0, 0.0, 0.9996, 400_000.0, -100_000.0, &series,
        );
        assert!(
            (rlat - lat).abs() < 1e-2 && (rlon - lon).abs() < 1e-2,
            "transverse mercator round trip at f32: ({lat}, {lon}) -> ({rlat}, {rlon})"
        );
    }
}

#[test]
fn auxiliary_latitudes_agree_at_f32() {
    // The fixed-count iterations are the thing to check here: `tauf` runs five
    // Newton steps and `lat_from_authalic_q` fifteen fixed-point steps, counts
    // chosen for full f64 precision. They must still converge at f32 rather than
    // oscillating or stalling.
    let e32 = 0.081_819_19_f32;
    let e64 = 0.081_819_190_842_621_f64;
    let mut i = -85;
    while i <= 85 {
        let phi32 = (i as f32).to_radians();
        let phi64 = f64::from(i).to_radians();

        agrees(
            authalic_q(phi32, e32),
            authalic_q(phi64, e64),
            32.0,
            "authalic q",
        );
        agrees(
            conformal_t(phi32, e32),
            conformal_t(phi64, e64),
            64.0,
            "conformal t",
        );
        agrees(m_factor(phi32, e32), m_factor(phi64, e64), 16.0, "m factor");
        agrees(
            taupf(phi32.tan(), e32),
            taupf(phi64.tan(), e64),
            64.0,
            "taupf",
        );

        // The inverses, as round trips within the carrier.
        let q32 = authalic_q(phi32, e32);
        let back = lat_from_authalic_q(q32, e32);
        assert!(
            (back - phi32).abs() < 1e-4,
            "authalic round trip at f32: {phi32} -> {back}"
        );
        let t32 = taupf(phi32.tan(), e32);
        let tan_back = tauf(t32, e32);
        assert!(
            (tan_back - phi32.tan()).abs() <= 1e-3 * phi32.tan().abs().max(1.0),
            "tauf round trip at f32: tan {} -> {tan_back}",
            phi32.tan()
        );
        i += 5;
    }

    // The sphere degeneracy has to hold at f32 too: e = 0 must not be a pole.
    assert!(
        authalic_q(0.5_f32, 0.0).is_finite(),
        "authalic_q must be finite on a sphere at f32 as well"
    );
    agrees(
        authalic_q(0.5_f32, 0.0),
        2.0 * 0.5_f64.sin(),
        16.0,
        "sphere authalic q",
    );
}

#[test]
fn planar_and_interval_kernels_agree_at_f32() {
    let pts: &[(f32, f32)] = &[
        (0.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (-3.5, 2.25),
        (1e4, -1e4),
    ];
    for &(ax, ay) in pts {
        for &(bx, by) in pts {
            agrees(
                distance_squared2d(ax, ay, bx, by),
                distance_squared2d(f64::from(ax), f64::from(ay), f64::from(bx), f64::from(by)),
                16.0,
                "distance squared",
            );
            for &(cx, cy) in pts {
                let s32 = orient2d(ax, ay, bx, by, cx, cy);
                let s64 = orient2d(
                    f64::from(ax),
                    f64::from(ay),
                    f64::from(bx),
                    f64::from(by),
                    f64::from(cx),
                    f64::from(cy),
                );
                // Signs agree only where the determinant is *well-conditioned*,
                // and the condition is not about how far apart the points are. It
                // is about how small the determinant is relative to the products
                // that formed it: `orient2d` is a difference of two products, so
                // the f32 error floor is roughly one f32 ulp of the larger
                // product, regardless of the spacing of the inputs.
                let magnitude = f64::from(ax).abs().max(f64::from(bx).abs())
                    * f64::from(ay).abs().max(f64::from(by).abs());
                let floor = 8.0 * ULP * magnitude.max(1.0);
                if s64.abs() > floor {
                    assert_eq!(
                        s32.partial_cmp(&0.0),
                        s64.partial_cmp(&0.0_f64),
                        "orient2d sign must agree where the determinant is above the \
                         f32 error floor {floor:e}: ({ax},{ay}) ({bx},{by}) ({cx},{cy}) \
                         gave {s32} vs {s64}"
                    );
                }
            }
            // Bounded against the *input* magnitude, not the result. Projecting a
            // point onto a segment forms products of the coordinates and then
            // subtracts, so the error floor is set by the largest intermediate —
            // and for a segment anchored at 1e4 that is 1e8, whatever the answer
            // happens to be. Scaling the bound by the result instead would demand
            // that a small distance between distant points be computed to a
            // precision the arithmetic never had.
            let got = point_segment_distance_squared2d(0.5, 0.5, ax, ay, bx, by);
            let want = point_segment_distance_squared2d(
                0.5,
                0.5,
                f64::from(ax),
                f64::from(ay),
                f64::from(bx),
                f64::from(by),
            );
            let span = f64::from(ax)
                .abs()
                .max(f64::from(ay).abs())
                .max(f64::from(bx).abs())
                .max(f64::from(by).abs())
                .max(1.0);
            let bound = 32.0 * ULP * (span * span).max(want.abs());
            assert!(
                (f64::from(got) - want).abs() <= bound,
                "point-segment distance squared: f32 {got:e} vs f64 {want:e} \
                 (error {:.3e} > {bound:.3e} for a span of {span:e})",
                (f64::from(got) - want).abs()
            );
        }
    }

    // Masks: for a host scalar `Compare::Mask` is `bool`, which is what lets these
    // kernels be used in ordinary control flow at all. Pinned because it is the
    // property the whole two-layer split rests on.
    assert!(interval::contains(0.0_f32, 1.0, 0.5));
    assert!(!interval::contains(0.0_f32, 1.0, 1.5));
    let (lo, hi) = interval::intersection(0.0_f32, 2.0, 1.0, 3.0);
    agrees(lo, 1.0, 1.0, "interval intersection low");
    agrees(hi, 2.0, 1.0, "interval intersection high");
}

#[test]
fn angular_normalization_agrees_at_f32() {
    let mut i = -20;
    while i <= 20 {
        let x32 = i as f32 * 0.7;
        let n32 = normalize_radians(x32);
        let n64 = normalize_radians(f64::from(x32));
        agrees(n32, n64, 32.0, "normalize radians");
        // The documented range is `(-pi, pi]`, not `[0, tau)` — the kernel does a
        // `rem_euclid` and then folds the upper half down. Worth pinning because
        // the two conventions are equally plausible and a caller that assumed the
        // other would be wrong by a full turn on half its inputs.
        assert!(
            n32 > -std::f32::consts::PI - 1e-6 && n32 <= std::f32::consts::PI + 1e-6,
            "normalized angle {n32} escaped (-pi, pi]"
        );
        i += 1;
    }
}

#[test]
fn f32_loses_an_orientation_sign_that_exact_arithmetic_keeps() {
    // The boundary of what a carrier-generic predicate can promise, as a worked
    // case rather than a caveat.
    //
    // For a = (10000, -10000), b = (1, 0), c = (0, 1) the determinant is
    //
    //     (1 - 10000)(1 + 10000) - (0 + 10000)(0 - 10000)
    //       = -99 999 999 + 100 000 000
    //       = 1
    //
    // *exactly* — but it is the difference of two products near 1e8. f32 carries
    // about seven digits, so everything below ~10 is already gone before the
    // subtraction happens: f32 returns 0 and reports three points as collinear
    // that are not.
    //
    // Note what does *not* rescue this. The points are far apart, the inputs are
    // exactly representable, and nothing is near-degenerate geometrically. The
    // conditioning of a determinant is not bounded by the separation of its
    // points, which is why "use well-separated input" is not a usable mitigation
    // and why `hausdorff-planar::space::r2::exact` computes signs in i128 instead.
    //
    // f64 survives this particular case — it has the digits to spare — but the
    // same construction scaled up defeats it too. The lesson is about the
    // technique, not about f32.
    let (ax, ay) = (10_000.0_f32, -10_000.0_f32);
    let (bx, by) = (1.0_f32, 0.0_f32);
    let (cx, cy) = (0.0_f32, 1.0_f32);

    let single = orient2d(ax, ay, bx, by, cx, cy);
    let double = orient2d(
        f64::from(ax),
        f64::from(ay),
        f64::from(bx),
        f64::from(by),
        f64::from(cx),
        f64::from(cy),
    );
    assert_eq!(single, 0.0, "f32 cancels the determinant to nothing");
    assert_eq!(double, 1.0, "the true determinant is exactly 1");

    // The exact integer predicate in `hausdorff-planar` is the answer to this, and
    // it gets the sign right at any scale. Reproduced here as the arithmetic
    // rather than the dependency, since this crate deliberately has no geometry
    // types: i128 products of the same differences cannot cancel.
    let (iax, iay) = (10_000i128, -10_000i128);
    let (ibx, iby) = (1i128, 0i128);
    let (icx, icy) = (0i128, 1i128);
    let exact = (ibx - iax) * (icy - iay) - (iby - iay) * (icx - iax);
    assert_eq!(exact, 1, "integer arithmetic keeps the sign f32 threw away");
}
