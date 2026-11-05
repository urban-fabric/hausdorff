use super::*;

// ---- Shared auxiliary-latitude helpers (Karney 2011) ----

/// `tan χ` from `tan φ`: the tangent of the **conformal latitude** χ given the
/// tangent of the geodetic latitude φ, on an ellipsoid of eccentricity `e`
/// (Karney 2011 `taupf`). The exact closed form (no series).
#[must_use]
pub fn taupf<T: Real>(tau: T, e: T) -> T {
    let one = T::one();
    let tau1 = tau.clone().hypot(one.clone());
    let sig = (e.clone() * (e * tau.clone() / tau1.clone()).atanh()).sinh();
    sig.clone().hypot(one) * tau - sig * tau1
}

/// The inverse of [`taupf`]: `tan φ` from `tan χ`, by Newton's method (Karney
/// 2011 `tauf`).
///
/// Runs a **fixed** five iterations — enough for full `f64` precision — rather
/// than breaking on a scalar tolerance, so the kernel lowers to a tensor
/// backend (which has no scalar `bool` to branch on) unchanged.
#[must_use]
pub fn tauf<T: Real>(taup: T, e: T) -> T {
    let one = T::one();
    let e2m = one.clone() - e.clone() * e.clone(); // 1 − e²
    let mut tau = taup.clone() / e2m.clone(); // initial guess
    for _ in 0..5 {
        let taupa = taupf(tau.clone(), e.clone());
        let dtau = (taup.clone() - taupa.clone())
            * (one.clone() + e2m.clone() * tau.clone() * tau.clone())
            / (e2m.clone() * tau.clone().hypot(one.clone()) * taupa.hypot(one.clone()));
        tau = tau + dtau;
    }
    tau
}

// ---- Auxiliary-latitude kernels (Snyder PP1395) ----
//
// The shared latitude transforms the conic / azimuthal / equal-area projection
// families are written in. Radians in/out; eccentricity is carried in the same
// lazy scalar as the coordinate. The inverse iterations run a **fixed** count
// rather than breaking on a scalar tolerance, so they lower to a tensor backend
// unchanged (they already converge to full `f64` precision well inside the bound).

/// `m = cosφ / √(1 − e² sin²φ)` (Snyder 14-15) — the radius of the parallel
/// scaled to the unit sphere; the prime-vertical-radius factor.
#[must_use]
pub fn m_factor<T: Real>(phi: T, e: T) -> T {
    let one = T::one();
    let e2 = e.clone() * e;
    let s = phi.clone().sin();
    phi.cos() / (one - e2 * s.clone() * s).sqrt()
}

/// The conformal isometric building block `t = tan(π/4 − φ/2) / [(1 − e sinφ)/(1 +
/// e sinφ)]^{e/2}` (Snyder 15-9) — `t = tan(π/4 − χ/2)`, so `χ = π/2 − 2·arctan t`.
#[must_use]
pub fn conformal_t<T: Real>(phi: T, e: T) -> T {
    let one = T::one();
    let half = T::half();
    let frac_pi_4 = T::frac_pi_4();
    let es = phi.clone().sin() * e.clone();
    let exponent = e * half.clone();
    (frac_pi_4 - phi * half).tan() / ((one.clone() - es.clone()) / (one + es)).powf(exponent)
}

/// Invert the conformal `t` back to geodetic latitude `φ` (radians) by Snyder's
/// fixed-point iteration (7-9). Its `e → 0` limit is the exact closed form
/// `φ = π/2 − 2·arctan t`.
#[must_use]
pub fn lat_from_conformal_t<T: Real>(t: T, e: T) -> T {
    let one = T::one();
    let two = T::two();
    let frac_pi_2 = T::frac_pi_2();
    let mut phi = frac_pi_2.clone() - two.clone() * t.clone().atan();
    for _ in 0..15 {
        let es = phi.clone().sin() * e.clone();
        let exponent = e.clone() * T::half();
        phi = frac_pi_2.clone()
            - two.clone()
                * (t.clone() * ((one.clone() - es.clone()) / (one.clone() + es)).powf(exponent))
                    .atan();
    }
    phi
}

/// `atanh(e·s)/e`, with the removable singularity at `e = 0` filled by its limit, `s`.
///
/// This is Snyder's `−(1/2e)·ln[(1 − e·s)/(1 + e·s)]` rewritten through
/// `ln[(1 + x)/(1 − x)] = 2·atanh(x)`. The rewrite is not cosmetic. The `1/e` factor
/// is a genuine pole at `e = 0`, and a sphere is a legitimate ellipsoid here
/// ([`Ellipsoid::SPHERE`](../../../hausdorff_base/earth/struct.Ellipsoid.html)), so the
/// original form returned `NaN` for every equal-area quantity on one: `inf · ln(1) =
/// inf · 0`. It is also more accurate at *small* non-zero `e`, where the original
/// takes the log of a ratio near 1 and then divides by a small number — losing about
/// `ε/e` relative precision (≈1e-8 at `e = 1e-8`) that `atanh`, built for exactly this
/// argument range, keeps.
fn atanh_over_e<T: Real>(s: T, e: T) -> T {
    let is_sphere = e.eq(&T::zero());
    // Sanitize the divisor *before* dividing rather than repairing the result after —
    // euclid's "double where". An eager elementwise carrier then never forms the `0/0`
    // in a sphere lane at all, instead of forming a NaN and discarding it. `half` is
    // an arbitrary safe stand-in that also keeps `atanh`'s argument inside (−1, 1) for
    // any |s| ≤ 1, so the dead lane cannot overflow either.
    let safe = T::select(&is_sphere, T::half(), e);
    T::select(&is_sphere, s.clone(), (safe.clone() * s).atanh() / safe)
}

/// The authalic-area function `q` (Snyder 3-12) — the equal-area analogue of the
/// conformal `t`, with `sin ξ = q / q(π/2)`. On a sphere this is `2·sinφ`, so
/// `sin ξ = sinφ` and the authalic latitude collapses to the geodetic one.
#[must_use]
pub fn authalic_q<T: Real>(phi: T, e: T) -> T {
    let one = T::one();
    let one_minus_e2 = one.clone() - e.clone() * e.clone();
    let s = phi.sin();
    let es = s.clone() * e.clone();
    let one_m_es2 = one - es.clone() * es;
    one_minus_e2 * (s.clone() / one_m_es2 + atanh_over_e(s, e))
}

/// Invert the authalic `q` back to geodetic latitude `φ` (radians) by Snyder's
/// iteration (3-16). The pole clamp `|q| ≥ q(π/2)` is left to the `f64` wrapper
/// (a scalar branch); this kernel runs the fixed-count fixed point.
#[must_use]
pub fn lat_from_authalic_q<T: Real>(q: T, e: T) -> T {
    let one = T::one();
    let two = T::two();
    let half = T::half();
    let one_minus_e2 = one.clone() - e.clone() * e.clone();
    let mut phi = (q.clone() * half).asin();
    for _ in 0..15 {
        let s = phi.clone().sin();
        let es = s.clone() * e.clone();
        let one_m_es2 = one.clone() - es.clone() * es;
        // `+ (1/2e)·ln[(1 − e·s)/(1 + e·s)]` is `− atanh(e·s)/e`; see `atanh_over_e`
        // for why the pole is not written literally.
        let delta = one_m_es2.clone() * one_m_es2.clone() / (two.clone() * phi.clone().cos())
            * (q.clone() / one_minus_e2.clone()
                - s.clone() / one_m_es2
                - atanh_over_e(s, e.clone()));
        phi = phi + delta;
    }
    phi
}

/// The longitude difference `lon − lon0`, wrapped into `[−180, 180]`.
///
/// Conic and azimuthal projections scale this difference by a cone constant `n`
/// and feed the result to `sin`/`cos`, so the *cone angle* is `n·Δλ`. If `Δλ`
/// arrives unwrapped it can reach ±360, and for a high-latitude conic (where `n`
/// approaches 1) the cone angle then passes ±180 — at which point the forward
/// projection wraps and the inverse's `atan2` recovers a **different meridian**.
/// The result is a finite, plausible, wrong coordinate: EPSG:3978 (Canada Atlas
/// Lambert, `lon0 = −95°`) mapped longitude 180° and recovered 140.35°.
///
/// Wrapping first makes `|Δλ| ≤ 180`, so `|n·Δλ| < 180` for any `n < 1` and the
/// projection is injective over the whole sphere rather than only within `180/n`
/// of its central meridian.
///
/// The wrap is **selected rather than computed unconditionally**, which is the
/// detail that matters. `(x + 180) mod 360 − 180` is algebraically the identity
/// for `|x| ≤ 180` but not *bit-exact*: at `x = 1e-10` the intermediate `180.0000000001`
/// has already lost the low bits, and the round trip returns `1.0000000827e-10`.
/// Applying it unconditionally would perturb every in-domain coordinate — and
/// shift results this crate pins against PROJ reference values — to repair a case
/// that only arises outside the domain. So the common path returns `raw`
/// untouched.
#[must_use]
pub fn central_lon_delta<T: Real>(lon_deg: T, lon0_deg: T) -> T {
    let raw = lon_deg - lon0_deg;
    let half = T::half_turn_degrees();
    let needs_wrap = raw.clone().abs().gt(&half);
    let wrapped = (raw.clone() + half.clone()).rem_euclid(T::full_turn_degrees()) - half;
    T::select(&needs_wrap, wrapped, raw)
}
