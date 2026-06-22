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

/// The authalic-area function `q` (Snyder 3-12) — the equal-area analogue of the
/// conformal `t`, with `sin ξ = q / q(π/2)`.
#[must_use]
pub fn authalic_q<T: Real>(phi: T, e: T) -> T {
    let one = T::one();
    let two = T::two();
    let e2 = e.clone() * e.clone();
    let one_minus_e2 = one.clone() - e2;
    let inv_2e = one.clone() / (two * e.clone());
    let es = phi.clone().sin() * e;
    let one_m_es2 = one.clone() - es.clone() * es.clone();
    one_minus_e2 * (phi.sin() / one_m_es2 - inv_2e * ((one.clone() - es.clone()) / (one + es)).ln())
}

/// Invert the authalic `q` back to geodetic latitude `φ` (radians) by Snyder's
/// iteration (3-16). The pole clamp `|q| ≥ q(π/2)` is left to the `f64` wrapper
/// (a scalar branch); this kernel runs the fixed-count fixed point.
#[must_use]
pub fn lat_from_authalic_q<T: Real>(q: T, e: T) -> T {
    let one = T::one();
    let two = T::two();
    let half = T::half();
    let e2 = e.clone() * e.clone();
    let one_minus_e2 = one.clone() - e2;
    let inv_2e = one.clone() / (two.clone() * e.clone());
    let mut phi = (q.clone() * half).asin();
    for _ in 0..15 {
        let es = phi.clone().sin() * e.clone();
        let one_m_es2 = one.clone() - es.clone() * es.clone();
        let delta = one_m_es2.clone() * one_m_es2.clone() / (two.clone() * phi.clone().cos())
            * (q.clone() / one_minus_e2.clone() - phi.clone().sin() / one_m_es2
                + inv_2e.clone() * ((one.clone() - es.clone()) / (one.clone() + es)).ln());
        phi = phi + delta;
    }
    phi
}
