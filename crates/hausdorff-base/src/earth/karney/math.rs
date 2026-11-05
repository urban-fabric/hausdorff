//! Scalar and angular primitives shared by the geodesic solver: the GeographicLib
//! tolerances, error-free transforms, polynomial/Clenshaw evaluation, degree-mode
//! trig, and the astroid root used to seed the nearly-antipodal case.
//!
//! The two pieces of *generic* arithmetic here — squaring and the error-free sum —
//! come from `euclid-numeric` rather than being reimplemented. [`sq`] stays as a
//! one-line forwarder because the formulas below are a line-by-line port and read
//! against the upstream source: `sq(calp0) * self.ep2` is recognizable where
//! `calp0.square() * self.ep2` is not. Everything else in this file is
//! geodesy-specific — degree-mode quadrant handling, the ellipsoid's `e·atanh(e·x)`,
//! the astroid seed — and belongs here, not in a general numerics crate.

use euclid_numeric::Real;

/// Newton-iteration budget before falling back to bisection alone.
pub(super) const MAXIT1: u64 = 20;
/// Total iteration budget (Newton + bisection); 53 = f64 mantissa digits.
pub(super) const MAXIT2: u64 = MAXIT1 + 53 + 10;

#[inline]
pub(super) fn tiny() -> f64 {
    f64::MIN_POSITIVE.sqrt()
}

#[inline]
pub(super) fn tol0() -> f64 {
    f64::EPSILON
}

#[inline]
pub(super) fn tol1() -> f64 {
    200.0 * tol0()
}

#[inline]
pub(super) fn tol2() -> f64 {
    tol0().sqrt()
}

#[inline]
pub(super) fn tolb() -> f64 {
    tol0() * tol2()
}

#[inline]
pub(super) fn xthresh() -> f64 {
    1000.0 * tol2()
}

#[inline]
pub(super) fn sq(x: f64) -> f64 {
    x.square()
}

/// Normalize the two-vector `(x, y)` to unit length in place.
pub(super) fn norm2(x: &mut f64, y: &mut f64) {
    let r = x.hypot(*y);
    *x /= r;
    *y /= r;
}

/// Error-free transformation of `u + v`: the rounded sum and its exact error.
///
/// `euclid-numeric`'s [`Real::two_sum`] is the same Knuth transform, and additionally
/// returns a zero error term when the sum is not finite. GeographicLib's form would
/// compute `inf - inf = NaN` there and poison every subsequent correction; the
/// angular difference below feeds its error term straight back in, so that NaN would
/// propagate into the result rather than staying in a discarded low-order bit.
#[inline]
pub(super) fn sum(u: f64, v: f64) -> (f64, f64) {
    u.two_sum(v)
}

/// Degree-`n` polynomial with coefficients `p[0..=n]`, highest power first.
///
/// An adapter over euclid's [`Real::polyval`], which takes the exact coefficient
/// slice. GeographicLib's call sites pass a degree alongside a longer table and
/// let the degree pick the window, so the windowing happens here rather than at
/// ten call sites — the same shape as [`sum`] above, which adapts
/// [`Real::two_sum`].
#[inline]
pub(super) fn polyval(n: usize, p: &[f64], x: f64) -> f64 {
    x.polyval(&p[..=n])
}

/// Sine and cosine of `x` in degrees, exact at the quadrant boundaries.
///
/// Now [`Real::sin_cos_degrees`]. Kept as a name here because the series code
/// reads as a transcription of GeographicLib, where this is `sincosd`.
#[inline]
pub(super) fn sincosd(x: f64) -> (f64, f64) {
    x.sin_cos_degrees()
}

/// Round a tiny angle (degrees) so values below 1/2⁵⁷ ≈ 0.7 pm underflow to
/// exactly 0, avoiding near-singular cases for non-zero-but-negligible inputs.
pub(super) fn ang_round(x: f64) -> f64 {
    let z = 1.0 / 16.0;
    let mut y = x.abs();
    // The compiler mustn't "simplify" z - (z - y) to y.
    if y < z {
        y = z - (z - y);
    }
    if x == 0.0 {
        0.0
    } else if x < 0.0 {
        -y
    } else {
        y
    }
}

/// Remainder of `x / y` in `[-y/2, y/2]`.
pub(super) fn remainder(x: f64, y: f64) -> f64 {
    let z = if x.is_finite() { x % y } else { f64::NAN };
    let z = if x == 0.0 { x } else { z };
    if z < -y / 2.0 {
        z + y
    } else if z < y / 2.0 {
        z
    } else {
        z - y
    }
}

/// Reduce an angle (degrees) to `(-180, 180]`.
pub(super) fn ang_normalize(x: f64) -> f64 {
    let y = remainder(x, 360.0);
    if y == -180.0 { 180.0 } else { y }
}

/// Replace latitudes outside `[-90, 90]` with NaN.
pub(super) fn lat_fix(x: f64) -> f64 {
    if x.abs() > 90.0 { f64::NAN } else { x }
}

/// `y - x`, reduced accurately to `[-180, 180]`, with the round-off error of
/// the reduction as the second component.
pub(super) fn ang_diff(x: f64, y: f64) -> (f64, f64) {
    let (d, t) = sum(ang_normalize(-x), ang_normalize(y));
    let d = ang_normalize(d);
    if d == 180.0 && t > 0.0 {
        sum(-180.0, t)
    } else {
        sum(d, t)
    }
}

/// `atan2(y, x)` in degrees, with the quadrant handled so results are exact
/// for axis-aligned arguments.
pub(super) fn atan2d(mut y: f64, mut x: f64) -> f64 {
    let mut q = if y.abs() > x.abs() {
        std::mem::swap(&mut x, &mut y);
        2
    } else {
        0
    };
    if x < 0.0 {
        q += 1;
        x = -x;
    }
    let mut ang = y.atan2(x).to_degrees();
    match q {
        0 => {}
        1 => ang = if y >= 0.0 { 180.0 - ang } else { -180.0 - ang },
        2 => ang = 90.0 - ang,
        3 => ang += -90.0,
        _ => unreachable!(),
    }
    ang
}

/// Evaluate the trigonometric series `Σ c[k]·sin(2k·x)` (`sinp`) or
/// `Σ c[k]·cos((2k−1)·x)` by Clenshaw summation. With `sinp` the `c` array is
/// 1-indexed (slot 0 unused).
pub(super) fn sin_cos_series(sinp: bool, sinx: f64, cosx: f64, c: &[f64]) -> f64 {
    let mut k = c.len();
    let mut n: i64 = k as i64 - i64::from(sinp);
    let ar: f64 = 2.0 * (cosx - sinx) * (cosx + sinx);
    let mut y1 = 0.0;
    let mut y0: f64 = if n & 1 != 0 {
        k -= 1;
        c[k]
    } else {
        0.0
    };
    n /= 2;
    while n > 0 {
        n -= 1;
        k -= 1;
        y1 = ar * y0 - y1 + c[k];
        k -= 1;
        y0 = ar * y1 - y0 + c[k];
    }
    if sinp {
        2.0 * sinx * cosx * y0
    } else {
        cosx * (y0 - y1)
    }
}

/// Solve `k⁴ + 2k³ − (x² + y² − 1)k² − 2y²k − y² = 0` for the positive root —
/// the astroid equation governing the nearly-antipodal starting guess.
pub(super) fn astroid(x: f64, y: f64) -> f64 {
    let p = sq(x);
    let q = sq(y);
    let r = (p + q - 1.0) / 6.0;
    if !(q == 0.0 && r <= 0.0) {
        let s = p * q / 4.0;
        let r2 = sq(r);
        let r3 = r * r2;
        let disc = s * (s + 2.0 * r3);
        let mut u = r;
        if disc >= 0.0 {
            let mut t3 = s + r3;
            t3 += if t3 < 0.0 { -disc.sqrt() } else { disc.sqrt() };
            let t = t3.cbrt();
            u += t + if t != 0.0 { r2 / t } else { 0.0 };
        } else {
            let ang = (-disc).sqrt().atan2(-(s + r3));
            u += 2.0 * r * (ang / 3.0).cos();
        }
        let v = (sq(u) + q).sqrt();
        let uv = if u < 0.0 { q / (v - u) } else { u + v };
        let w = (uv - q) / (2.0 * v);
        uv / ((uv + sq(w)).sqrt() + w)
    } else {
        0.0
    }
}

/// `+1`/`-1` if the edge `lon1 → lon2` crosses the prime meridian going east /
/// west (`±0` counts as positive longitude), else `0` — the pole-encirclement
/// counter of the polygon-area accumulator.
pub(super) fn transit(lon1: f64, lon2: f64) -> i64 {
    let (lon12, _) = ang_diff(lon1, lon2);
    let lon1 = ang_normalize(lon1);
    let lon2 = ang_normalize(lon2);
    if lon12 > 0.0 && ((lon1 < 0.0 && lon2 >= 0.0) || (lon1 > 0.0 && lon2 == 0.0)) {
        1
    } else if lon12 < 0.0 && lon1 >= 0.0 && lon2 < 0.0 {
        -1
    } else {
        0
    }
}

/// `e·atanh(e·x)` for an ellipsoid with signed `es = sign(f)·e` — handles both
/// oblate (`atanh`) and prolate (`atan`) cases.
pub(super) fn eatanhe(x: f64, es: f64) -> f64 {
    if es > 0.0 {
        es * (es * x).atanh()
    } else {
        -es * (es * x).atan()
    }
}
