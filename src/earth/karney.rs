//! Ellipsoidal geodesics — Karney's algorithms for the direct and inverse
//! geodesic problems on an ellipsoid of revolution, accurate to round-off
//! (~15 nm on Earth) for any pair of points, including near-antipodal ones
//! where Vincenty's method fails to converge.
//!
//! Ported from `geographiclib-rs` 0.2.7 (MIT), itself a port of Charles F. F.
//! Karney's GeographicLib (see the crate `NOTICE`); the algorithms are from
//! Karney, *Algorithms for geodesics*, J. Geodesy 87, 43–55 (2013). The port
//! drops the capability masks, arc-mode positioning, and the reduced-length /
//! geodesic-scale / area outputs — [`Karney::inverse`] and [`Karney::direct`]
//! always solve for distance and azimuths, nothing else. `sincosd` uses
//! GeographicLib's `fmod`-based quadrant reduction (the Python/JS form) rather
//! than `remquo`, which hausdorff has no dependency for.
//!
//! The public surface is on [`Ellipsoid`](crate::earth::Ellipsoid)
//! (`distance` / `inverse` / `destination` / `direct`); this module is the
//! kernel. Differential-tested against `geographiclib-rs` (see
//! `tests/geodesic.rs`).

use std::f64::consts::PI;
use std::sync::LazyLock;

use crate::earth::ellipsoid::Ellipsoid;

/// Order of the series expansions (GeographicLib's compile-time default).
const GEODESIC_ORDER: usize = 6;
/// The `C1`/`C1'`/`C2` Fourier coefficient arrays are 1-indexed: slot 0 unused.
const CARR_SIZE: usize = GEODESIC_ORDER + 1;
/// Number of entries in the `C3` coefficient triangle.
const NC3X: usize = 15;

/// Newton-iteration budget before falling back to bisection alone.
const MAXIT1: u64 = 20;
/// Total iteration budget (Newton + bisection); 53 = f64 mantissa digits.
const MAXIT2: u64 = MAXIT1 + 53 + 10;

#[inline]
fn tiny() -> f64 {
    f64::MIN_POSITIVE.sqrt()
}

#[inline]
fn tol0() -> f64 {
    f64::EPSILON
}

#[inline]
fn tol1() -> f64 {
    200.0 * tol0()
}

#[inline]
fn tol2() -> f64 {
    tol0().sqrt()
}

#[inline]
fn tolb() -> f64 {
    tol0() * tol2()
}

#[inline]
fn xthresh() -> f64 {
    1000.0 * tol2()
}

#[inline]
fn sq(x: f64) -> f64 {
    x * x
}

/// Normalize the two-vector `(x, y)` to unit length in place.
fn norm2(x: &mut f64, y: &mut f64) {
    let r = x.hypot(*y);
    *x /= r;
    *y /= r;
}

/// Error-free transformation of `u + v`: the rounded sum and its exact error.
fn sum(u: f64, v: f64) -> (f64, f64) {
    let s = u + v;
    let up = s - v;
    let vpp = s - up;
    let up = up - u;
    let vpp = vpp - v;
    (s, -(up + vpp))
}

/// Evaluate the degree-`n` polynomial with coefficients `p[0..=n]` (highest
/// power first) at `x` by Horner's rule.
fn polyval(n: usize, p: &[f64], x: f64) -> f64 {
    let mut y = p[0];
    for val in &p[1..=n] {
        y = y * x + val;
    }
    y
}

/// Round a tiny angle (degrees) so values below 1/2⁵⁷ ≈ 0.7 pm underflow to
/// exactly 0, avoiding near-singular cases for non-zero-but-negligible inputs.
fn ang_round(x: f64) -> f64 {
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
fn remainder(x: f64, y: f64) -> f64 {
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
fn ang_normalize(x: f64) -> f64 {
    let y = remainder(x, 360.0);
    if y == -180.0 { 180.0 } else { y }
}

/// Replace latitudes outside `[-90, 90]` with NaN.
fn lat_fix(x: f64) -> f64 {
    if x.abs() > 90.0 { f64::NAN } else { x }
}

/// `y - x`, reduced accurately to `[-180, 180]`, with the round-off error of
/// the reduction as the second component.
fn ang_diff(x: f64, y: f64) -> (f64, f64) {
    let (d, t) = sum(ang_normalize(-x), ang_normalize(y));
    let d = ang_normalize(d);
    if d == 180.0 && t > 0.0 {
        sum(-180.0, t)
    } else {
        sum(d, t)
    }
}

/// Sine and cosine of `x` in degrees, exact at the quadrant boundaries
/// (GeographicLib's `fmod`-based reduction).
fn sincosd(x: f64) -> (f64, f64) {
    let r = if x.is_finite() { x % 360.0 } else { f64::NAN };
    let q = if r.is_nan() {
        0
    } else {
        (r / 90.0).round() as i64
    };
    let r = (r - 90.0 * q as f64).to_radians();
    let (s, c) = r.sin_cos();
    let (mut sinx, cosx) = match q.rem_euclid(4) {
        0 => (s, c),
        1 => (c, -s),
        2 => (-s, -c),
        3 => (-c, s),
        _ => unreachable!(),
    };
    // IEEE 754 special values (F.10.1.12/13): -0 cosines become +0; a zero
    // sine keeps the sign of the argument.
    let cosx = cosx + 0.0;
    if sinx == 0.0 {
        sinx = sinx.copysign(x);
    }
    (sinx, cosx)
}

/// `atan2(y, x)` in degrees, with the quadrant handled so results are exact
/// for axis-aligned arguments.
fn atan2d(mut y: f64, mut x: f64) -> f64 {
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
fn sin_cos_series(sinp: bool, sinx: f64, cosx: f64, c: &[f64]) -> f64 {
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
fn astroid(x: f64, y: f64) -> f64 {
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

/// `A₁(ε) − 1`: the mean-to-true distance scale of the series.
fn a1m1f(eps: f64) -> f64 {
    const COEFF: [f64; 5] = [1.0, 4.0, 64.0, 0.0, 256.0];
    let m = GEODESIC_ORDER / 2;
    let t = polyval(m, &COEFF, sq(eps)) / COEFF[m + 1];
    (t + eps) / (1.0 - eps)
}

/// The `C₁` Fourier coefficients (distance series).
fn c1f(eps: f64, c: &mut [f64; CARR_SIZE]) {
    const COEFF: [f64; 18] = [
        -1.0, 6.0, -16.0, 32.0, -9.0, 64.0, -128.0, 2048.0, 9.0, -16.0, 768.0, 3.0, -5.0, 512.0,
        -7.0, 1280.0, -7.0, 2048.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    // Indexed loop as in the reference port (the iterator form benched slower
    // there).
    #[allow(clippy::needless_range_loop)]
    for l in 1..=GEODESIC_ORDER {
        let m = (GEODESIC_ORDER - l) / 2;
        c[l] = d * polyval(m, &COEFF[o..], eps2) / COEFF[o + m + 1];
        o += m + 2;
        d *= eps;
    }
}

/// The `C₁'` Fourier coefficients (the inverse of the distance series, for the
/// direct problem).
fn c1pf(eps: f64, c: &mut [f64; CARR_SIZE]) {
    const COEFF: [f64; 18] = [
        205.0, -432.0, 768.0, 1536.0, 4005.0, -4736.0, 3840.0, 12288.0, -225.0, 116.0, 384.0,
        -7173.0, 2695.0, 7680.0, 3467.0, 7680.0, 38081.0, 61440.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    // Indexed loop as in the reference port (the iterator form benched slower
    // there).
    #[allow(clippy::needless_range_loop)]
    for l in 1..=GEODESIC_ORDER {
        let m = (GEODESIC_ORDER - l) / 2;
        c[l] = d * polyval(m, &COEFF[o..], eps2) / COEFF[o + m + 1];
        o += m + 2;
        d *= eps;
    }
}

/// `A₂(ε) − 1`: the scale of the reduced-length series.
fn a2m1f(eps: f64) -> f64 {
    const COEFF: [f64; 5] = [-11.0, -28.0, -192.0, 0.0, 256.0];
    let m = GEODESIC_ORDER / 2;
    let t = polyval(m, &COEFF, sq(eps)) / COEFF[m + 1];
    (t - eps) / (1.0 + eps)
}

/// The `C₂` Fourier coefficients (reduced-length series).
fn c2f(eps: f64, c: &mut [f64; CARR_SIZE]) {
    const COEFF: [f64; 18] = [
        1.0, 2.0, 16.0, 32.0, 35.0, 64.0, 384.0, 2048.0, 15.0, 80.0, 768.0, 7.0, 35.0, 512.0, 63.0,
        1280.0, 77.0, 2048.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    // Indexed loop as in the reference port (the iterator form benched slower
    // there).
    #[allow(clippy::needless_range_loop)]
    for l in 1..=GEODESIC_ORDER {
        let m = (GEODESIC_ORDER - l) / 2;
        c[l] = d * polyval(m, &COEFF[o..], eps2) / COEFF[o + m + 1];
        o += m + 2;
        d *= eps;
    }
}

/// Coefficients of the `A₃` polynomial in the third flattening `n`.
const COEFF_A3: [f64; 18] = [
    -3.0, 128.0, -2.0, -3.0, 64.0, -1.0, -3.0, -1.0, 16.0, 3.0, -1.0, -2.0, 8.0, 1.0, -1.0, 2.0,
    1.0, 1.0,
];

/// Coefficients of the `C₃` polynomials in the third flattening `n`.
const COEFF_C3: [f64; 45] = [
    3.0, 128.0, 2.0, 5.0, 128.0, -1.0, 3.0, 3.0, 64.0, -1.0, 0.0, 1.0, 8.0, -1.0, 1.0, 4.0, 5.0,
    256.0, 1.0, 3.0, 128.0, -3.0, -2.0, 3.0, 64.0, 1.0, -3.0, 2.0, 32.0, 7.0, 512.0, -10.0, 9.0,
    384.0, 5.0, -9.0, 5.0, 192.0, 7.0, 512.0, -14.0, 7.0, 512.0, 21.0, 2560.0,
];

/// The geodesic solver for one ellipsoid: the ellipsoid's derived constants
/// plus the `A₃`/`C₃` series coefficients, which depend only on the third
/// flattening and so are computed once per ellipsoid.
pub(crate) struct Karney {
    a: f64,
    f: f64,
    f1: f64,
    ep2: f64,
    n: f64,
    b: f64,
    etol2: f64,
    a3x: [f64; GEODESIC_ORDER],
    c3x: [f64; NC3X],
}

/// The WGS84 solver, built once — the common case ([`Ellipsoid::WGS84`] and
/// the `space::r2` `Geodesic` metric) skips re-deriving the series
/// coefficients per call.
static WGS84: LazyLock<Karney> = LazyLock::new(|| Karney::new(Ellipsoid::WGS84));

/// The cached WGS84 solver.
pub(crate) fn wgs84() -> &'static Karney {
    &WGS84
}

impl Karney {
    pub(crate) fn new(ellipsoid: Ellipsoid) -> Self {
        let a = ellipsoid.a;
        let f = ellipsoid.f;
        let f1 = 1.0 - f;
        let ep2 = ellipsoid.ep2();
        let n = ellipsoid.n();
        let b = ellipsoid.b();
        let etol2 = 0.1 * tol2() / (f.abs().max(0.001) * (1.0 - f / 2.0).min(1.0) / 2.0).sqrt();

        let mut a3x = [0.0; GEODESIC_ORDER];
        let mut o = 0;
        for (k, j) in (0..GEODESIC_ORDER).rev().enumerate() {
            let m = j.min(GEODESIC_ORDER - j - 1);
            a3x[k] = polyval(m, &COEFF_A3[o..], n) / COEFF_A3[o + m + 1];
            o += m + 2;
        }

        let mut c3x = [0.0; NC3X];
        let mut o = 0;
        let mut k = 0;
        for l in 1..GEODESIC_ORDER {
            for j in (l..GEODESIC_ORDER).rev() {
                let m = j.min(GEODESIC_ORDER - j - 1);
                c3x[k] = polyval(m, &COEFF_C3[o..], n) / COEFF_C3[o + m + 1];
                k += 1;
                o += m + 2;
            }
        }

        Karney {
            a,
            f,
            f1,
            ep2,
            n,
            b,
            etol2,
            a3x,
            c3x,
        }
    }

    /// `A₃(ε)`: the longitude-series scale.
    fn a3f(&self, eps: f64) -> f64 {
        polyval(GEODESIC_ORDER - 1, &self.a3x, eps)
    }

    /// The `C₃` Fourier coefficients (longitude series) at `eps`.
    fn c3f(&self, eps: f64, c: &mut [f64; GEODESIC_ORDER]) {
        let mut mult = 1.0;
        let mut o = 0;
        #[allow(clippy::needless_range_loop)]
        for l in 1..GEODESIC_ORDER {
            let m = GEODESIC_ORDER - l - 1;
            mult *= eps;
            c[l] = mult * polyval(m, &self.c3x[o..], eps);
            o += m + 1;
        }
    }

    /// The length integrals between two points on a geodesic, in units of `b`:
    /// the distance `s12b` (if `want_distance`) and the reduced length `m12b`
    /// with the series constant `m0` (if `want_reduced`); unsolicited outputs
    /// are NaN. The bool pair replaces GeographicLib's capability mask.
    #[allow(clippy::too_many_arguments)]
    fn lengths(
        &self,
        eps: f64,
        sig12: f64,
        ssig1: f64,
        csig1: f64,
        dn1: f64,
        ssig2: f64,
        csig2: f64,
        dn2: f64,
        want_distance: bool,
        want_reduced: bool,
        c1a: &mut [f64; CARR_SIZE],
        c2a: &mut [f64; CARR_SIZE],
    ) -> (f64, f64, f64) {
        let mut s12b = f64::NAN;
        let mut m12b = f64::NAN;
        let mut m0 = f64::NAN;

        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut m0x = 0.0;
        let mut j12 = 0.0;
        if want_distance || want_reduced {
            a1 = a1m1f(eps);
            c1f(eps, c1a);
            if want_reduced {
                a2 = a2m1f(eps);
                c2f(eps, c2a);
                m0x = a1 - a2;
                a2 += 1.0;
            }
            a1 += 1.0;
        }
        if want_distance {
            let b1 =
                sin_cos_series(true, ssig2, csig2, c1a) - sin_cos_series(true, ssig1, csig1, c1a);
            s12b = a1 * (sig12 + b1);
            if want_reduced {
                let b2 = sin_cos_series(true, ssig2, csig2, c2a)
                    - sin_cos_series(true, ssig1, csig1, c2a);
                j12 = m0x * sig12 + (a1 * b1 - a2 * b2);
            }
        } else if want_reduced {
            for l in 1..=GEODESIC_ORDER {
                c2a[l] = a1 * c1a[l] - a2 * c2a[l];
            }
            j12 = m0x * sig12
                + (sin_cos_series(true, ssig2, csig2, c2a)
                    - sin_cos_series(true, ssig1, csig1, c2a));
        }
        if want_reduced {
            m0 = m0x;
            m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * j12;
        }
        (s12b, m12b, m0)
    }

    /// The starting guess for the inverse problem's Newton iteration: returns
    /// `(sig12, salp1, calp1, salp2, calp2, dnm)`. `sig12 ≥ 0` means the
    /// short-line case was solved exactly and no iteration is needed; the
    /// nearly-antipodal case seeds `salp1`/`calp1` from the astroid equation.
    #[allow(clippy::too_many_arguments)]
    fn inverse_start(
        &self,
        sbet1: f64,
        cbet1: f64,
        dn1: f64,
        sbet2: f64,
        cbet2: f64,
        dn2: f64,
        lam12: f64,
        slam12: f64,
        clam12: f64,
        c1a: &mut [f64; CARR_SIZE],
        c2a: &mut [f64; CARR_SIZE],
    ) -> (f64, f64, f64, f64, f64, f64) {
        let mut sig12 = -1.0;
        let mut salp2 = f64::NAN;
        let mut calp2 = f64::NAN;
        let mut dnm = f64::NAN;

        let mut somg12: f64;
        let mut comg12: f64;

        let sbet12 = sbet2 * cbet1 - cbet2 * sbet1;
        let cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
        let sbet12a = sbet2 * cbet1 + cbet2 * sbet1;

        let shortline = cbet12 >= 0.0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
        if shortline {
            let mut sbetm2 = sq(sbet1 + sbet2);
            sbetm2 /= sbetm2 + sq(cbet1 + cbet2);
            dnm = (1.0 + self.ep2 * sbetm2).sqrt();
            let omg12 = lam12 / (self.f1 * dnm);
            somg12 = omg12.sin();
            comg12 = omg12.cos();
        } else {
            somg12 = slam12;
            comg12 = clam12;
        }

        let mut salp1 = cbet2 * somg12;
        let mut calp1 = if comg12 >= 0.0 {
            sbet12 + cbet2 * sbet1 * sq(somg12) / (1.0 + comg12)
        } else {
            sbet12a - cbet2 * sbet1 * sq(somg12) / (1.0 - comg12)
        };

        let ssig12 = salp1.hypot(calp1);
        let csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

        if shortline && ssig12 < self.etol2 {
            salp2 = cbet1 * somg12;
            calp2 = sbet12
                - cbet1
                    * sbet2
                    * (if comg12 >= 0.0 {
                        sq(somg12) / (1.0 + comg12)
                    } else {
                        1.0 - comg12
                    });
            norm2(&mut salp2, &mut calp2);
            sig12 = ssig12.atan2(csig12);
        } else if self.n.abs() > 0.1
            || csig12 >= 0.0
            || ssig12 >= 6.0 * self.n.abs() * PI * sq(cbet1)
        {
            // The default starting guess (salp1/calp1 above) is good enough.
        } else {
            // Nearly antipodal: solve the astroid problem for the guess.
            let x: f64;
            let y: f64;
            let betscale: f64;
            let lamscale: f64;
            let lam12x = (-slam12).atan2(-clam12);
            if self.f >= 0.0 {
                let k2 = sq(sbet1) * self.ep2;
                let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
                lamscale = self.f * cbet1 * self.a3f(eps) * PI;
                betscale = lamscale * cbet1;
                x = lam12x / lamscale;
                y = sbet12a / betscale;
            } else {
                let cbet12a = cbet2 * cbet1 - sbet2 * sbet1;
                let bet12a = sbet12a.atan2(cbet12a);
                let (_, m12b, m0) = self.lengths(
                    self.n,
                    PI + bet12a,
                    sbet1,
                    -cbet1,
                    dn1,
                    sbet2,
                    cbet2,
                    dn2,
                    false,
                    true,
                    c1a,
                    c2a,
                );
                x = -1.0 + m12b / (cbet1 * cbet2 * m0 * PI);
                betscale = if x < -0.01 {
                    sbet12a / x
                } else {
                    -self.f * sq(cbet1) * PI
                };
                lamscale = betscale / cbet1;
                y = lam12x / lamscale;
            }
            if y > -tol1() && x > -1.0 - xthresh() {
                if self.f >= 0.0 {
                    salp1 = (-x).min(1.0);
                    calp1 = -(1.0 - sq(salp1)).sqrt();
                } else {
                    calp1 = x.max(if x > -tol1() { 0.0 } else { -1.0 });
                    salp1 = (1.0 - sq(calp1)).sqrt();
                }
            } else {
                let k = astroid(x, y);
                let omg12a = lamscale
                    * if self.f >= 0.0 {
                        -x * k / (1.0 + k)
                    } else {
                        -y * (1.0 + k) / k
                    };
                somg12 = omg12a.sin();
                comg12 = -(omg12a.cos());
                salp1 = cbet2 * somg12;
                calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1.0 - comg12);
            }
        }

        if salp1 > 0.0 || salp1.is_nan() {
            norm2(&mut salp1, &mut calp1);
        } else {
            salp1 = 1.0;
            calp1 = 0.0;
        }
        (sig12, salp1, calp1, salp2, calp2, dnm)
    }

    /// One evaluation of the inverse problem for a trial `alp1`: the longitude
    /// mismatch `v = λ₁₂(alp1) − λ₁₂` (the Newton residual) and the geometry
    /// at that azimuth. Returns
    /// `(v, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, dv)`, where
    /// `dv` is `dv/dalp1` when `diffp` (else NaN).
    #[allow(clippy::too_many_arguments)]
    fn lambda12(
        &self,
        sbet1: f64,
        cbet1: f64,
        dn1: f64,
        sbet2: f64,
        cbet2: f64,
        dn2: f64,
        salp1: f64,
        mut calp1: f64,
        slam120: f64,
        clam120: f64,
        diffp: bool,
        c1a: &mut [f64; CARR_SIZE],
        c2a: &mut [f64; CARR_SIZE],
        c3a: &mut [f64; GEODESIC_ORDER],
    ) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        if sbet1 == 0.0 && calp1 == 0.0 {
            // Break the degeneracy of an equatorial line: any azimuth works.
            calp1 = -tiny();
        }
        let salp0 = salp1 * cbet1;
        let calp0 = calp1.hypot(salp1 * sbet1);

        let mut ssig1 = sbet1;
        let somg1 = salp0 * sbet1;
        let mut csig1 = calp1 * cbet1;
        let comg1 = calp1 * cbet1;
        norm2(&mut ssig1, &mut csig1);

        let salp2 = if cbet2 != cbet1 { salp0 / cbet2 } else { salp1 };
        let calp2 = if cbet2 != cbet1 || sbet2.abs() != -sbet1 {
            (sq(calp1 * cbet1)
                + if cbet1 < -sbet1 {
                    (cbet2 - cbet1) * (cbet1 + cbet2)
                } else {
                    (sbet1 - sbet2) * (sbet1 + sbet2)
                })
            .sqrt()
                / cbet2
        } else {
            calp1.abs()
        };
        let mut ssig2 = sbet2;
        let somg2 = salp0 * sbet2;
        let mut csig2 = calp2 * cbet2;
        let comg2 = calp2 * cbet2;
        norm2(&mut ssig2, &mut csig2);

        let sig12 = ((csig1 * ssig2 - ssig1 * csig2).max(0.0)).atan2(csig1 * csig2 + ssig1 * ssig2);
        let somg12 = (comg1 * somg2 - somg1 * comg2).max(0.0);
        let comg12 = comg1 * comg2 + somg1 * somg2;
        let eta = (somg12 * clam120 - comg12 * slam120).atan2(comg12 * clam120 + somg12 * slam120);

        let k2 = sq(calp0) * self.ep2;
        let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
        self.c3f(eps, c3a);
        let b312 =
            sin_cos_series(true, ssig2, csig2, c3a) - sin_cos_series(true, ssig1, csig1, c3a);
        let domg12 = -self.f * self.a3f(eps) * salp0 * (sig12 + b312);
        let lam12 = eta + domg12;

        let dlam12 = if diffp {
            if calp2 == 0.0 {
                -2.0 * self.f1 * dn1 / sbet1
            } else {
                let (_, m12b, _) = self.lengths(
                    eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, false, true, c1a, c2a,
                );
                m12b * self.f1 / (calp2 * cbet2)
            }
        } else {
            f64::NAN
        };
        (
            lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, dlam12,
        )
    }

    /// The inverse geodesic problem: the shortest-path distance in metres and
    /// the forward azimuths (degrees, clockwise from north) at both endpoints,
    /// for two points in latitude/longitude degrees. Returns
    /// `(s12, azi1, azi2)`; NaN for a non-finite or out-of-range latitude.
    pub(crate) fn inverse(&self, lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> (f64, f64, f64) {
        let (mut lon12, mut lon12s) = ang_diff(lon1, lon2);
        let mut lonsign = if lon12 >= 0.0 { 1.0 } else { -1.0 };

        // Canonicalize: east-going (lonsign), |lat1| ≥ |lat2| (swapp),
        // lat1 ≤ 0 (latsign); undone on the azimuths at the end.
        lon12 = lonsign * ang_round(lon12);
        lon12s = ang_round((180.0 - lon12) - lonsign * lon12s);
        let lam12 = lon12.to_radians();
        let (slam12, clam12) = if lon12 > 90.0 {
            let (s, c) = sincosd(lon12s);
            (s, -c)
        } else {
            sincosd(lon12)
        };
        let mut lat1 = ang_round(lat_fix(lat1));
        let mut lat2 = ang_round(lat_fix(lat2));

        let swapp = if lat1.abs() < lat2.abs() { -1.0 } else { 1.0 };
        if swapp < 0.0 {
            lonsign *= -1.0;
            std::mem::swap(&mut lat1, &mut lat2);
        }
        let latsign = if lat1 < 0.0 { 1.0 } else { -1.0 };
        lat1 *= latsign;
        lat2 *= latsign;

        let (mut sbet1, mut cbet1) = sincosd(lat1);
        sbet1 *= self.f1;
        norm2(&mut sbet1, &mut cbet1);
        cbet1 = cbet1.max(tiny());

        let (mut sbet2, mut cbet2) = sincosd(lat2);
        sbet2 *= self.f1;
        norm2(&mut sbet2, &mut cbet2);
        cbet2 = cbet2.max(tiny());

        // Tighten near-symmetric configurations so the strict inequalities the
        // algorithm relies on hold exactly.
        if cbet1 < -sbet1 {
            if cbet2 == cbet1 {
                sbet2 = if sbet2 < 0.0 { sbet1 } else { -sbet1 };
            }
        } else if sbet2.abs() == -sbet1 {
            cbet2 = cbet1;
        }

        let dn1 = (1.0 + self.ep2 * sq(sbet1)).sqrt();
        let dn2 = (1.0 + self.ep2 * sq(sbet2)).sqrt();

        let mut c1a = [0.0; CARR_SIZE];
        let mut c2a = [0.0; CARR_SIZE];
        let mut c3a = [0.0; GEODESIC_ORDER];

        let mut meridian = lat1 == -90.0 || slam12 == 0.0;
        let mut calp1 = 0.0;
        let mut salp1 = 0.0;
        let mut calp2 = 0.0;
        let mut salp2 = 0.0;
        let mut ssig1 = 0.0;
        let mut csig1 = 0.0;
        let mut ssig2 = 0.0;
        let mut csig2 = 0.0;
        let mut s12x = 0.0;

        if meridian {
            // The geodesic runs along the meridian (possibly over a pole).
            calp1 = clam12;
            salp1 = slam12;
            calp2 = 1.0;
            salp2 = 0.0;

            ssig1 = sbet1;
            csig1 = calp1 * cbet1;
            ssig2 = sbet2;
            csig2 = calp2 * cbet2;

            let sig12 =
                ((csig1 * ssig2 - ssig1 * csig2).max(0.0)).atan2(csig1 * csig2 + ssig1 * ssig2);
            let (s12b, m12b, _) = self.lengths(
                self.n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, true, true, &mut c1a, &mut c2a,
            );
            s12x = s12b;

            if sig12 < 1.0 || m12b >= 0.0 {
                if sig12 < 3.0 * tiny() {
                    s12x = 0.0;
                }
                s12x *= self.b;
            } else {
                // The reduced length went negative: the meridian path is no
                // longer the shortest; fall through to the general solver.
                meridian = false;
            }
        }

        if !meridian && sbet1 == 0.0 && (self.f <= 0.0 || lon12s >= self.f * 180.0) {
            // Both points on the equator, geodesic stays on it.
            calp1 = 0.0;
            calp2 = 0.0;
            salp1 = 1.0;
            salp2 = 1.0;
            s12x = self.a * lam12;
        } else if !meridian {
            let (sig12i, salp1i, calp1i, salp2i, calp2i, dnm) = self.inverse_start(
                sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, &mut c1a, &mut c2a,
            );
            let mut sig12 = sig12i;
            salp1 = salp1i;
            calp1 = calp1i;
            salp2 = salp2i;
            calp2 = calp2i;

            if sig12 >= 0.0 {
                // Short line: the starting guess is the answer.
                s12x = sig12 * self.b * dnm;
            } else {
                // Newton's method on alp1, bracketed by bisection.
                let mut tripn = false;
                let mut tripb = false;
                let mut salp1a = tiny();
                let mut calp1a = 1.0;
                let mut salp1b = tiny();
                let mut calp1b = -1.0;
                let mut eps = 0.0;
                for numit in 0..MAXIT2 {
                    let (v, salp2i, calp2i, sig12i, ssig1i, csig1i, ssig2i, csig2i, epsi, dv) =
                        self.lambda12(
                            sbet1,
                            cbet1,
                            dn1,
                            sbet2,
                            cbet2,
                            dn2,
                            salp1,
                            calp1,
                            slam12,
                            clam12,
                            numit < MAXIT1,
                            &mut c1a,
                            &mut c2a,
                            &mut c3a,
                        );
                    salp2 = salp2i;
                    calp2 = calp2i;
                    sig12 = sig12i;
                    ssig1 = ssig1i;
                    csig1 = csig1i;
                    ssig2 = ssig2i;
                    csig2 = csig2i;
                    eps = epsi;

                    if tripb || v.abs() < (if tripn { 8.0 } else { 1.0 }) * tol0() || v.is_nan() {
                        break;
                    }
                    // Maintain the bracket [alp1a, alp1b] around the root.
                    if v > 0.0 && (numit > MAXIT1 || calp1 / salp1 > calp1b / salp1b) {
                        salp1b = salp1;
                        calp1b = calp1;
                    } else if v < 0.0 && (numit > MAXIT1 || calp1 / salp1 < calp1a / salp1a) {
                        salp1a = salp1;
                        calp1a = calp1;
                    }
                    if numit < MAXIT1 && dv > 0.0 {
                        let dalp1 = -v / dv;
                        let (sdalp1, cdalp1) = dalp1.sin_cos();
                        let nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
                        if nsalp1 > 0.0 && dalp1.abs() < PI {
                            calp1 = calp1 * cdalp1 - salp1 * sdalp1;
                            salp1 = nsalp1;
                            norm2(&mut salp1, &mut calp1);
                            // Near convergence Newton may oscillate at
                            // round-off; loosen the exit threshold.
                            tripn = v.abs() <= 16.0 * tol0();
                            continue;
                        }
                    }
                    // Newton stepped outside the bracket (or its derivative is
                    // unusable): bisect instead.
                    salp1 = (salp1a + salp1b) / 2.0;
                    calp1 = (calp1a + calp1b) / 2.0;
                    norm2(&mut salp1, &mut calp1);
                    tripn = false;
                    tripb = (salp1a - salp1).abs() + (calp1a - calp1) < tolb()
                        || (salp1 - salp1b).abs() + (calp1 - calp1b) < tolb();
                }
                let (s12b, _, _) = self.lengths(
                    eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, true, false, &mut c1a,
                    &mut c2a,
                );
                s12x = s12b * self.b;
            }
        }

        // Undo the canonical swaps/sign flips on the azimuths.
        if swapp < 0.0 {
            std::mem::swap(&mut salp2, &mut salp1);
            std::mem::swap(&mut calp2, &mut calp1);
        }
        salp1 *= swapp * lonsign;
        calp1 *= swapp * latsign;
        salp2 *= swapp * lonsign;
        calp2 *= swapp * latsign;

        (s12x, atan2d(salp1, calp1), atan2d(salp2, calp2))
    }

    /// The direct geodesic problem: from a point in latitude/longitude
    /// degrees, travel `s12` metres along the geodesic with initial azimuth
    /// `azi1` (degrees, clockwise from north). Returns `(lat2, lon2, azi2)`
    /// with `azi2` the forward azimuth at the destination; NaN for a
    /// non-finite or out-of-range latitude.
    pub(crate) fn direct(&self, lat1: f64, lon1: f64, azi1: f64, s12: f64) -> (f64, f64, f64) {
        let azi1 = ang_normalize(azi1);
        let (salp1, calp1) = sincosd(ang_round(azi1));
        let lat1 = lat_fix(lat1);

        let (mut sbet1, mut cbet1) = sincosd(ang_round(lat1));
        sbet1 *= self.f1;
        norm2(&mut sbet1, &mut cbet1);
        cbet1 = cbet1.max(tiny());
        let salp0 = salp1 * cbet1;
        let calp0 = calp1.hypot(salp1 * sbet1);

        let mut ssig1 = sbet1;
        let somg1 = salp0 * sbet1;
        let mut csig1 = if sbet1 != 0.0 || calp1 != 0.0 {
            cbet1 * calp1
        } else {
            1.0
        };
        let comg1 = csig1;
        norm2(&mut ssig1, &mut csig1);

        let k2 = sq(calp0) * self.ep2;
        let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);

        let a1m1 = a1m1f(eps);
        let mut c1a = [0.0; CARR_SIZE];
        c1f(eps, &mut c1a);
        let b11 = sin_cos_series(true, ssig1, csig1, &c1a);
        let (s, c) = b11.sin_cos();
        // tau1 = sig1 + B11: the "mean" arc the distance series inverts over.
        let stau1 = ssig1 * c + csig1 * s;
        let ctau1 = csig1 * c - ssig1 * s;

        let mut c1pa = [0.0; CARR_SIZE];
        c1pf(eps, &mut c1pa);

        let mut c3a = [0.0; GEODESIC_ORDER];
        self.c3f(eps, &mut c3a);
        let a3c = -self.f * salp0 * self.a3f(eps);
        let b31 = sin_cos_series(true, ssig1, csig1, &c3a);

        // Invert the distance series: s12 → sig12, via tau and the C1' series.
        let tau12 = s12 / (self.b * (1.0 + a1m1));
        let (s, c) = tau12.sin_cos();
        let b12 = -sin_cos_series(true, stau1 * c + ctau1 * s, ctau1 * c - stau1 * s, &c1pa);
        let mut sig12 = tau12 - (b12 - b11);
        let mut ssig12 = sig12.sin();
        let mut csig12 = sig12.cos();
        if self.f.abs() > 0.01 {
            // One Newton correction recovers full accuracy for |f| > 1/100.
            let ssig2 = ssig1 * csig12 + csig1 * ssig12;
            let csig2 = csig1 * csig12 - ssig1 * ssig12;
            let b12 = sin_cos_series(true, ssig2, csig2, &c1a);
            let serr = (1.0 + a1m1) * (sig12 + (b12 - b11)) - s12 / self.b;
            sig12 -= serr / (1.0 + k2 * sq(ssig2)).sqrt();
            ssig12 = sig12.sin();
            csig12 = sig12.cos();
        }

        let ssig2 = ssig1 * csig12 + csig1 * ssig12;
        let mut csig2 = csig1 * csig12 - ssig1 * ssig12;
        let sbet2 = calp0 * ssig2;
        let mut cbet2 = salp0.hypot(calp0 * csig2);
        if cbet2 == 0.0 {
            // Over a pole: nudge off the singularity.
            cbet2 = tiny();
            csig2 = tiny();
        }
        let salp2 = salp0;
        let calp2 = calp0 * csig2;

        let somg2 = salp0 * ssig2;
        let comg2 = csig2;
        let omg12 = (somg2 * comg1 - comg2 * somg1).atan2(comg2 * comg1 + somg2 * somg1);
        let lam12 = omg12 + a3c * (sig12 + (sin_cos_series(true, ssig2, csig2, &c3a) - b31));
        let lon12 = lam12.to_degrees();
        let lon2 = ang_normalize(ang_normalize(lon1) + ang_normalize(lon12));
        let lat2 = atan2d(sbet2, self.f1 * cbet2);
        let azi2 = atan2d(salp2, calp2);
        (lat2, lon2, azi2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sincosd_is_exact_at_quadrants() {
        assert_eq!(sincosd(0.0), (0.0, 1.0));
        assert_eq!(sincosd(90.0), (1.0, 0.0));
        assert_eq!(sincosd(180.0), (0.0, -1.0));
        assert_eq!(sincosd(-90.0), (-1.0, 0.0));
        assert_eq!(sincosd(720.0), (0.0, 1.0));
    }

    #[test]
    fn equator_arc_is_exact() {
        // Along the equator the geodesic is the equator itself: one degree of
        // longitude is exactly a/√(1) scaled by the radian measure.
        let k = Karney::new(Ellipsoid::WGS84);
        let (s12, azi1, azi2) = k.inverse(0.0, 0.0, 0.0, 1.0);
        assert!((azi1 - 90.0).abs() < 1e-13 && (azi2 - 90.0).abs() < 1e-13);
        let expected = Ellipsoid::WGS84.a * 1f64.to_radians();
        assert!((s12 - expected).abs() < 1e-6, "{s12} vs {expected}");
    }

    #[test]
    fn sphere_matches_haversine() {
        // With f = 0 the Karney machinery degenerates to great circles.
        let k = Karney::new(Ellipsoid::SPHERE);
        let a = crate::earth::GeoPoint::from_lat_lon(48.85, 2.35);
        let b = crate::earth::GeoPoint::from_lat_lon(-33.87, 151.21);
        let (s12, _, _) = k.inverse(a.lat, a.lon, b.lat, b.lon);
        let hav = crate::earth::geodesic::haversine_on_m(a, b, Ellipsoid::SPHERE.a);
        assert!((s12 - hav).abs() < 1e-6, "{s12} vs {hav}");
    }

    #[test]
    fn direct_inverts_inverse() {
        let k = wgs84();
        let (lat1, lon1, lat2, lon2) = (40.6398, -73.7787, 51.4775, -0.4614); // JFK → LHR
        let (s12, azi1, _) = k.inverse(lat1, lon1, lat2, lon2);
        let (dlat, dlon, _) = k.direct(lat1, lon1, azi1, s12);
        assert!((dlat - lat2).abs() < 1e-12, "lat {dlat}");
        assert!((dlon - lon2).abs() < 1e-12, "lon {dlon}");
    }

    #[test]
    fn nearly_antipodal_converges() {
        // The classic Vincenty-failure regime; Karney must converge and the
        // distance must be slightly under half the ellipsoid circumference.
        let k = wgs84();
        let (s12, azi1, azi2) = k.inverse(0.5, 0.0, -0.4999, 179.9999);
        assert!(s12.is_finite() && azi1.is_finite() && azi2.is_finite());
        let half_meridian = PI * Ellipsoid::WGS84.b(); // lower bound scale
        assert!(s12 > half_meridian, "{s12}");
        assert!(s12 < PI * Ellipsoid::WGS84.a, "{s12}");
    }
}
