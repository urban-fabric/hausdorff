//! The geodesic solver: the [`Karney`] per-ellipsoid state (derived constants
//! plus the `A₃`/`C₃`/`C₄` series coefficients) and the inverse/direct drivers —
//! `inverse`, `gen_inverse`, the Newton/bisection iteration, the polygon-area
//! accumulator, and `line` (which hands off to [`LineKernel`](super::line::LineKernel)).
//!
//! The driver logic is split by problem: the per-ellipsoid state and the
//! series accessors live here, the inverse helpers (`lengths`, `inverse_start`,
//! `lambda12`) in [`inverse_helpers`], the inverse drivers and polygon-area
//! accumulator in [`inverse`], and the direct problem (`direct`/`line`) in
//! [`direct`].

mod direct;
mod inverse;
mod inverse_helpers;

use std::sync::LazyLock;

use crate::earth::ellipsoid::Ellipsoid;

use super::math::{eatanhe, polyval, sq, tol2};
use super::series::{COEFF_A3, COEFF_C3, COEFF_C4, GEODESIC_ORDER, NC3X, NC4X};

/// The geodesic solver for one ellipsoid: the ellipsoid's derived constants
/// plus the `A₃`/`C₃` series coefficients, which depend only on the third
/// flattening and so are computed once per ellipsoid.
pub(crate) struct Karney {
    pub(super) a: f64,
    pub(super) f: f64,
    pub(super) f1: f64,
    pub(super) e2: f64,
    pub(super) ep2: f64,
    pub(super) n: f64,
    pub(super) b: f64,
    /// Authalic radius squared — half the ellipsoid's total area is `2π·c2`.
    pub(super) c2: f64,
    pub(super) etol2: f64,
    pub(super) a3x: [f64; GEODESIC_ORDER],
    pub(super) c3x: [f64; NC3X],
    pub(super) c4x: [f64; NC4X],
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
        let e2 = ellipsoid.e2();
        let ep2 = ellipsoid.ep2();
        let n = ellipsoid.n();
        let b = ellipsoid.b();
        let c2 = (sq(a)
            + sq(b)
                * (if e2 == 0.0 {
                    1.0
                } else {
                    eatanhe(1.0, (if f < 0.0 { -1.0 } else { 1.0 }) * e2.abs().sqrt()) / e2
                }))
            / 2.0;
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

        let mut c4x = [0.0; NC4X];
        let mut o = 0;
        let mut k = 0;
        for l in 0..GEODESIC_ORDER {
            for j in (l..GEODESIC_ORDER).rev() {
                let m = GEODESIC_ORDER - j - 1;
                c4x[k] = polyval(m, &COEFF_C4[o..], n) / COEFF_C4[o + m + 1];
                k += 1;
                o += m + 2;
            }
        }

        Karney {
            a,
            f,
            f1,
            e2,
            ep2,
            n,
            b,
            c2,
            etol2,
            a3x,
            c3x,
            c4x,
        }
    }

    /// `A₃(ε)`: the longitude-series scale.
    pub(super) fn a3f(&self, eps: f64) -> f64 {
        polyval(GEODESIC_ORDER - 1, &self.a3x, eps)
    }

    /// The `C₃` Fourier coefficients (longitude series) at `eps`.
    pub(super) fn c3f(&self, eps: f64, c: &mut [f64; GEODESIC_ORDER]) {
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

    /// The `C₄` Fourier coefficients (area series) at `eps`. Unlike the other
    /// series this one is 0-indexed (used with `sinp = false`).
    pub(super) fn c4f(&self, eps: f64, c: &mut [f64; GEODESIC_ORDER]) {
        let mut mult = 1.0;
        let mut o = 0;
        #[allow(clippy::needless_range_loop)]
        for l in 0..GEODESIC_ORDER {
            let m = GEODESIC_ORDER - l - 1;
            c[l] = mult * polyval(m, &self.c4x[o..], eps);
            o += m + 1;
            mult *= eps;
        }
    }
}
