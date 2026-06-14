//! [`LineKernel`] — the precomputed origin-and-azimuth state of the direct
//! problem, ready to place a point at any distance ([`LineKernel::position`])
//! without re-deriving the series. Built by [`Karney::line`](super::solver::Karney::line).

use super::math::{ang_normalize, atan2d, sin_cos_series, sq, tiny};
use super::series::{CARR_SIZE, GEODESIC_ORDER};

/// A geodesic line: the precomputed origin-and-azimuth state of the direct
/// problem, ready to place points at any distance ([`position`](Self::position))
/// without re-deriving the series. The kernel behind
/// [`crate::earth::GeodesicLine`].
pub(crate) struct LineKernel {
    pub(super) f: f64,
    pub(super) f1: f64,
    pub(super) b: f64,
    pub(super) lon1: f64,
    pub(super) k2: f64,
    pub(super) salp0: f64,
    pub(super) calp0: f64,
    pub(super) ssig1: f64,
    pub(super) csig1: f64,
    pub(super) somg1: f64,
    pub(super) comg1: f64,
    pub(super) stau1: f64,
    pub(super) ctau1: f64,
    pub(super) a1m1: f64,
    pub(super) b11: f64,
    pub(super) a3c: f64,
    pub(super) b31: f64,
    pub(super) c1a: [f64; CARR_SIZE],
    pub(super) c1pa: [f64; CARR_SIZE],
    pub(super) c3a: [f64; GEODESIC_ORDER],
}

impl LineKernel {
    /// The point `s12` metres along the line (negative is backwards), as
    /// `(lat2, lon2, azi2)` with `azi2` the forward azimuth there.
    pub(crate) fn position(&self, s12: f64) -> (f64, f64, f64) {
        let (ssig1, csig1) = (self.ssig1, self.csig1);
        let (somg1, comg1) = (self.somg1, self.comg1);
        let (salp0, calp0) = (self.salp0, self.calp0);
        let (a1m1, b11, a3c, b31, k2) = (self.a1m1, self.b11, self.a3c, self.b31, self.k2);
        let (stau1, ctau1) = (self.stau1, self.ctau1);
        let (c1a, c1pa, c3a) = (self.c1a, self.c1pa, self.c3a);
        let lon1 = self.lon1;

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
