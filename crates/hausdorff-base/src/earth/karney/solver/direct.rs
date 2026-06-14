//! The direct-problem entry points: [`direct`](Karney::direct) and
//! [`line`](Karney::line), the origin-dependent set-up that builds the
//! [`LineKernel`](crate::earth::karney::line::LineKernel) a caller placing many
//! points along one geodesic pays for once.

use super::Karney;
use crate::earth::karney::line::LineKernel;
use crate::earth::karney::math::{
    ang_normalize, ang_round, lat_fix, norm2, sin_cos_series, sincosd, sq, tiny,
};
use crate::earth::karney::series::{CARR_SIZE, GEODESIC_ORDER, a1m1f, c1f, c1pf};

impl Karney {
    /// The direct geodesic problem: from a point in latitude/longitude
    /// degrees, travel `s12` metres along the geodesic with initial azimuth
    /// `azi1` (degrees, clockwise from north). Returns `(lat2, lon2, azi2)`
    /// with `azi2` the forward azimuth at the destination; NaN for a
    /// non-finite or out-of-range latitude.
    pub(crate) fn direct(&self, lat1: f64, lon1: f64, azi1: f64, s12: f64) -> (f64, f64, f64) {
        self.line(lat1, lon1, azi1).position(s12)
    }

    /// Set up the geodesic line from a point with the given initial azimuth —
    /// the origin-dependent half of the direct problem, factored out so a
    /// caller placing many points along one geodesic
    /// ([`LineKernel::position`]) pays it once.
    pub(crate) fn line(&self, lat1: f64, lon1: f64, azi1: f64) -> LineKernel {
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

        LineKernel {
            f: self.f,
            f1: self.f1,
            b: self.b,
            lon1,
            k2,
            salp0,
            calp0,
            ssig1,
            csig1,
            somg1,
            comg1,
            stau1,
            ctau1,
            a1m1,
            b11,
            a3c,
            b31,
            c1a,
            c1pa,
            c3a,
        }
    }
}
