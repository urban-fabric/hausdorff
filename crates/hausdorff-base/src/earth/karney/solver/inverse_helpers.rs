//! Private helpers for the inverse driver: the geodesic length integrals
//! ([`lengths`](Karney::lengths)), the Newton-iteration starting guess
//! ([`inverse_start`](Karney::inverse_start)), and the per-trial longitude
//! mismatch ([`lambda12`](Karney::lambda12)) that the Newton residual is built
//! from.

use std::f64::consts::PI;

use super::Karney;
use crate::earth::karney::math::{astroid, norm2, sin_cos_series, sq, tiny, tol1, xthresh};
use crate::earth::karney::series::{CARR_SIZE, GEODESIC_ORDER, a1m1f, a2m1f, c1f, c2f};

impl Karney {
    /// The length integrals between two points on a geodesic, in units of `b`:
    /// the distance `s12b` (if `want_distance`) and the reduced length `m12b`
    /// with the series constant `m0` (if `want_reduced`); unsolicited outputs
    /// are NaN. The bool pair replaces GeographicLib's capability mask.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn lengths(
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
    pub(super) fn inverse_start(
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
    /// `(v, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dv)`,
    /// where `dv` is `dv/dalp1` when `diffp` (else NaN) and `domg12` is the
    /// ω−λ correction the area path needs.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn lambda12(
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
    ) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
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
            lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12,
        )
    }
}
