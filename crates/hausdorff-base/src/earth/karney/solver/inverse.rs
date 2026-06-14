//! The inverse-problem drivers: [`inverse`](Karney::inverse) and the shared
//! [`gen_inverse`](Karney::gen_inverse) with its Newton/bisection iteration and
//! geodesic-area accumulation, plus the polygon-area accumulator
//! ([`ring_perimeter_area`](Karney::ring_perimeter_area)) built on top.

use std::f64::consts::{FRAC_1_SQRT_2, PI};

use super::Karney;
use crate::earth::karney::math::{
    MAXIT1, MAXIT2, ang_diff, ang_round, atan2d, lat_fix, norm2, sin_cos_series, sincosd, sq, tiny,
    tol0, tolb, transit,
};
use crate::earth::karney::series::{CARR_SIZE, GEODESIC_ORDER};

impl Karney {
    /// The inverse geodesic problem: the shortest-path distance in metres and
    /// the forward azimuths (degrees, clockwise from north) at both endpoints,
    /// for two points in latitude/longitude degrees. Returns
    /// `(s12, azi1, azi2)`; NaN for a non-finite or out-of-range latitude.
    pub(crate) fn inverse(&self, lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> (f64, f64, f64) {
        let (s12, salp1, calp1, salp2, calp2, _) = self.gen_inverse(lat1, lon1, lat2, lon2, false);
        (s12, atan2d(salp1, calp1), atan2d(salp2, calp2))
    }

    /// The inverse problem with the geodesic area: `(s12, S12)` where `S12` is
    /// the signed area (m²) between the geodesic and the equator — the
    /// edge-wise quantity [`ring_perimeter_area`](Self::ring_perimeter_area)
    /// accumulates.
    fn inverse_with_area(&self, lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> (f64, f64) {
        let (s12, _, _, _, _, s12_area) = self.gen_inverse(lat1, lon1, lat2, lon2, true);
        (s12, s12_area)
    }

    /// The shared inverse driver. Returns
    /// `(s12, salp1, calp1, salp2, calp2, S12)`; `S12` (the area between the
    /// geodesic and the equator, clockwise-positive before the final sign
    /// fix-up) is NaN unless `want_area`.
    fn gen_inverse(
        &self,
        lat1: f64,
        lon1: f64,
        lat2: f64,
        lon2: f64,
        want_area: bool,
    ) -> (f64, f64, f64, f64, f64, f64) {
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
        // ω₁₂ bookkeeping for the area path; somg12 > 1 marks "not yet set"
        // (fall back to omg12 at the end, as in the reference).
        let mut somg12 = 2.0;
        let mut comg12 = 0.0;
        let mut omg12 = 0.0;

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
            omg12 = lam12 / self.f1;
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
                omg12 = lam12 / (self.f1 * dnm);
            } else {
                // Newton's method on alp1, bracketed by bisection.
                let mut tripn = false;
                let mut tripb = false;
                let mut salp1a = tiny();
                let mut calp1a = 1.0;
                let mut salp1b = tiny();
                let mut calp1b = -1.0;
                let mut eps = 0.0;
                let mut domg12 = 0.0;
                for numit in 0..MAXIT2 {
                    let (
                        v,
                        salp2i,
                        calp2i,
                        sig12i,
                        ssig1i,
                        csig1i,
                        ssig2i,
                        csig2i,
                        epsi,
                        domg12i,
                        dv,
                    ) = self.lambda12(
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
                    domg12 = domg12i;

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
                if want_area {
                    let (sdomg12, cdomg12) = domg12.sin_cos();
                    somg12 = slam12 * cdomg12 - clam12 * sdomg12;
                    comg12 = clam12 * cdomg12 + slam12 * sdomg12;
                }
            }
        }

        let mut s12_area = f64::NAN;
        if want_area {
            // The area S12 between the geodesic and the equator: the C4 series
            // integral plus the spherical-excess term c2·α₁₂, computed in the
            // canonical frame and sign-fixed with the swaps below.
            let salp0 = salp1 * cbet1;
            let calp0 = calp1.hypot(salp1 * sbet1);
            s12_area = if calp0 != 0.0 && salp0 != 0.0 {
                let mut ssig1 = sbet1;
                let mut csig1 = calp1 * cbet1;
                let mut ssig2 = sbet2;
                let mut csig2 = calp2 * cbet2;
                let k2 = sq(calp0) * self.ep2;
                let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
                let a4 = sq(self.a) * calp0 * salp0 * self.e2;
                norm2(&mut ssig1, &mut csig1);
                norm2(&mut ssig2, &mut csig2);
                let mut c4a = [0.0; GEODESIC_ORDER];
                self.c4f(eps, &mut c4a);
                let b41 = sin_cos_series(false, ssig1, csig1, &c4a);
                let b42 = sin_cos_series(false, ssig2, csig2, &c4a);
                a4 * (b42 - b41)
            } else {
                0.0
            };

            if !meridian && somg12 > 1.0 {
                somg12 = omg12.sin();
                comg12 = omg12.cos();
            }

            let alp12 = if !meridian && comg12 > -FRAC_1_SQRT_2 && sbet2 - sbet1 < 1.75 {
                // Geodesic runs through the vertex: the spherical-excess
                // half-angle form is the numerically stable one here.
                let domg12 = 1.0 + comg12;
                let dbet1 = 1.0 + cbet1;
                let dbet2 = 1.0 + cbet2;
                2.0 * (somg12 * (sbet1 * dbet2 + sbet2 * dbet1))
                    .atan2(domg12 * (sbet1 * sbet2 + dbet1 * dbet2))
            } else {
                let mut salp12 = salp2 * calp1 - calp2 * salp1;
                let mut calp12 = calp2 * calp1 + salp2 * salp1;
                // alp12 = 0 with calp12 < 0 means alp12 = ±π: pick the sign
                // consistent with alp1 so a half-planet edge doesn't flip.
                if salp12 == 0.0 && calp12 < 0.0 {
                    salp12 = tiny() * calp1;
                    calp12 = -1.0;
                }
                salp12.atan2(calp12)
            };
            s12_area += self.c2 * alp12;
            s12_area *= swapp * lonsign * latsign;
            s12_area += 0.0;
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

        (s12x, salp1, calp1, salp2, calp2, s12_area)
    }

    /// The ellipsoid's total surface area, m².
    pub(crate) fn total_area(&self) -> f64 {
        4.0 * PI * self.c2
    }

    /// Perimeter (metres) and signed area (m²) of the polygon whose vertices
    /// are `ring` as `(lat, lon)` degrees — GeographicLib's `PolygonArea`
    /// accumulator. The ring closes itself (a trailing duplicate of the first
    /// vertex contributes a zero-length edge, so closed rings are fine too).
    ///
    /// `ccw_positive` picks the winding that counts as positive. With
    /// `signed`, a ring wound the other way comes back negative (wrong for
    /// rings enclosing more than half the planet); without it, the area is
    /// always the enclosed area in `[0, total)`.
    pub(crate) fn ring_perimeter_area(
        &self,
        ring: &[(f64, f64)],
        ccw_positive: bool,
        signed: bool,
    ) -> (f64, f64) {
        let Some((&(first_lat, first_lon), rest)) = ring.split_first() else {
            return (0.0, 0.0);
        };
        let mut perimeter = 0.0;
        let mut area = 0.0;
        let mut crossings: i64 = 0;
        let (mut last_lat, mut last_lon) = (first_lat, first_lon);
        for &(lat, lon) in rest.iter().chain(std::iter::once(&(first_lat, first_lon))) {
            let (s12, s12_area) = self.inverse_with_area(last_lat, last_lon, lat, lon);
            perimeter += s12;
            area += s12_area;
            crossings += transit(last_lon, lon);
            (last_lat, last_lon) = (lat, lon);
        }

        // Fold the accumulated S12 into a polygon area: account for prime-
        // meridian crossings (pole encirclement), fix the winding convention,
        // and reduce modulo the ellipsoid's total area.
        let area0 = 4.0 * PI * self.c2;
        let mut area = area % area0;
        if crossings % 2 != 0 {
            area += if area < 0.0 {
                area0 / 2.0
            } else {
                -area0 / 2.0
            };
        }
        // S12 accumulates in the clockwise sense.
        if ccw_positive {
            area = -area;
        }
        if signed {
            if area > area0 / 2.0 {
                area -= area0;
            } else if area <= -area0 / 2.0 {
                area += area0;
            }
        } else if area < 0.0 {
            area += area0;
        }
        (perimeter, area)
    }
}
