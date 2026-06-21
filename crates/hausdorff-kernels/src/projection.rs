//! Map-projection kernels over any [`LazyReal`] scalar — geographic
//! (lat/lon degrees) ↔ projected planar metres, written once so the same
//! source runs on `f64`, `f32`, and tensor/expression backends.
//!
//! These are the backend-neutral peers of the concrete `Projection` impls in
//! `hausdorff-earth`: that crate keeps the `Ellipsoid`/`GeoPoint`/`Coord`-typed
//! API and delegates the arithmetic here, exactly as the earth geodesic
//! wrappers delegate to [`crate::earth`]. The per-point coordinate is the lazy
//! scalar `T`; an ellipsoid's constants (the eccentricity, the Krüger series,
//! the false origin) are `f64` — they depend only on the datum, not on the
//! point being transformed — and are lifted per use with
//! [`scalar_like`](LazyReal::scalar_like), the same pattern [`crate::earth`]
//! uses for [`EARTH_RADIUS_M`](crate::earth::EARTH_RADIUS_M).
//!
//! Provided: [`WebMercator`](web_mercator_forward) (spherical Mercator) and the
//! Karney–Krüger [`TransverseMercator`](transverse_mercator_forward) (the basis
//! of UTM). Following Karney 2011, *Transverse Mercator with an accuracy of a
//! few nanometers*, and Snyder, *Map Projections — A Working Manual* (USGS
//! PP1395, 1987).

use core::f64::consts::{FRAC_PI_2, FRAC_PI_4};

use euclid::numeric::LazyReal;

use crate::earth::LatLon;

// ---- Web Mercator (spherical) ----

/// Web Mercator (EPSG:3857) forward: geographic lat/lon **degrees** to projected
/// `(x, y)` metres on a sphere of `radius` metres.
///
/// `y = R·ψ(φ)` with `ψ` the isometric latitude `log((1 + sinφ)/cosφ)` — the
/// exact spherical-Mercator form (`riemann` `Hausdorff/Rhumb.isoLat`).
#[must_use]
pub fn web_mercator_forward<T: LazyReal>(lat_deg: T, lon_deg: T, radius: T) -> (T, T) {
    let one = lat_deg.scalar_like(1.0);
    let (s, c) = lat_deg.to_radians().sin_cos();
    let x = radius.clone() * lon_deg.to_radians();
    let y = radius * ((one + s) / c).ln();
    (x, y)
}

/// Web Mercator inverse: projected `(x, y)` metres back to geographic lat/lon
/// **degrees**, on a sphere of `radius` metres.
///
/// `φ = arcsin(tanh(y/R))` — the provably exact inverse (`riemann`
/// `Hausdorff/Mercator.mercator_inverse`), since `tanh(ψ(φ)) = sinφ`.
#[must_use]
pub fn web_mercator_inverse<T: LazyReal>(x: T, y: T, radius: T) -> (T, T) {
    let lat_deg = (y / radius.clone()).tanh().asin().to_degrees();
    let lon_deg = (x / radius).to_degrees();
    (lat_deg, lon_deg)
}

/// Project a [`LatLon`] provider through [`web_mercator_forward`] — the entry
/// point a tensor/columnar lat/lon source uses to project en masse.
#[must_use]
pub fn web_mercator_project_at<T, P>(p: &P, radius: T) -> (T, T)
where
    T: LazyReal,
    P: LatLon<T>,
{
    web_mercator_forward(p.lat_deg(), p.lon_deg(), radius)
}

// ---- Transverse Mercator (Karney–Krüger, ellipsoidal) ----

/// Karney–Krüger Transverse Mercator series constants for one ellipsoid: the
/// rectifying radius `A`, the forward (`α`) and inverse (`β`) coefficient sets
/// (order `n⁶`), and the eccentricity `e`. They depend only on the ellipsoid,
/// so they are computed once in `f64` and lifted per point with `scalar_like`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KrugerSeries {
    /// Rectifying radius `A = a/(1+n)·(1 + n²/4 + n⁴/64 + n⁶/256)`.
    pub a_rect: f64,
    /// Forward coefficients α₁..α₆.
    pub alpha: [f64; 6],
    /// Inverse coefficients β₁..β₆.
    pub beta: [f64; 6],
    /// First eccentricity `e = √(f(2−f))`.
    pub e: f64,
}

impl KrugerSeries {
    /// The series for an ellipsoid given its semi-major axis `a` and flattening
    /// `f` (order `n⁶`, Karney 2011).
    #[must_use]
    pub fn for_ellipsoid(a: f64, f: f64) -> Self {
        let n = f / (2.0 - f); // third flattening
        let (n2, n3, n4, n5, n6) = (n * n, n * n * n, n.powi(4), n.powi(5), n.powi(6));
        let a_rect = a / (1.0 + n) * (1.0 + n2 / 4.0 + n4 / 64.0 + n6 / 256.0);
        let alpha = [
            n / 2.0 - 2.0 / 3.0 * n2 + 5.0 / 16.0 * n3 + 41.0 / 180.0 * n4 - 127.0 / 288.0 * n5
                + 7891.0 / 37800.0 * n6,
            13.0 / 48.0 * n2 - 3.0 / 5.0 * n3 + 557.0 / 1440.0 * n4 + 281.0 / 630.0 * n5
                - 1_983_433.0 / 1_935_360.0 * n6,
            61.0 / 240.0 * n3 - 103.0 / 140.0 * n4
                + 15_061.0 / 26_880.0 * n5
                + 167_603.0 / 181_440.0 * n6,
            49_561.0 / 161_280.0 * n4 - 179.0 / 168.0 * n5 + 6_601_661.0 / 7_257_600.0 * n6,
            34_729.0 / 80_640.0 * n5 - 3_418_889.0 / 1_995_840.0 * n6,
            212_378_941.0 / 319_334_400.0 * n6,
        ];
        let beta = [
            n / 2.0 - 2.0 / 3.0 * n2 + 37.0 / 96.0 * n3 - 1.0 / 360.0 * n4 - 81.0 / 512.0 * n5
                + 96_199.0 / 604_800.0 * n6,
            1.0 / 48.0 * n2 + 1.0 / 15.0 * n3 - 437.0 / 1440.0 * n4 + 46.0 / 105.0 * n5
                - 1_118_711.0 / 3_870_720.0 * n6,
            17.0 / 480.0 * n3 - 37.0 / 840.0 * n4 - 209.0 / 4480.0 * n5 + 5569.0 / 90_720.0 * n6,
            4397.0 / 161_280.0 * n4 - 11.0 / 504.0 * n5 - 830_251.0 / 7_257_600.0 * n6,
            4583.0 / 161_280.0 * n5 - 108_847.0 / 3_991_680.0 * n6,
            20_648_693.0 / 638_668_800.0 * n6,
        ];
        Self {
            a_rect,
            alpha,
            beta,
            e: (f * (2.0 - f)).sqrt(),
        }
    }

    /// ξ at the natural origin (latitude `lat_origin_deg`, on the central
    /// meridian) — the rectifying-latitude offset subtracted so northings are
    /// measured from the origin parallel. An ellipsoid+origin constant in `f64`.
    fn xi_origin(&self, lat_origin_deg: f64) -> f64 {
        if lat_origin_deg == 0.0 {
            return 0.0;
        }
        let xip = taupf(lat_origin_deg.to_radians().tan(), self.e).atan();
        let mut xi = xip;
        for (j, &a_j) in self.alpha.iter().enumerate() {
            xi += a_j * (2.0 * (j as f64 + 1.0) * xip).sin();
        }
        xi
    }
}

/// Transverse Mercator forward: geographic lat/lon **degrees** to projected
/// `(easting, northing)` metres, via the Karney–Krüger series. `series` carries
/// the ellipsoid constants; the remaining `f64` parameters are the central
/// meridian, latitude of origin, central-meridian scale `k0`, and false origin.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn transverse_mercator_forward<T: LazyReal>(
    lat_deg: T,
    lon_deg: T,
    central_meridian: f64,
    lat_origin: f64,
    k0: f64,
    false_easting: f64,
    false_northing: f64,
    series: &KrugerSeries,
) -> (T, T) {
    let e = lat_deg.scalar_like(series.e);
    // Conformal latitude (as its tangent τ′), then the Gauss–Schreiber
    // intermediate (ξ′, η′) about the central meridian.
    let lam = (lon_deg - lat_deg.scalar_like(central_meridian)).to_radians();
    let (slam, clam) = lam.sin_cos();
    let taup = taupf(lat_deg.clone().to_radians().tan(), e);
    let xip = taup.clone().atan2(clam.clone());
    let etap = (slam / taup.hypot(clam)).asinh();
    // Krüger series ξ′,η′ → ξ,η.
    let (mut xi, mut eta) = (xip.clone(), etap.clone());
    for (j, &a_j) in series.alpha.iter().enumerate() {
        let two_j = 2.0 * (j as f64 + 1.0);
        let a_j = lat_deg.scalar_like(a_j);
        let two_xip = lat_deg.scalar_like(two_j) * xip.clone();
        let two_etap = lat_deg.scalar_like(two_j) * etap.clone();
        xi = xi + a_j.clone() * two_xip.clone().sin() * two_etap.clone().cosh();
        eta = eta + a_j * two_xip.cos() * two_etap.sinh();
    }
    let scale = lat_deg.scalar_like(k0 * series.a_rect);
    let xi_origin = lat_deg.scalar_like(series.xi_origin(lat_origin));
    let easting = lat_deg.scalar_like(false_easting) + scale.clone() * eta;
    let northing = lat_deg.scalar_like(false_northing) + scale * (xi - xi_origin);
    (easting, northing)
}

/// Transverse Mercator inverse: projected `(easting, northing)` metres back to
/// geographic lat/lon **degrees**. Inverse-pair to
/// [`transverse_mercator_forward`].
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn transverse_mercator_inverse<T: LazyReal>(
    easting: T,
    northing: T,
    central_meridian: f64,
    lat_origin: f64,
    k0: f64,
    false_easting: f64,
    false_northing: f64,
    series: &KrugerSeries,
) -> (T, T) {
    let e = easting.scalar_like(series.e);
    let scale = easting.scalar_like(k0 * series.a_rect);
    let xi = easting.scalar_like(series.xi_origin(lat_origin))
        + (northing - easting.scalar_like(false_northing)) / scale.clone();
    let eta = (easting.clone() - easting.scalar_like(false_easting)) / scale;
    // Krüger series ξ,η → ξ′,η′.
    let (mut xip, mut etap) = (xi.clone(), eta.clone());
    for (j, &b_j) in series.beta.iter().enumerate() {
        let two_j = 2.0 * (j as f64 + 1.0);
        let b_j = easting.scalar_like(b_j);
        let two_xi = easting.scalar_like(two_j) * xi.clone();
        let two_eta = easting.scalar_like(two_j) * eta.clone();
        xip = xip - b_j.clone() * two_xi.clone().sin() * two_eta.clone().cosh();
        etap = etap - b_j * two_xi.cos() * two_eta.sinh();
    }
    // Invert the Gauss–Schreiber step: recover λ and the conformal latitude
    // (as τ′), then the geodetic latitude via `tauf`.
    let (sxip, cxip) = xip.sin_cos();
    let sinh_etap = etap.sinh();
    let lam = sinh_etap.clone().atan2(cxip.clone());
    let taup = sxip / sinh_etap.hypot(cxip);
    let lat_deg = tauf(taup, e).atan().to_degrees();
    let lon_deg = easting.scalar_like(central_meridian) + lam.to_degrees();
    (lat_deg, lon_deg)
}

/// Project a [`LatLon`] provider through [`transverse_mercator_forward`] — the
/// entry point a tensor/columnar lat/lon source uses to project en masse.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn transverse_mercator_project_at<T, P>(
    p: &P,
    central_meridian: f64,
    lat_origin: f64,
    k0: f64,
    false_easting: f64,
    false_northing: f64,
    series: &KrugerSeries,
) -> (T, T)
where
    T: LazyReal,
    P: LatLon<T>,
{
    transverse_mercator_forward(
        p.lat_deg(),
        p.lon_deg(),
        central_meridian,
        lat_origin,
        k0,
        false_easting,
        false_northing,
        series,
    )
}

// ---- Shared auxiliary-latitude helpers (Karney 2011) ----

/// `tan χ` from `tan φ`: the tangent of the **conformal latitude** χ given the
/// tangent of the geodetic latitude φ, on an ellipsoid of eccentricity `e`
/// (Karney 2011 `taupf`). The exact closed form (no series).
#[must_use]
pub fn taupf<T: LazyReal>(tau: T, e: T) -> T {
    let one = tau.scalar_like(1.0);
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
pub fn tauf<T: LazyReal>(taup: T, e: T) -> T {
    let one = taup.scalar_like(1.0);
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
pub fn m_factor<T: LazyReal>(phi: T, e: T) -> T {
    let one = phi.scalar_like(1.0);
    let e2 = e.clone() * e;
    let s = phi.clone().sin();
    phi.cos() / (one - e2 * s.clone() * s).sqrt()
}

/// The conformal isometric building block `t = tan(π/4 − φ/2) / [(1 − e sinφ)/(1 +
/// e sinφ)]^{e/2}` (Snyder 15-9) — `t = tan(π/4 − χ/2)`, so `χ = π/2 − 2·arctan t`.
#[must_use]
pub fn conformal_t<T: LazyReal>(phi: T, e: T) -> T {
    let one = phi.scalar_like(1.0);
    let half = phi.scalar_like(0.5);
    let frac_pi_4 = phi.scalar_like(FRAC_PI_4);
    let es = phi.clone().sin() * e.clone();
    let exponent = e * half.clone();
    (frac_pi_4 - phi * half).tan() / ((one.clone() - es.clone()) / (one + es)).powf(exponent)
}

/// Invert the conformal `t` back to geodetic latitude `φ` (radians) by Snyder's
/// fixed-point iteration (7-9). Its `e → 0` limit is the exact closed form
/// `φ = π/2 − 2·arctan t`.
#[must_use]
pub fn lat_from_conformal_t<T: LazyReal>(t: T, e: T) -> T {
    let one = t.scalar_like(1.0);
    let two = t.scalar_like(2.0);
    let frac_pi_2 = t.scalar_like(FRAC_PI_2);
    let mut phi = frac_pi_2.clone() - two.clone() * t.clone().atan();
    for _ in 0..15 {
        let es = phi.clone().sin() * e.clone();
        let exponent = e.clone() * phi.scalar_like(0.5);
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
pub fn authalic_q<T: LazyReal>(phi: T, e: T) -> T {
    let one = phi.scalar_like(1.0);
    let two = phi.scalar_like(2.0);
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
pub fn lat_from_authalic_q<T: LazyReal>(q: T, e: T) -> T {
    let one = q.scalar_like(1.0);
    let two = q.scalar_like(2.0);
    let half = q.scalar_like(0.5);
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

// ---- Conic / azimuthal / equal-area projection families ----
//
// Each family's per-projection constants (cone constant, scale, origin radius)
// depend only on the ellipsoid and standard parallels, so they are computed once
// in `f64` (a `*Constants` struct, the peer of `KrugerSeries`) and lifted per
// point with `scalar_like`; only the coordinate is the lazy `T`. The constants
// constructors reuse the auxiliary-latitude kernels above at `T = f64`.

/// Lambert Conformal Conic (2SP) constants: eccentricity, cone constant `n`,
/// scale `F`, and origin radius `ρ0`, plus the semi-major axis and false origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LccConstants {
    pub a: f64,
    pub e: f64,
    pub n: f64,
    pub f: f64,
    pub rho0: f64,
    pub lon0: f64,
    pub false_easting: f64,
    pub false_northing: f64,
}

impl LccConstants {
    /// Solve the cone constant `n`, scale `F`, and origin radius `ρ0` from the two
    /// standard parallels (Snyder pp. 107–109).
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        a: f64,
        e: f64,
        lat1: f64,
        lat2: f64,
        lat0: f64,
        lon0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        let (phi1, phi2) = (lat1.to_radians(), lat2.to_radians());
        let (m1, m2) = (m_factor(phi1, e), m_factor(phi2, e));
        let (t0, t1, t2) = (
            conformal_t(lat0.to_radians(), e),
            conformal_t(phi1, e),
            conformal_t(phi2, e),
        );
        let n = if (phi1 - phi2).abs() < 1e-12 {
            phi1.sin()
        } else {
            (m1.ln() - m2.ln()) / (t1.ln() - t2.ln())
        };
        let f = m1 / (n * t1.powf(n));
        let rho0 = a * f * t0.powf(n);
        Self {
            a,
            e,
            n,
            f,
            rho0,
            lon0,
            false_easting,
            false_northing,
        }
    }
}

/// Lambert Conformal Conic forward: geographic lat/lon **degrees** → `(x, y)` m.
#[must_use]
pub fn lcc_forward<T: LazyReal>(lat_deg: T, lon_deg: T, c: &LccConstants) -> (T, T) {
    let e = lat_deg.scalar_like(c.e);
    let t = conformal_t(lat_deg.clone().to_radians(), e);
    let rho = lat_deg.scalar_like(c.a * c.f) * t.powf(lat_deg.scalar_like(c.n));
    let theta = (lon_deg - lat_deg.scalar_like(c.lon0)).to_radians() * lat_deg.scalar_like(c.n);
    let (st, ct) = theta.sin_cos();
    let x = lat_deg.scalar_like(c.false_easting) + rho.clone() * st;
    let y = lat_deg.scalar_like(c.false_northing + c.rho0) - rho * ct;
    (x, y)
}

/// Lambert Conformal Conic inverse: `(x, y)` m → geographic lat/lon **degrees**.
#[must_use]
pub fn lcc_inverse<T: LazyReal>(x_in: T, y_in: T, c: &LccConstants) -> (T, T) {
    let e = x_in.scalar_like(c.e);
    let x = x_in.clone() - x_in.scalar_like(c.false_easting);
    let y0 = x_in.scalar_like(c.rho0) - (y_in - x_in.scalar_like(c.false_northing));
    let nsign = c.n.signum();
    let rho = x.clone().hypot(y0.clone()).copysign(x_in.scalar_like(c.n));
    let theta = (x.clone() * x_in.scalar_like(nsign)).atan2(y0 * x_in.scalar_like(nsign));
    let t = (rho / x_in.scalar_like(c.a * c.f)).powf(x_in.scalar_like(1.0 / c.n));
    let lat = lat_from_conformal_t(t, e).to_degrees();
    let lon = (theta / x_in.scalar_like(c.n)).to_degrees() + x_in.scalar_like(c.lon0);
    (lat, lon)
}

/// Polar Stereographic (ellipsoidal) constants: eccentricity, aspect sign `s`
/// (north `+1` / south `−1`), and the standard-parallel `tc`/`mc`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarStereoConstants {
    pub a: f64,
    pub e: f64,
    pub s: f64,
    pub tc: f64,
    pub mc: f64,
    pub lon0: f64,
    pub false_easting: f64,
    pub false_northing: f64,
}

impl PolarStereoConstants {
    /// Constants pinned to the latitude of true scale `lat_ts` (Snyder pp. 160–161).
    #[must_use]
    pub fn new(
        a: f64,
        e: f64,
        lat_ts: f64,
        lon0: f64,
        north: bool,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        let s = if north { 1.0 } else { -1.0 };
        let phic = s * lat_ts.to_radians();
        Self {
            a,
            e,
            s,
            tc: conformal_t(phic, e),
            mc: m_factor(phic, e),
            lon0,
            false_easting,
            false_northing,
        }
    }
}

/// Polar Stereographic forward: geographic lat/lon **degrees** → `(x, y)` m.
#[must_use]
pub fn polar_stereographic_forward<T: LazyReal>(
    lat_deg: T,
    lon_deg: T,
    c: &PolarStereoConstants,
) -> (T, T) {
    let e = lat_deg.scalar_like(c.e);
    let phi = lat_deg.clone().to_radians() * lat_deg.scalar_like(c.s);
    let rho = lat_deg.scalar_like(c.a * c.mc / c.tc) * conformal_t(phi, e);
    let dlon = (lon_deg - lat_deg.scalar_like(c.lon0)).to_radians();
    let (sdl, cdl) = dlon.sin_cos();
    let x = lat_deg.scalar_like(c.false_easting) + rho.clone() * sdl;
    let y = lat_deg.scalar_like(c.false_northing) - rho * cdl * lat_deg.scalar_like(c.s);
    (x, y)
}

/// Polar Stereographic inverse: `(x, y)` m → geographic lat/lon **degrees**.
#[must_use]
pub fn polar_stereographic_inverse<T: LazyReal>(
    x_in: T,
    y_in: T,
    c: &PolarStereoConstants,
) -> (T, T) {
    let e = x_in.scalar_like(c.e);
    let x = x_in.clone() - x_in.scalar_like(c.false_easting);
    let y = y_in - x_in.scalar_like(c.false_northing);
    let rho = x.clone().hypot(y.clone());
    let t = rho * x_in.scalar_like(c.tc) / x_in.scalar_like(c.a * c.mc);
    let lat = (lat_from_conformal_t(t, e) * x_in.scalar_like(c.s)).to_degrees();
    let lon = x_in.scalar_like(c.lon0) + x.atan2(y * x_in.scalar_like(-c.s)).to_degrees();
    (lat, lon)
}

/// Albers Equal-Area Conic constants: eccentricity, cone constant `n`, the
/// constant `C`, and origin radius `ρ0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlbersConstants {
    pub a: f64,
    pub e: f64,
    pub n: f64,
    pub cc: f64,
    pub rho0: f64,
    pub lon0: f64,
    pub false_easting: f64,
    pub false_northing: f64,
}

impl AlbersConstants {
    /// Solve the cone constant `n`, `C`, and origin radius `ρ0` from the two
    /// standard parallels (Snyder pp. 101–102).
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        a: f64,
        e: f64,
        lat1: f64,
        lat2: f64,
        lat0: f64,
        lon0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        let (phi1, phi2) = (lat1.to_radians(), lat2.to_radians());
        let (m1, m2) = (m_factor(phi1, e), m_factor(phi2, e));
        let (q0, q1, q2) = (
            authalic_q(lat0.to_radians(), e),
            authalic_q(phi1, e),
            authalic_q(phi2, e),
        );
        let n = if (phi1 - phi2).abs() < 1e-12 {
            phi1.sin()
        } else {
            (m1 * m1 - m2 * m2) / (q2 - q1)
        };
        let cc = m1 * m1 + n * q1;
        let rho0 = a * (cc - n * q0).sqrt() / n;
        Self {
            a,
            e,
            n,
            cc,
            rho0,
            lon0,
            false_easting,
            false_northing,
        }
    }
}

/// Albers Equal-Area forward: geographic lat/lon **degrees** → `(x, y)` m.
#[must_use]
pub fn albers_forward<T: LazyReal>(lat_deg: T, lon_deg: T, c: &AlbersConstants) -> (T, T) {
    let e = lat_deg.scalar_like(c.e);
    let q = authalic_q(lat_deg.clone().to_radians(), e);
    let rho = lat_deg.scalar_like(c.a)
        * (lat_deg.scalar_like(c.cc) - lat_deg.scalar_like(c.n) * q).sqrt()
        / lat_deg.scalar_like(c.n);
    let theta = (lon_deg - lat_deg.scalar_like(c.lon0)).to_radians() * lat_deg.scalar_like(c.n);
    let (st, ct) = theta.sin_cos();
    let x = lat_deg.scalar_like(c.false_easting) + rho.clone() * st;
    let y = lat_deg.scalar_like(c.false_northing + c.rho0) - rho * ct;
    (x, y)
}

/// Albers Equal-Area inverse: `(x, y)` m → geographic lat/lon **degrees**.
#[must_use]
pub fn albers_inverse<T: LazyReal>(x_in: T, y_in: T, c: &AlbersConstants) -> (T, T) {
    let e = x_in.scalar_like(c.e);
    let x = x_in.clone() - x_in.scalar_like(c.false_easting);
    let y0 = x_in.scalar_like(c.rho0) - (y_in - x_in.scalar_like(c.false_northing));
    let nsign = c.n.signum();
    let rho = x.clone().hypot(y0.clone()).copysign(x_in.scalar_like(c.n));
    let theta = (x.clone() * x_in.scalar_like(nsign)).atan2(y0 * x_in.scalar_like(nsign));
    let rho_na = rho * x_in.scalar_like(c.n) / x_in.scalar_like(c.a);
    let q = (x_in.scalar_like(c.cc) - rho_na.clone() * rho_na) / x_in.scalar_like(c.n);
    let lat = lat_from_authalic_q(q, e).to_degrees();
    let lon = (theta / x_in.scalar_like(c.n)).to_degrees() + x_in.scalar_like(c.lon0);
    (lat, lon)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual {actual} != expected {expected} within {tolerance}"
        );
    }

    const WGS84_A: f64 = 6_378_137.0;
    const WGS84_F: f64 = 1.0 / 298.257_223_563;

    /// The Transverse Mercator parameters of UTM zone `zone` (north).
    fn utm(zone: u8) -> (f64, KrugerSeries) {
        let central_meridian = f64::from(i32::from(zone) * 6 - 183);
        (
            central_meridian,
            KrugerSeries::for_ellipsoid(WGS84_A, WGS84_F),
        )
    }

    #[test]
    fn transverse_mercator_central_meridian_equator_is_false_origin() {
        // On the central meridian at the equator (UTM lat_origin = 0) the
        // projection lands exactly on the false origin: (500000, 0).
        let (cm, ks) = utm(31); // zone 31N → central meridian 3°E
        let (e, n) = transverse_mercator_forward(0.0_f64, cm, cm, 0.0, 0.9996, 500_000.0, 0.0, &ks);
        assert_close(e, 500_000.0, 1e-6);
        assert_close(n, 0.0, 1e-6);
    }

    #[test]
    fn transverse_mercator_round_trips() {
        let (cm, ks) = utm(31);
        for &(lat, lon) in &[(0.0, 3.0), (48.5, 2.0), (-33.9, 4.5), (60.0, 1.0)] {
            let (e, n) =
                transverse_mercator_forward(lat, lon, cm, 0.0, 0.9996, 500_000.0, 0.0, &ks);
            let (blat, blon) =
                transverse_mercator_inverse(e, n, cm, 0.0, 0.9996, 500_000.0, 0.0, &ks);
            assert_close(blat, lat, 1e-9);
            assert_close(blon, lon, 1e-9);
        }
    }

    #[test]
    fn web_mercator_round_trips_and_matches_closed_form() {
        let r = WGS84_A;
        // y at the equator is 0; x is R·λ.
        let (x, y) = web_mercator_forward(0.0_f64, 45.0, r);
        assert_close(x, r * 45.0_f64.to_radians(), 1e-6);
        assert_close(y, 0.0, 1e-9);
        for &(lat, lon) in &[(0.0, 45.0), (40.0, -75.0), (-12.5, 130.0)] {
            let (x, y) = web_mercator_forward(lat, lon, r);
            let (blat, blon) = web_mercator_inverse(x, y, r);
            assert_close(blat, lat, 1e-9);
            assert_close(blon, lon, 1e-9);
        }
    }

    #[test]
    fn projections_run_for_f32() {
        // Genericity over `LazyReal`: the same kernels run at `f32` (hence on
        // any lazy backend — `f64`/`f32`/tensor).
        let ks32 = KrugerSeries::for_ellipsoid(WGS84_A, WGS84_F);
        let (e, n) =
            transverse_mercator_forward(48.5_f32, 2.0_f32, 3.0, 0.0, 0.9996, 500_000.0, 0.0, &ks32);
        assert!(e.is_finite() && n.is_finite());
        let (x, y) = web_mercator_forward(40.0_f32, -75.0_f32, WGS84_A as f32);
        assert!(x.is_finite() && y.is_finite());
    }

    #[test]
    fn lat_lon_entry_matches_scalar() {
        let r = WGS84_A;
        let p = (40.0_f64, -75.0_f64); // (lat, lon) implements LatLon
        let viatuple = web_mercator_project_at(&p, r);
        let direct = web_mercator_forward(40.0_f64, -75.0_f64, r);
        assert_close(viatuple.0, direct.0, 1e-12);
        assert_close(viatuple.1, direct.1, 1e-12);
    }

    // Clarke 1866 — the ellipsoid of Snyder's worked numerical examples.
    const CLARKE_A: f64 = 6_378_206.4;
    const CLARKE_F: f64 = 1.0 / 294.978_698_2;

    fn ecc(f: f64) -> f64 {
        (f * (2.0 - f)).sqrt()
    }

    #[test]
    fn lcc_matches_snyder_and_round_trips() {
        // Snyder PP1395 p. 295: parallels 33°/45°, origin (23°, −96°),
        // point (35°, −75°) → (1 894 410.9, 1 564 649.5) m.
        let c = LccConstants::new(CLARKE_A, ecc(CLARKE_F), 33.0, 45.0, 23.0, -96.0, 0.0, 0.0);
        let (x, y) = lcc_forward(35.0_f64, -75.0, &c);
        assert_close(x, 1_894_410.9, 0.2);
        assert_close(y, 1_564_649.5, 0.2);
        let (lat, lon) = lcc_inverse(x, y, &c);
        assert_close(lat, 35.0, 1e-9);
        assert_close(lon, -75.0, 1e-9);
    }

    #[test]
    fn albers_matches_snyder_and_round_trips() {
        // Snyder PP1395 p. 292: parallels 29.5°/45.5°, origin (23°, −96°),
        // point (35°, −75°) → (1 885 472.7, 1 535 925.0) m.
        let c = AlbersConstants::new(CLARKE_A, ecc(CLARKE_F), 29.5, 45.5, 23.0, -96.0, 0.0, 0.0);
        let (x, y) = albers_forward(35.0_f64, -75.0, &c);
        assert_close(x, 1_885_472.7, 0.2);
        assert_close(y, 1_535_925.0, 0.2);
        let (lat, lon) = albers_inverse(x, y, &c);
        assert_close(lat, 35.0, 1e-9);
        assert_close(lon, -75.0, 1e-9);
    }

    #[test]
    fn polar_stereographic_round_trips() {
        let c = PolarStereoConstants::new(
            WGS84_A,
            ecc(WGS84_F),
            71.0,
            0.0,
            true,
            2_000_000.0,
            2_000_000.0,
        );
        for &(lat, lon) in &[(85.0_f64, 30.0), (72.0, -120.0), (89.0, 175.0)] {
            let (x, y) = polar_stereographic_forward(lat, lon, &c);
            let (blat, blon) = polar_stereographic_inverse(x, y, &c);
            assert_close(blat, lat, 1e-7);
            assert_close(blon, lon, 1e-6);
        }
    }

    #[test]
    fn auxiliary_latitudes_round_trip() {
        let e = ecc(WGS84_F);
        for &lat in &[0.0_f64, 12.3, 45.0, -60.0, 80.0] {
            let phi = lat.to_radians();
            assert_close(
                lat_from_conformal_t(conformal_t(phi, e), e).to_degrees(),
                lat,
                1e-9,
            );
            assert_close(
                lat_from_authalic_q(authalic_q(phi, e), e).to_degrees(),
                lat,
                1e-9,
            );
        }
    }

    #[test]
    fn families_run_for_f32() {
        // Genericity over `LazyReal`: the family kernels run at `f32` (constants
        // stay `f64`, the coordinate is the lazy scalar) — so they lower to any
        // lazy backend, tensors included.
        let e = ecc(WGS84_F);
        let lcc = LccConstants::new(WGS84_A, e, 33.0, 45.0, 39.0, -96.0, 0.0, 0.0);
        let (x, y) = lcc_forward(40.0_f32, -100.0_f32, &lcc);
        assert!(x.is_finite() && y.is_finite());
        let alb = AlbersConstants::new(WGS84_A, e, 29.5, 45.5, 23.0, -96.0, 0.0, 0.0);
        let (ax, ay) = albers_forward(35.0_f32, -100.0_f32, &alb);
        assert!(ax.is_finite() && ay.is_finite());
        let ps = PolarStereoConstants::new(WGS84_A, e, 71.0, 0.0, true, 0.0, 0.0);
        let (px, py) = polar_stereographic_forward(85.0_f32, 30.0_f32, &ps);
        assert!(px.is_finite() && py.is_finite());
        let e32 = e as f32;
        assert!(m_factor(0.7_f32, e32).is_finite());
        assert!(conformal_t(0.7_f32, e32).is_finite());
        assert!(authalic_q(0.7_f32, e32).is_finite());
    }
}
