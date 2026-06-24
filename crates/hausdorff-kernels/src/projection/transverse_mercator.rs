use super::*;

// ---- Transverse Mercator (Karney–Krüger, ellipsoidal) ----

/// Karney–Krüger Transverse Mercator series constants for one ellipsoid: the
/// rectifying radius `A`, the forward (`α`) and inverse (`β`) coefficient sets
/// (order `n⁶`), and the eccentricity `e`. They depend only on the ellipsoid,
/// so they are computed once in `f64` and can be lifted into a typed constants
/// view once per backend.
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

    /// Lift the ellipsoid and projection parameters into a backend scalar once,
    /// so repeated point projections can reuse the same typed constants.
    #[must_use]
    pub fn lift<T: Real>(
        &self,
        central_meridian: f64,
        lat_origin: f64,
        k0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> KrugerSeriesT<T> {
        KrugerSeriesT {
            e: T::scalar(self.e),
            alpha: self.alpha.map(T::scalar),
            beta: self.beta.map(T::scalar),
            two_j: core::array::from_fn(|j| T::scalar(2.0 * (j as f64 + 1.0))),
            scale: T::scalar(k0 * self.a_rect),
            xi_origin: T::scalar(self.xi_origin(lat_origin)),
            central_meridian: T::scalar(central_meridian),
            false_easting: T::scalar(false_easting),
            false_northing: T::scalar(false_northing),
        }
    }
}

/// Backend-typed Transverse Mercator constants.
///
/// This bundles the lifted Krüger series with the projection origin and false
/// origin. Construct with [`KrugerSeries::lift`] from the f64 definition, or
/// build directly when a backend owns its constants natively.
#[derive(Clone, Debug)]
pub struct KrugerSeriesT<T> {
    pub e: T,
    pub alpha: [T; 6],
    pub beta: [T; 6],
    pub two_j: [T; 6],
    pub scale: T,
    pub xi_origin: T,
    pub central_meridian: T,
    pub false_easting: T,
    pub false_northing: T,
}

/// Transverse Mercator forward: geographic lat/lon **degrees** to projected
/// `(easting, northing)` metres, via the Karney–Krüger series. `series` carries
/// the ellipsoid constants; the remaining `f64` parameters are the central
/// meridian, latitude of origin, central-meridian scale `k0`, and false origin.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn transverse_mercator_forward<T: Real>(
    lat_deg: T,
    lon_deg: T,
    central_meridian: f64,
    lat_origin: f64,
    k0: f64,
    false_easting: f64,
    false_northing: f64,
    series: &KrugerSeries,
) -> (T, T) {
    let c = series.lift(
        central_meridian,
        lat_origin,
        k0,
        false_easting,
        false_northing,
    );
    transverse_mercator_forward_typed(lat_deg, lon_deg, &c)
}

/// Transverse Mercator forward with already-typed constants.
#[must_use]
pub fn transverse_mercator_forward_typed<T: Real>(
    lat_deg: T,
    lon_deg: T,
    c: &KrugerSeriesT<T>,
) -> (T, T) {
    // Conformal latitude (as its tangent τ′), then the Gauss–Schreiber
    // intermediate (ξ′, η′) about the central meridian.
    let lam = (lon_deg - c.central_meridian.clone()).to_radians();
    let (slam, clam) = lam.sin_cos();
    let taup = taupf(lat_deg.clone().to_radians().tan(), c.e.clone());
    let xip = taup.clone().atan2(clam.clone());
    let etap = (slam / taup.hypot(clam)).asinh();
    // Krüger series ξ′,η′ → ξ,η.
    let (mut xi, mut eta) = (xip.clone(), etap.clone());
    for (a_j, two_j) in c.alpha.iter().cloned().zip(c.two_j.iter().cloned()) {
        let two_xip = two_j.clone() * xip.clone();
        let two_etap = two_j * etap.clone();
        xi = xi + a_j.clone() * two_xip.clone().sin() * two_etap.clone().cosh();
        eta = eta + a_j * two_xip.cos() * two_etap.sinh();
    }
    let easting = c.false_easting.clone() + c.scale.clone() * eta;
    let northing = c.false_northing.clone() + c.scale.clone() * (xi - c.xi_origin.clone());
    (easting, northing)
}

/// Transverse Mercator inverse: projected `(easting, northing)` metres back to
/// geographic lat/lon **degrees**. Inverse-pair to
/// [`transverse_mercator_forward`].
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn transverse_mercator_inverse<T: Real>(
    easting: T,
    northing: T,
    central_meridian: f64,
    lat_origin: f64,
    k0: f64,
    false_easting: f64,
    false_northing: f64,
    series: &KrugerSeries,
) -> (T, T) {
    let c = series.lift(
        central_meridian,
        lat_origin,
        k0,
        false_easting,
        false_northing,
    );
    transverse_mercator_inverse_typed(easting, northing, &c)
}

/// Transverse Mercator inverse with already-typed constants.
#[must_use]
pub fn transverse_mercator_inverse_typed<T: Real>(
    easting: T,
    northing: T,
    c: &KrugerSeriesT<T>,
) -> (T, T) {
    let xi = c.xi_origin.clone() + (northing - c.false_northing.clone()) / c.scale.clone();
    let eta = (easting.clone() - c.false_easting.clone()) / c.scale.clone();
    // Krüger series ξ,η → ξ′,η′.
    let (mut xip, mut etap) = (xi.clone(), eta.clone());
    for (b_j, two_j) in c.beta.iter().cloned().zip(c.two_j.iter().cloned()) {
        let two_xi = two_j.clone() * xi.clone();
        let two_eta = two_j * eta.clone();
        xip = xip - b_j.clone() * two_xi.clone().sin() * two_eta.clone().cosh();
        etap = etap - b_j * two_xi.cos() * two_eta.sinh();
    }
    // Invert the Gauss–Schreiber step: recover λ and the conformal latitude
    // (as τ′), then the geodetic latitude via `tauf`.
    let (sxip, cxip) = xip.sin_cos();
    let sinh_etap = etap.sinh();
    let lam = sinh_etap.clone().atan2(cxip.clone());
    let taup = sxip / sinh_etap.hypot(cxip);
    let lat_deg = tauf(taup, c.e.clone()).atan().to_degrees();
    let lon_deg = c.central_meridian.clone() + lam.to_degrees();
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
    T: Real,
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

/// Project a [`LatLon`] provider through [`transverse_mercator_forward_typed`].
#[must_use]
pub fn transverse_mercator_project_at_typed<T, P>(p: &P, c: &KrugerSeriesT<T>) -> (T, T)
where
    T: Real,
    P: LatLon<T>,
{
    transverse_mercator_forward_typed(p.lat_deg(), p.lon_deg(), c)
}
