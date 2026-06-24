use super::*;

// ---- Conic / azimuthal / equal-area projection families ----
//
// Each family's per-projection constants (cone constant, scale, origin radius)
// depend only on the ellipsoid and standard parallels, so they are computed once
// in `f64` (a `*Constants` struct, the peer of `KrugerSeries`) and can be lifted
// into a typed constants view once per backend. The f64 convenience functions
// still lift internally; the `_typed` entry points let tensor/expression callers
// reuse pre-lifted constants across many coordinates.

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

    /// Lift these f64 projection parameters into a backend scalar once, so
    /// repeated point projections can reuse the same typed constants.
    #[must_use]
    pub fn lift<T: Real>(&self) -> LccConstantsT<T> {
        LccConstantsT {
            a_f: T::scalar(self.a * self.f),
            e: T::scalar(self.e),
            n: T::scalar(self.n),
            inv_n: T::scalar(1.0 / self.n),
            rho0: T::scalar(self.rho0),
            lon0: T::scalar(self.lon0),
            false_easting: T::scalar(self.false_easting),
            false_northing: T::scalar(self.false_northing),
            false_northing_plus_rho0: T::scalar(self.false_northing + self.rho0),
            nsign: T::scalar(self.n.signum()),
        }
    }
}

/// Backend-typed Lambert Conformal Conic constants.
///
/// Construct with [`LccConstants::lift`] from the f64 projection definition, or
/// build directly when a backend owns its constants natively.
#[derive(Clone, Debug)]
pub struct LccConstantsT<T> {
    pub a_f: T,
    pub e: T,
    pub n: T,
    pub inv_n: T,
    pub rho0: T,
    pub lon0: T,
    pub false_easting: T,
    pub false_northing: T,
    pub false_northing_plus_rho0: T,
    pub nsign: T,
}

/// Lambert Conformal Conic forward: geographic lat/lon **degrees** → `(x, y)` m.
#[must_use]
pub fn lcc_forward<T: Real>(lat_deg: T, lon_deg: T, c: &LccConstants) -> (T, T) {
    let c = c.lift();
    lcc_forward_typed(lat_deg, lon_deg, &c)
}

/// Lambert Conformal Conic forward with already-typed constants.
#[must_use]
pub fn lcc_forward_typed<T: Real>(lat_deg: T, lon_deg: T, c: &LccConstantsT<T>) -> (T, T) {
    let t = conformal_t(lat_deg.clone().to_radians(), c.e.clone());
    let rho = c.a_f.clone() * t.powf(c.n.clone());
    let theta = (lon_deg - c.lon0.clone()).to_radians() * c.n.clone();
    let (st, ct) = theta.sin_cos();
    let x = c.false_easting.clone() + rho.clone() * st;
    let y = c.false_northing_plus_rho0.clone() - rho * ct;
    (x, y)
}

/// Lambert Conformal Conic inverse: `(x, y)` m → geographic lat/lon **degrees**.
#[must_use]
pub fn lcc_inverse<T: Real>(x_in: T, y_in: T, c: &LccConstants) -> (T, T) {
    let c = c.lift();
    lcc_inverse_typed(x_in, y_in, &c)
}

/// Lambert Conformal Conic inverse with already-typed constants.
#[must_use]
pub fn lcc_inverse_typed<T: Real>(x_in: T, y_in: T, c: &LccConstantsT<T>) -> (T, T) {
    let x = x_in.clone() - c.false_easting.clone();
    let y0 = c.rho0.clone() - (y_in - c.false_northing.clone());
    let rho = x.clone().hypot(y0.clone()).copysign(c.n.clone());
    let theta = (x.clone() * c.nsign.clone()).atan2(y0 * c.nsign.clone());
    let t = (rho / c.a_f.clone()).powf(c.inv_n.clone());
    let lat = lat_from_conformal_t(t, c.e.clone()).to_degrees();
    let lon = (theta / c.n.clone()).to_degrees() + c.lon0.clone();
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

    /// Lift these f64 projection parameters into a backend scalar once, so
    /// repeated point projections can reuse the same typed constants.
    #[must_use]
    pub fn lift<T: Real>(&self) -> PolarStereoConstantsT<T> {
        PolarStereoConstantsT {
            e: T::scalar(self.e),
            s: T::scalar(self.s),
            neg_s: T::scalar(-self.s),
            a_mc_over_tc: T::scalar(self.a * self.mc / self.tc),
            tc_over_a_mc: T::scalar(self.tc / (self.a * self.mc)),
            lon0: T::scalar(self.lon0),
            false_easting: T::scalar(self.false_easting),
            false_northing: T::scalar(self.false_northing),
        }
    }
}

/// Backend-typed Polar Stereographic constants.
#[derive(Clone, Debug)]
pub struct PolarStereoConstantsT<T> {
    pub e: T,
    pub s: T,
    pub neg_s: T,
    pub a_mc_over_tc: T,
    pub tc_over_a_mc: T,
    pub lon0: T,
    pub false_easting: T,
    pub false_northing: T,
}

/// Polar Stereographic forward: geographic lat/lon **degrees** → `(x, y)` m.
#[must_use]
pub fn polar_stereographic_forward<T: Real>(
    lat_deg: T,
    lon_deg: T,
    c: &PolarStereoConstants,
) -> (T, T) {
    let c = c.lift();
    polar_stereographic_forward_typed(lat_deg, lon_deg, &c)
}

/// Polar Stereographic forward with already-typed constants.
#[must_use]
pub fn polar_stereographic_forward_typed<T: Real>(
    lat_deg: T,
    lon_deg: T,
    c: &PolarStereoConstantsT<T>,
) -> (T, T) {
    let phi = lat_deg.clone().to_radians() * c.s.clone();
    let rho = c.a_mc_over_tc.clone() * conformal_t(phi, c.e.clone());
    let dlon = (lon_deg - c.lon0.clone()).to_radians();
    let (sdl, cdl) = dlon.sin_cos();
    let x = c.false_easting.clone() + rho.clone() * sdl;
    let y = c.false_northing.clone() - rho * cdl * c.s.clone();
    (x, y)
}

/// Polar Stereographic inverse: `(x, y)` m → geographic lat/lon **degrees**.
#[must_use]
pub fn polar_stereographic_inverse<T: Real>(x_in: T, y_in: T, c: &PolarStereoConstants) -> (T, T) {
    let c = c.lift();
    polar_stereographic_inverse_typed(x_in, y_in, &c)
}

/// Polar Stereographic inverse with already-typed constants.
#[must_use]
pub fn polar_stereographic_inverse_typed<T: Real>(
    x_in: T,
    y_in: T,
    c: &PolarStereoConstantsT<T>,
) -> (T, T) {
    let x = x_in.clone() - c.false_easting.clone();
    let y = y_in - c.false_northing.clone();
    let rho = x.clone().hypot(y.clone());
    let t = rho * c.tc_over_a_mc.clone();
    let lat = (lat_from_conformal_t(t, c.e.clone()) * c.s.clone()).to_degrees();
    let lon = c.lon0.clone() + x.atan2(y * c.neg_s.clone()).to_degrees();
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

    /// Lift these f64 projection parameters into a backend scalar once, so
    /// repeated point projections can reuse the same typed constants.
    #[must_use]
    pub fn lift<T: Real>(&self) -> AlbersConstantsT<T> {
        AlbersConstantsT {
            a: T::scalar(self.a),
            e: T::scalar(self.e),
            n: T::scalar(self.n),
            cc: T::scalar(self.cc),
            rho0: T::scalar(self.rho0),
            lon0: T::scalar(self.lon0),
            false_easting: T::scalar(self.false_easting),
            false_northing: T::scalar(self.false_northing),
            false_northing_plus_rho0: T::scalar(self.false_northing + self.rho0),
            nsign: T::scalar(self.n.signum()),
        }
    }
}

/// Backend-typed Albers Equal-Area constants.
#[derive(Clone, Debug)]
pub struct AlbersConstantsT<T> {
    pub a: T,
    pub e: T,
    pub n: T,
    pub cc: T,
    pub rho0: T,
    pub lon0: T,
    pub false_easting: T,
    pub false_northing: T,
    pub false_northing_plus_rho0: T,
    pub nsign: T,
}

/// Albers Equal-Area forward: geographic lat/lon **degrees** → `(x, y)` m.
#[must_use]
pub fn albers_forward<T: Real>(lat_deg: T, lon_deg: T, c: &AlbersConstants) -> (T, T) {
    let c = c.lift();
    albers_forward_typed(lat_deg, lon_deg, &c)
}

/// Albers Equal-Area forward with already-typed constants.
#[must_use]
pub fn albers_forward_typed<T: Real>(lat_deg: T, lon_deg: T, c: &AlbersConstantsT<T>) -> (T, T) {
    let q = authalic_q(lat_deg.clone().to_radians(), c.e.clone());
    let rho = c.a.clone() * (c.cc.clone() - c.n.clone() * q).sqrt() / c.n.clone();
    let theta = (lon_deg - c.lon0.clone()).to_radians() * c.n.clone();
    let (st, ct) = theta.sin_cos();
    let x = c.false_easting.clone() + rho.clone() * st;
    let y = c.false_northing_plus_rho0.clone() - rho * ct;
    (x, y)
}

/// Albers Equal-Area inverse: `(x, y)` m → geographic lat/lon **degrees**.
#[must_use]
pub fn albers_inverse<T: Real>(x_in: T, y_in: T, c: &AlbersConstants) -> (T, T) {
    let c = c.lift();
    albers_inverse_typed(x_in, y_in, &c)
}

/// Albers Equal-Area inverse with already-typed constants.
#[must_use]
pub fn albers_inverse_typed<T: Real>(x_in: T, y_in: T, c: &AlbersConstantsT<T>) -> (T, T) {
    let x = x_in.clone() - c.false_easting.clone();
    let y0 = c.rho0.clone() - (y_in - c.false_northing.clone());
    let rho = x.clone().hypot(y0.clone()).copysign(c.n.clone());
    let theta = (x.clone() * c.nsign.clone()).atan2(y0 * c.nsign.clone());
    let rho_na = rho * c.n.clone() / c.a.clone();
    let q = (c.cc.clone() - rho_na.clone() * rho_na) / c.n.clone();
    let lat = lat_from_authalic_q(q, c.e.clone()).to_degrees();
    let lon = (theta / c.n.clone()).to_degrees() + c.lon0.clone();
    (lat, lon)
}
