//! The conic, azimuthal, and local projection families — Lambert Conformal
//! Conic, Polar Stereographic, Albers Equal-Area, and the local tangent plane —
//! split out of the parent [`super`] projection module. They share the parent's
//! conformal/authalic helpers (`conformal_t`, `m_factor`, …) via `use super::*`.

use super::*;

/// **Lambert Conformal Conic (2SP)** — the conformal conic of national/regional
/// grids (Snyder pp. 107–109), pinned to two standard parallels where scale is
/// true. Parameterized by ellipsoid, the two standard parallels, the false
/// origin (latitude/longitude), and the false easting/northing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LambertConformalConic {
    pub ellipsoid: Ellipsoid,
    /// First standard parallel, degrees.
    pub lat1: f64,
    /// Second standard parallel, degrees.
    pub lat2: f64,
    /// Latitude of the false origin, degrees.
    pub lat0: f64,
    /// Central meridian (longitude of origin), degrees.
    pub lon0: f64,
    pub false_easting: f64,
    pub false_northing: f64,
}

impl LambertConformalConic {
    /// The cone constant `n`, scale `F`, and origin radius `ρ0` shared by the
    /// forward and inverse maps.
    fn constants(&self) -> (f64, f64, f64, f64) {
        let e = self.ellipsoid.e2().sqrt();
        let (phi1, phi2) = (self.lat1.to_radians(), self.lat2.to_radians());
        let (m1, m2) = (m_factor(phi1, e), m_factor(phi2, e));
        let (t0, t1, t2) = (
            conformal_t(self.lat0.to_radians(), e),
            conformal_t(phi1, e),
            conformal_t(phi2, e),
        );
        let n = if (phi1 - phi2).abs() < 1e-12 {
            phi1.sin()
        } else {
            (m1.ln() - m2.ln()) / (t1.ln() - t2.ln())
        };
        let f = m1 / (n * t1.powf(n));
        let rho0 = self.ellipsoid.a * f * t0.powf(n);
        (e, n, f, rho0)
    }
}

impl Projection for LambertConformalConic {
    fn project(&self, p: GeoPoint) -> Coord {
        let (e, n, f, rho0) = self.constants();
        let t = conformal_t(p.lat.to_radians(), e);
        let rho = self.ellipsoid.a * f * t.powf(n);
        let theta = n * (p.lon - self.lon0).to_radians();
        Coord::new(
            self.false_easting + rho * theta.sin(),
            self.false_northing + rho0 - rho * theta.cos(),
        )
    }

    fn inverse(&self, c: Coord) -> GeoPoint {
        let (e, n, f, rho0) = self.constants();
        let x = c.x - self.false_easting;
        let y0 = rho0 - (c.y - self.false_northing);
        let rho = (x * x + y0 * y0).sqrt().copysign(n);
        let theta = (x * n.signum()).atan2(y0 * n.signum());
        let t = (rho / (self.ellipsoid.a * f)).powf(1.0 / n);
        let lat = lat_from_conformal_t(t, e).to_degrees();
        let lon = (theta / n).to_degrees() + self.lon0;
        GeoPoint::from_lat_lon(lat, lon)
    }
}

/// **Polar Stereographic** (ellipsoidal) — the conformal azimuthal projection of
/// the polar regions (Snyder pp. 160–161), pinned to a *latitude of true scale*
/// `lat_ts` (the standard parallel). `north` selects the aspect (north or south
/// pole). EPSG:3995 (Arctic) and EPSG:3031 (Antarctic) are instances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarStereographic {
    pub ellipsoid: Ellipsoid,
    /// Latitude of true scale (standard parallel), degrees.
    pub lat_ts: f64,
    /// Longitude straight down from the pole on the map, degrees.
    pub lon0: f64,
    /// North-pole aspect if `true`, south-pole if `false`.
    pub north: bool,
    pub false_easting: f64,
    pub false_northing: f64,
}

impl Projection for PolarStereographic {
    fn project(&self, p: GeoPoint) -> Coord {
        let e = self.ellipsoid.e2().sqrt();
        let s = if self.north { 1.0 } else { -1.0 };
        let phi = s * p.lat.to_radians();
        let phic = s * self.lat_ts.to_radians();
        let (tc, mc) = (conformal_t(phic, e), m_factor(phic, e));
        let rho = self.ellipsoid.a * mc * conformal_t(phi, e) / tc;
        let dlon = (p.lon - self.lon0).to_radians();
        Coord::new(
            self.false_easting + rho * dlon.sin(),
            self.false_northing - s * rho * dlon.cos(),
        )
    }

    fn inverse(&self, c: Coord) -> GeoPoint {
        let e = self.ellipsoid.e2().sqrt();
        let s = if self.north { 1.0 } else { -1.0 };
        let phic = s * self.lat_ts.to_radians();
        let (tc, mc) = (conformal_t(phic, e), m_factor(phic, e));
        let x = c.x - self.false_easting;
        let y = c.y - self.false_northing;
        let rho = x.hypot(y);
        let t = rho * tc / (self.ellipsoid.a * mc);
        let lat = (s * lat_from_conformal_t(t, e)).to_degrees();
        let lon = self.lon0 + x.atan2(-s * y).to_degrees();
        GeoPoint::from_lat_lon(lat, lon)
    }
}

/// The authalic-area function `q` (Snyder 3-12) — the equal-area analogue of the
/// conformal `t`; the integrand behind the Albers / cylindrical-equal-area maps.
fn authalic_q(phi: f64, e: f64) -> f64 {
    let es = e * phi.sin();
    (1.0 - e * e) * (phi.sin() / (1.0 - es * es) - (1.0 / (2.0 * e)) * ((1.0 - es) / (1.0 + es)).ln())
}

/// Invert the authalic `q` back to geodetic latitude by Snyder's iteration (3-16).
fn lat_from_authalic_q(q: f64, e: f64) -> f64 {
    let q_pole = authalic_q(std::f64::consts::FRAC_PI_2, e);
    if q.abs() >= q_pole - 1e-12 {
        return std::f64::consts::FRAC_PI_2.copysign(q);
    }
    let mut phi = (q / 2.0).asin();
    for _ in 0..15 {
        let es = e * phi.sin();
        let delta = (1.0 - es * es).powi(2) / (2.0 * phi.cos())
            * (q / (1.0 - e * e) - phi.sin() / (1.0 - es * es)
                + (1.0 / (2.0 * e)) * ((1.0 - es) / (1.0 + es)).ln());
        phi += delta;
        if delta.abs() < 1e-14 {
            break;
        }
    }
    phi
}

/// **Albers Equal-Area Conic** (ellipsoidal) — the equal-area conic of choropleth
/// and statistical maps (Snyder pp. 101–102), pinned to two standard parallels.
/// EPSG:5070 (Albers CONUS) is an instance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlbersEqualArea {
    pub ellipsoid: Ellipsoid,
    pub lat1: f64,
    pub lat2: f64,
    pub lat0: f64,
    pub lon0: f64,
    pub false_easting: f64,
    pub false_northing: f64,
}

impl AlbersEqualArea {
    /// `(e, n, C, ρ0)` — the cone constant, the constant `C`, and the origin
    /// radius shared by the forward and inverse maps.
    fn constants(&self) -> (f64, f64, f64, f64) {
        let e = self.ellipsoid.e2().sqrt();
        let (phi1, phi2) = (self.lat1.to_radians(), self.lat2.to_radians());
        let (m1, m2) = (m_factor(phi1, e), m_factor(phi2, e));
        let (q0, q1, q2) = (
            authalic_q(self.lat0.to_radians(), e),
            authalic_q(phi1, e),
            authalic_q(phi2, e),
        );
        let n = if (phi1 - phi2).abs() < 1e-12 {
            phi1.sin()
        } else {
            (m1 * m1 - m2 * m2) / (q2 - q1)
        };
        let cc = m1 * m1 + n * q1;
        let rho0 = self.ellipsoid.a * (cc - n * q0).sqrt() / n;
        (e, n, cc, rho0)
    }
}

impl Projection for AlbersEqualArea {
    fn project(&self, p: GeoPoint) -> Coord {
        let (e, n, cc, rho0) = self.constants();
        let q = authalic_q(p.lat.to_radians(), e);
        let rho = self.ellipsoid.a * (cc - n * q).sqrt() / n;
        let theta = n * (p.lon - self.lon0).to_radians();
        Coord::new(
            self.false_easting + rho * theta.sin(),
            self.false_northing + rho0 - rho * theta.cos(),
        )
    }

    fn inverse(&self, c: Coord) -> GeoPoint {
        let (e, n, cc, rho0) = self.constants();
        let x = c.x - self.false_easting;
        let y0 = rho0 - (c.y - self.false_northing);
        let rho = x.hypot(y0).copysign(n);
        let theta = (x * n.signum()).atan2(y0 * n.signum());
        let q = (cc - (rho * n / self.ellipsoid.a).powi(2)) / n;
        let lat = lat_from_authalic_q(q, e).to_degrees();
        let lon = (theta / n).to_degrees() + self.lon0;
        GeoPoint::from_lat_lon(lat, lon)
    }
}

/// A **local tangent-plane** (topocentric ENU) projection — the standard local
/// engineering / site-level coordinate system. Anchored at a geodetic origin, it
/// maps a geographic point to its **East/North displacement in metres** on the
/// plane tangent at the origin (the Up component is dropped). Exact at the
/// origin; the planar distortion grows with distance, so it is for local-area
/// work (surveys, construction sites, robotics) rather than wide regions.
///
/// Unlike the global map projections this is **parameterized by a local origin**,
/// not a standard grid. It is built on the proof-backed ENU frame
/// ([`EnuFrame`](hausdorff_base::earth::EnuFrame)): the `up` axis is the geodetic
/// surface normal and the east/north axes are orthonormal (`riemann`
/// `Hausdorff/Bearing`), and geodetic↔ECEF↔ENU are exact, so the projection
/// round-trips to machine precision.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalTangentPlane {
    frame: EnuFrame,
}

impl LocalTangentPlane {
    /// A local tangent plane anchored at `origin` on `ellipsoid` (origin at
    /// ellipsoidal height 0).
    #[must_use]
    pub fn new(ellipsoid: Ellipsoid, origin: GeoPoint) -> Self {
        Self {
            frame: ellipsoid.enu_frame(origin, 0.0),
        }
    }

    /// The anchored ENU frame (for full 3-D east/north/up access).
    #[must_use]
    pub fn frame(&self) -> EnuFrame {
        self.frame
    }
}

impl Projection for LocalTangentPlane {
    fn project(&self, p: GeoPoint) -> Coord {
        let enu = self.frame.to_enu(p, 0.0);
        Coord::new(enu.east, enu.north)
    }

    fn inverse(&self, c: Coord) -> GeoPoint {
        // Solve for the `up` that lands the point back on the ellipsoid surface
        // (height 0), so this is a true inverse of `project` on the surface, not
        // a drop onto the tangent plane. `up ≈ height` near the point, so the
        // Newton step `up -= h` converges in a few iterations.
        let mut up = 0.0;
        for _ in 0..8 {
            let (geo, h) = self.frame.from_enu(Enu {
                east: c.x,
                north: c.y,
                up,
            });
            if h.abs() < 1e-9 {
                return geo;
            }
            up -= h;
        }
        self.frame
            .from_enu(Enu {
                east: c.x,
                north: c.y,
                up,
            })
            .0
    }
}

impl ConformalProjection for LambertConformalConic {
    /// `k = n·ρ / (a·m(φ))` — unity at the two standard parallels (Snyder 14-16).
    fn point_scale(&self, lat_deg: f64) -> f64 {
        let (e, n, f, _rho0) = self.constants();
        let phi = lat_deg.to_radians();
        let rho = self.ellipsoid.a * f * conformal_t(phi, e).powf(n);
        rho * n / (self.ellipsoid.a * m_factor(phi, e))
    }
}

impl ConformalProjection for PolarStereographic {
    /// `k = mc·t(φ) / (tc·m(φ))` — unity at the latitude of true scale.
    fn point_scale(&self, lat_deg: f64) -> f64 {
        let e = self.ellipsoid.e2().sqrt();
        let s = if self.north { 1.0 } else { -1.0 };
        let phi = s * lat_deg.to_radians();
        let phic = s * self.lat_ts.to_radians();
        let (tc, mc) = (conformal_t(phic, e), m_factor(phic, e));
        mc * conformal_t(phi, e) / (tc * m_factor(phi, e))
    }
}
