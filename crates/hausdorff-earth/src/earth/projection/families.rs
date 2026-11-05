//! The conic, azimuthal, and local projection families — Lambert Conformal
//! Conic, Polar Stereographic, Albers Equal-Area, and the local tangent plane —
//! split out of the parent [`super`] projection module. The conformal/authalic
//! kernels come from `hausdorff-base`'s `earth::latitude` (beside the
//! `Ellipsoid::{conformal,authalic}_latitude` methods); `m_factor` comes from the
//! parent via `use super::*`.

use super::*;
use hausdorff_base::earth::latitude::conformal_t;
use hausdorff_kernels::projection::{
    AlbersConstants, LccConstants, PolarStereoConstants, albers_forward, albers_inverse,
    lcc_forward, lcc_inverse, polar_stereographic_forward, polar_stereographic_inverse,
};

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
    /// The Snyder cone constant `n`, scale `F`, and origin radius `ρ0` — solved by
    /// the backend-neutral [`LccConstants`] kernel (shared by the forward, inverse,
    /// and point-scale maps).
    fn kernel_constants(&self) -> LccConstants {
        LccConstants::new(
            self.ellipsoid.a,
            self.ellipsoid.e2().sqrt(),
            self.lat1,
            self.lat2,
            self.lat0,
            self.lon0,
            self.false_easting,
            self.false_northing,
        )
    }
}

impl Projection for LambertConformalConic {
    fn project(&self, p: GeoPoint) -> Coord {
        let (x, y) = lcc_forward(p.lat, p.lon, &self.kernel_constants());
        Coord::new(x, y)
    }

    /// Proof-backed by `riemann` `Hausdorff/ConicProjection`: the closed-form
    /// inverse recovers the radius `ρ` and cone angle `θ` exactly from the forward
    /// map (`conic_inv_rho`/`conic_inv_theta` via `hypot_rho`/`atan2_rho`), and
    /// `radius_recover` recovers the conformal parameter — so `inverse` is a true
    /// left inverse on the conic geometry (the ellipsoid latitude series aside).
    fn inverse(&self, c: Coord) -> GeoPoint {
        let (lat, lon) = lcc_inverse(c.x, c.y, &self.kernel_constants());
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

impl PolarStereographic {
    /// The standard-parallel constants — solved by the backend-neutral
    /// [`PolarStereoConstants`] kernel.
    fn kernel_constants(&self) -> PolarStereoConstants {
        PolarStereoConstants::new(
            self.ellipsoid.a,
            self.ellipsoid.e2().sqrt(),
            self.lat_ts,
            self.lon0,
            self.north,
            self.false_easting,
            self.false_northing,
        )
    }
}

impl Projection for PolarStereographic {
    fn project(&self, p: GeoPoint) -> Coord {
        let (x, y) = polar_stereographic_forward(p.lat, p.lon, &self.kernel_constants());
        Coord::new(x, y)
    }

    fn inverse(&self, c: Coord) -> GeoPoint {
        let (lat, lon) = polar_stereographic_inverse(c.x, c.y, &self.kernel_constants());
        GeoPoint::from_lat_lon(lat, lon)
    }
}

/// **Albers Equal-Area Conic** (ellipsoidal) — the equal-area conic of choropleth
/// and statistical maps (Snyder pp. 101–102), pinned to two standard parallels.
/// EPSG:5070 (Albers CONUS) is an instance.
///
/// **Equal-area is proof-backed** by `riemann` `Hausdorff/Authalic`: the cone
/// `ρ² = (C − 2n·q)/n²` is affine in the authalic area `q` (`albersRhoSq_hasDerivAt`),
/// so the plane area element `(n/2)·d(ρ²)` equals the spheroid's zonal density
/// `q′(φ)` (`albers_area_scale`) — area scale ≡ 1, exact (not just the numeric
/// graticule check).
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
    /// `(n, C, ρ0)` — the cone constant, the constant `C`, and the origin radius,
    /// solved by the backend-neutral [`AlbersConstants`] kernel.
    fn kernel_constants(&self) -> AlbersConstants {
        AlbersConstants::new(
            self.ellipsoid.a,
            self.ellipsoid.e2().sqrt(),
            self.lat1,
            self.lat2,
            self.lat0,
            self.lon0,
            self.false_easting,
            self.false_northing,
        )
    }
}

impl Projection for AlbersEqualArea {
    fn project(&self, p: GeoPoint) -> Coord {
        let (x, y) = albers_forward(p.lat, p.lon, &self.kernel_constants());
        Coord::new(x, y)
    }

    fn inverse(&self, c: Coord) -> GeoPoint {
        let (lat, lon) = albers_inverse(c.x, c.y, &self.kernel_constants());
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
/// ([`EnuFrame`]): the `up` axis is the geodetic
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
        let k = self.kernel_constants();
        let phi = lat_deg.to_radians();
        let rho = k.a * k.f * conformal_t(phi, k.e).powf(k.n);
        rho * k.n / (k.a * m_factor(phi, k.e))
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
