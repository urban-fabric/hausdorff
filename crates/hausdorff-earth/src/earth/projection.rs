//! `projection` — interchange between the geographic/spherical world and
//! projected **planar** coordinates via standard map projections.
//!
//! A [`Projection`](crate::earth::projection::Projection) maps WGS84 geographic
//! [`GeoPoint`]s (lat/lon — the hub shared with the spherical [`hausdorff_sphere::space::s2`]
//! side through [`hausdorff_sphere::space::s2::S2Point::from_geo_point`]) to and from planar
//! [`Coord`](hausdorff_planar::space::r2::Coord)s in metres, so the `r2` planar algorithms
//! can run on real projected data and the results map back to the sphere.
//!
//! Projections are parameterized by an [`Ellipsoid`](crate::earth::projection::Ellipsoid)
//! (so they are datum- and planet-flexible, mirroring
//! [`crate::earth::sphere::Sphere`] on the metric side) and follow Snyder,
//! *Map Projections — A Working Manual* (USGS PP1395, 1987). Provided:
//! [`WebMercator`](crate::earth::projection::WebMercator) (EPSG:3857) and
//! [`TransverseMercator`](crate::earth::projection::TransverseMercator) (the basis
//! of [`Utm`](crate::earth::projection::Utm)). The trivial `x = lon`, `y = lat`
//! relabelling — *not* a projection — stays in [`crate::convert`].
//!
//! Whole geometries project vertex-by-vertex through
//! [`Projection::project_geometry`](crate::earth::projection::Projection::project_geometry) /
//! [`Projection::inverse_geometry`](crate::earth::projection::Projection::inverse_geometry)
//! (via [`MapCoords`](hausdorff_planar::space::r2::MapCoords)), so lon/lat `r2` geometry
//! can be carried into projected metres, computed on, and carried back.

use crate::earth::GeoPoint;
use hausdorff_base::earth::{Enu, EnuFrame};
use hausdorff_planar::space::r2::{Coord, MapCoords};

pub use crate::earth::ellipsoid::Ellipsoid;

mod families;
mod tissot;
pub use families::{AlbersEqualArea, LambertConformalConic, LocalTangentPlane, PolarStereographic};
pub use tissot::Tissot;

#[cfg(test)]
mod tests;

/// A map projection: geographic ↔ projected planar metres.
pub trait Projection {
    /// Geographic → planar metres.
    fn project(&self, p: GeoPoint) -> Coord;
    /// Planar metres → geographic.
    fn inverse(&self, c: Coord) -> GeoPoint;

    /// Project a lon/lat planar coordinate (`x = lon`, `y = lat` degrees — the
    /// [`crate::convert`] relabelling) to projected planar metres.
    fn project_coord(&self, c: Coord) -> Coord {
        self.project(GeoPoint::from_lon_lat(c.x, c.y))
    }

    /// Inverse-project planar metres back to a lon/lat planar coordinate
    /// (`x = lon`, `y = lat` degrees).
    fn inverse_coord(&self, c: Coord) -> Coord {
        let p = self.inverse(c);
        Coord::new(p.lon, p.lat)
    }

    /// Project a whole geometry, vertex by vertex, from lon/lat degrees
    /// (`x = lon`, `y = lat`) into projected planar metres — after which the
    /// `r2` algorithms compute in metres. Inverse-pair to
    /// [`inverse_geometry`](Self::inverse_geometry).
    ///
    /// Only vertices are mapped: an edge stays the straight segment between its
    /// projected endpoints, so a long edge is *not* densified along the
    /// projected curve (and a [`Rect`](hausdorff_planar::space::r2::Rect) maps by its two
    /// corners — see [`MapCoords`]).
    fn project_geometry<G: MapCoords>(&self, g: &G) -> G
    where
        Self: Sized,
    {
        g.map_coords(|c| self.project_coord(c))
    }

    /// Inverse-project a whole geometry, vertex by vertex, from projected
    /// planar metres back to lon/lat degrees (`x = lon`, `y = lat`).
    /// Inverse-pair to [`project_geometry`](Self::project_geometry).
    fn inverse_geometry<G: MapCoords>(&self, g: &G) -> G
    where
        Self: Sized,
    {
        g.map_coords(|c| self.inverse_coord(c))
    }
}

/// Web Mercator (spherical Mercator) — the projection of web map tiles. Planar
/// units are metres; valid for `|lat| ≲ 85.06°`. EPSG:3857 uses the WGS84
/// semi-major axis as the sphere radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WebMercator {
    /// Sphere radius, metres.
    pub radius: f64,
}

impl WebMercator {
    /// The standard web-tile projection (EPSG:3857).
    pub const EPSG_3857: WebMercator = WebMercator {
        radius: 6_378_137.0,
    };
}

impl Projection for WebMercator {
    fn project(&self, p: GeoPoint) -> Coord {
        // y = R·ψ(φ), ψ the isometric latitude `log((1 + sinφ)/cosφ)` — the
        // exact spherical-Mercator form, `riemann` `Hausdorff/Rhumb.isoLat`.
        let (s, c) = p.lat.to_radians().sin_cos();
        Coord::new(self.radius * p.lon.to_radians(), self.radius * ((1.0 + s) / c).ln())
    }

    fn inverse(&self, c: Coord) -> GeoPoint {
        // φ = arcsin(tanh(y/R)) — the **provably exact** inverse (`riemann`
        // `Hausdorff/Mercator.mercator_inverse`: forward∘inverse = id), since
        // `tanh(ψ(φ)) = sinφ` (`Rhumb.tanh_isoLat`).
        GeoPoint::from_lat_lon(
            (c.y / self.radius).tanh().asin().to_degrees(),
            (c.x / self.radius).to_degrees(),
        )
    }
}

/// Transverse Mercator — the cylindrical-transverse projection underlying UTM.
/// Parameterized by ellipsoid, central meridian, scale, and false origin, so it
/// also expresses national grids.
///
/// Implemented with the **Karney–Krüger** series to order `n⁶` (Karney 2011,
/// *Transverse Mercator with an accuracy of a few nanometers*) — accurate to
/// sub-millimetre across the **whole** zone, not just near the central meridian
/// (where the lower-order Snyder series degrades). Still a series, so it stays on
/// the numerically-validated side of the `f64 = ℝ` proof boundary (cf. the module
/// note); only the *spherical* limit is exact-provable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseMercator {
    pub ellipsoid: Ellipsoid,
    /// Central meridian, degrees.
    pub central_meridian: f64,
    /// Latitude of natural origin, degrees — the parallel northings are measured
    /// from. UTM uses the equator (`0.0`); national grids (e.g. the British
    /// National Grid at 49°N) use a non-zero origin.
    pub lat_origin: f64,
    /// Scale factor on the central meridian.
    pub k0: f64,
    pub false_easting: f64,
    pub false_northing: f64,
}

impl TransverseMercator {
    /// The number of UTM zones (1..=60).
    pub const UTM_ZONES: u8 = 60;

    /// The Transverse Mercator for a UTM `zone` (1..=60) on the WGS84 ellipsoid.
    ///
    /// In debug builds an out-of-range `zone` panics; prefer [`try_utm`](Self::try_utm)
    /// to validate untrusted input.
    #[must_use]
    pub fn utm(zone: u8, north: bool) -> Self {
        debug_assert!(
            (1..=Self::UTM_ZONES).contains(&zone),
            "UTM zone {zone} out of range 1..=60"
        );
        Self {
            ellipsoid: Ellipsoid::WGS84,
            central_meridian: f64::from(i32::from(zone) * 6 - 183),
            lat_origin: 0.0,
            k0: 0.9996,
            false_easting: 500_000.0,
            false_northing: if north { 0.0 } else { 10_000_000.0 },
        }
    }

    /// The Transverse Mercator for a UTM `zone`, or `None` if `zone` is not in
    /// `1..=60`. The checked form of [`utm`](Self::utm).
    #[must_use]
    pub fn try_utm(zone: u8, north: bool) -> Option<Self> {
        (1..=Self::UTM_ZONES)
            .contains(&zone)
            .then(|| Self::utm(zone, north))
    }
}

impl Projection for TransverseMercator {
    fn project(&self, p: GeoPoint) -> Coord {
        let ks = self.kruger();
        // Conformal latitude (as its tangent τ′), then the Gauss–Schreiber
        // intermediate (ξ′, η′) about the central meridian.
        let lam = (p.lon - self.central_meridian).to_radians();
        let (slam, clam) = lam.sin_cos();
        let taup = taupf(p.lat.to_radians().tan(), ks.e);
        let xip = taup.atan2(clam);
        let etap = (slam / taup.hypot(clam)).asinh();
        // Krüger series ξ′,η′ → ξ,η.
        let (mut xi, mut eta) = (xip, etap);
        for (j, &a_j) in ks.alpha.iter().enumerate() {
            let two_j = 2.0 * (j as f64 + 1.0);
            xi += a_j * (two_j * xip).sin() * (two_j * etap).cosh();
            eta += a_j * (two_j * xip).cos() * (two_j * etap).sinh();
        }
        let scale = self.k0 * ks.a_rect;
        Coord::new(
            self.false_easting + scale * eta,
            self.false_northing + scale * (xi - ks.xi_origin(self.lat_origin)),
        )
    }

    fn inverse(&self, coord: Coord) -> GeoPoint {
        let ks = self.kruger();
        let scale = self.k0 * ks.a_rect;
        let xi = ks.xi_origin(self.lat_origin) + (coord.y - self.false_northing) / scale;
        let eta = (coord.x - self.false_easting) / scale;
        // Krüger series ξ,η → ξ′,η′.
        let (mut xip, mut etap) = (xi, eta);
        for (j, &b_j) in ks.beta.iter().enumerate() {
            let two_j = 2.0 * (j as f64 + 1.0);
            xip -= b_j * (two_j * xi).sin() * (two_j * eta).cosh();
            etap -= b_j * (two_j * xi).cos() * (two_j * eta).sinh();
        }
        // Invert the Gauss–Schreiber step: recover λ and the conformal latitude
        // (as τ′), then the geodetic latitude via `tauf`.
        let (sxip, cxip) = xip.sin_cos();
        let sinh_etap = etap.sinh();
        let lam = sinh_etap.atan2(cxip);
        let taup = sxip / sinh_etap.hypot(cxip);
        let lat = tauf(taup, ks.e).atan();
        GeoPoint::from_lat_lon(
            lat.to_degrees(),
            self.central_meridian + lam.to_degrees(),
        )
    }
}

/// Karney–Krüger Transverse Mercator series constants for one ellipsoid: the
/// rectifying radius `A`, the forward (`α`) and inverse (`β`) coefficient sets
/// (order `n⁶`), and the eccentricity `e`. Depend only on the ellipsoid.
struct KrugerSeries {
    /// Rectifying radius `A = a/(1+n)·(1 + n²/4 + n⁴/64 + n⁶/256)`.
    a_rect: f64,
    /// Forward coefficients α₁..α₆.
    alpha: [f64; 6],
    /// Inverse coefficients β₁..β₆.
    beta: [f64; 6],
    /// First eccentricity.
    e: f64,
}

impl KrugerSeries {
    /// ξ at the natural origin (latitude `lat_origin`, on the central meridian) —
    /// the rectifying-latitude offset subtracted so northings are measured from
    /// the origin parallel. At λ = 0 the Gauss–Schreiber η′ is 0.
    fn xi_origin(&self, lat_origin: f64) -> f64 {
        if lat_origin == 0.0 {
            return 0.0;
        }
        let xip = taupf(lat_origin.to_radians().tan(), self.e).atan();
        let mut xi = xip;
        for (j, &a_j) in self.alpha.iter().enumerate() {
            xi += a_j * (2.0 * (j as f64 + 1.0) * xip).sin();
        }
        xi
    }
}

impl TransverseMercator {
    /// The Karney–Krüger series constants for this projection's ellipsoid.
    fn kruger(&self) -> KrugerSeries {
        let f = self.ellipsoid.f;
        let n = f / (2.0 - f); // third flattening
        let (n2, n3, n4, n5, n6) = (n * n, n * n * n, n.powi(4), n.powi(5), n.powi(6));
        let a_rect =
            self.ellipsoid.a / (1.0 + n) * (1.0 + n2 / 4.0 + n4 / 64.0 + n6 / 256.0);
        // Krüger forward (α) and inverse (β) coefficients, order n⁶ (Karney 2011).
        let alpha = [
            n / 2.0 - 2.0 / 3.0 * n2 + 5.0 / 16.0 * n3 + 41.0 / 180.0 * n4
                - 127.0 / 288.0 * n5
                + 7891.0 / 37800.0 * n6,
            13.0 / 48.0 * n2 - 3.0 / 5.0 * n3 + 557.0 / 1440.0 * n4 + 281.0 / 630.0 * n5
                - 1_983_433.0 / 1_935_360.0 * n6,
            61.0 / 240.0 * n3 - 103.0 / 140.0 * n4 + 15_061.0 / 26_880.0 * n5
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
        KrugerSeries {
            a_rect,
            alpha,
            beta,
            e: self.ellipsoid.e2().sqrt(),
        }
    }
}

/// A UTM coordinate: zone (1..=60), hemisphere, and easting/northing in metres
/// (WGS84). A thin convenience over [`TransverseMercator::utm`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Utm {
    pub zone: u8,
    pub north: bool,
    pub easting: f64,
    pub northing: f64,
}

impl Utm {
    /// The UTM zone (1..=60) for a longitude in degrees.
    #[must_use]
    pub fn zone_for_lon(lon: f64) -> u8 {
        (((lon + 180.0) / 6.0).floor() as i64).rem_euclid(60) as u8 + 1
    }

    /// The easting/northing as a planar [`Coord`].
    #[must_use]
    pub fn coord(self) -> Coord {
        Coord::new(self.easting, self.northing)
    }

    /// Project a geographic point to UTM (zone chosen from its longitude).
    #[must_use]
    pub fn from_geo(p: GeoPoint) -> Utm {
        let zone = Self::zone_for_lon(p.lon);
        let north = p.lat >= 0.0;
        let c = TransverseMercator::utm(zone, north).project(p);
        Utm {
            zone,
            north,
            easting: c.x,
            northing: c.y,
        }
    }

    /// Inverse-project this UTM coordinate back to geographic.
    #[must_use]
    pub fn to_geo(self) -> GeoPoint {
        TransverseMercator::utm(self.zone, self.north).inverse(self.coord())
    }
}

// ---- Shared ellipsoidal-projection helpers (Snyder PP1395) ----
//
// The auxiliary-latitude kernels (`conformal_t`/`lat_from_conformal_t`,
// `authalic_q`/`lat_from_authalic_q`, the meridian arc) live in
// `hausdorff-base`'s `earth::latitude` — the canonical home beside the
// `Ellipsoid::{conformal,authalic,rectifying}_latitude` methods they back — and
// are imported by [`families`] / used through the `Ellipsoid` methods here.

/// `m = cosφ / √(1 − e² sin²φ)` (Snyder 14-15) — the radius of the parallel
/// scaled to the unit sphere; the prime-vertical-radius factor.
fn m_factor(phi: f64, e: f64) -> f64 {
    let s = phi.sin();
    phi.cos() / (1.0 - e * e * s * s).sqrt()
}

/// `tan χ` from `tan φ`: the tangent of the **conformal latitude** χ given the
/// tangent of the geodetic latitude φ, on an ellipsoid of eccentricity `e`
/// (Karney 2011 `taupf`). The exact closed form (no series).
fn taupf(tau: f64, e: f64) -> f64 {
    let tau1 = tau.hypot(1.0);
    let sig = (e * (e * tau / tau1).atanh()).sinh();
    sig.hypot(1.0) * tau - sig * tau1
}

/// The inverse of [`taupf`]: `tan φ` from `tan χ`, by Newton's method (Karney
/// 2011 `tauf`) — a handful of iterations to full `f64` precision.
fn tauf(taup: f64, e: f64) -> f64 {
    let e2m = 1.0 - e * e; // 1 − e²
    let mut tau = taup / e2m; // initial guess
    let tol = f64::EPSILON.sqrt() / 10.0;
    let stol = tol * taup.abs().max(1.0);
    for _ in 0..5 {
        let taupa = taupf(tau, e);
        let dtau = (taup - taupa) * (1.0 + e2m * tau * tau)
            / (e2m * tau.hypot(1.0) * taupa.hypot(1.0));
        tau += dtau;
        if dtau.abs() < stol {
            break;
        }
    }
    tau
}

/// A **conformal** (angle-preserving) projection — a property defined through the
/// trait: at every point the scale factor is the *same in all directions*
/// (isotropic), so local shapes and angles are preserved (at the cost of area).
/// [`point_scale`](Self::point_scale) is that single factor. Mercator, the
/// conformal conics, and the stereographic projections are conformal; equal-area
/// projections (Albers) are not.
pub trait ConformalProjection: Projection {
    /// The **isotropic point scale factor** at geodetic latitude `lat_deg`: the
    /// local ratio of projected distance to true distance, equal along the
    /// meridian and the parallel (that equality *is* conformality).
    ///
    /// The meridian = parallel scale equality is **proof-backed at the Jacobian
    /// level** (not just the scale value) by `riemann`: the keystone
    /// `Hausdorff/Rhumb.hasDerivAt_isoLat` (`dψ/dφ = sec φ`) gives
    /// `Hausdorff/Mercator.mercator_conformal` (Web Mercator's two principal scale
    /// factors, from the actual partial derivatives, are equal) and its mirror
    /// `Hausdorff/Stereographic.{hasDerivAt_log_conformalT, stereographic_conformal}`
    /// for the stereographic / conformal-conic family.
    fn point_scale(&self, lat_deg: f64) -> f64;
}

impl ConformalProjection for WebMercator {
    /// `sec φ = 1/cos φ` — **proof-backed** by `riemann`
    /// `Hausdorff/Mercator.mercator_point_scale` (`cosh ψ = 1/cos φ`, identical
    /// along meridians and parallels).
    fn point_scale(&self, lat_deg: f64) -> f64 {
        1.0 / lat_deg.to_radians().cos()
    }
}


/// An **equal-area** (area-preserving) projection — the property, defined through
/// the trait, that every region maps to a region of *proportional* area (a
/// constant area scale). It is the conjugate of conformality: an equal-area map
/// preserves areas at the cost of angles. Albers Equal-Area is the instance here;
/// the authalic latitude `q` is exactly what makes its Jacobian area-preserving.
pub trait EqualAreaProjection: Projection {
    /// The constant **area scale factor** (projected area ÷ true area) the
    /// projection holds fixed — `1` for a true-scale equal-area map.
    fn area_scale(&self) -> f64 {
        1.0
    }
}

impl EqualAreaProjection for AlbersEqualArea {}

