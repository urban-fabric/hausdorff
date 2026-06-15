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

use std::f64::consts::PI;

use crate::earth::GeoPoint;
use hausdorff_base::earth::{Enu, EnuFrame};
use hausdorff_planar::space::r2::{Coord, MapCoords};

pub use crate::earth::ellipsoid::Ellipsoid;

mod families;
pub use families::{AlbersEqualArea, LambertConformalConic, LocalTangentPlane, PolarStereographic};

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
/// also expresses national grids. Snyder's series, accurate to ~1 mm within a
/// few degrees of the central meridian.
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
        let a = self.ellipsoid.a;
        let e2 = self.ellipsoid.e2();
        let ep2 = e2 / (1.0 - e2);

        let phi = p.lat.to_radians();
        let (sphi, cphi, tphi) = (phi.sin(), phi.cos(), phi.tan());
        let n = a / (1.0 - e2 * sphi * sphi).sqrt();
        let t = tphi * tphi;
        let c = ep2 * cphi * cphi;
        let aa = (p.lon - self.central_meridian).to_radians() * cphi;
        let m = meridian_arc(phi, a, e2) - meridian_arc(self.lat_origin.to_radians(), a, e2);

        let easting = self.false_easting
            + self.k0
                * n
                * (aa
                    + (1.0 - t + c) * aa.powi(3) / 6.0
                    + (5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * ep2) * aa.powi(5) / 120.0);
        let northing = self.false_northing
            + self.k0
                * (m + n
                    * tphi
                    * (aa * aa / 2.0
                        + (5.0 - t + 9.0 * c + 4.0 * c * c) * aa.powi(4) / 24.0
                        + (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * ep2) * aa.powi(6)
                            / 720.0));
        Coord::new(easting, northing)
    }

    fn inverse(&self, coord: Coord) -> GeoPoint {
        let a = self.ellipsoid.a;
        let e2 = self.ellipsoid.e2();
        let ep2 = e2 / (1.0 - e2);

        let x = coord.x - self.false_easting;
        let y = coord.y - self.false_northing;

        let m = meridian_arc(self.lat_origin.to_radians(), a, e2) + y / self.k0;
        let mu = m / (a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2.powi(3) / 256.0));
        let e1 = (1.0 - (1.0 - e2).sqrt()) / (1.0 + (1.0 - e2).sqrt());
        let phi1 = mu
            + (3.0 * e1 / 2.0 - 27.0 * e1.powi(3) / 32.0) * (2.0 * mu).sin()
            + (21.0 * e1 * e1 / 16.0 - 55.0 * e1.powi(4) / 32.0) * (4.0 * mu).sin()
            + (151.0 * e1.powi(3) / 96.0) * (6.0 * mu).sin()
            + (1097.0 * e1.powi(4) / 512.0) * (8.0 * mu).sin();

        let (sphi1, cphi1, tphi1) = (phi1.sin(), phi1.cos(), phi1.tan());
        let t1 = tphi1 * tphi1;
        let c1 = ep2 * cphi1 * cphi1;
        let n1 = a / (1.0 - e2 * sphi1 * sphi1).sqrt();
        let r1 = a * (1.0 - e2) / (1.0 - e2 * sphi1 * sphi1).powf(1.5);
        let d = x / (n1 * self.k0);

        let lat = phi1
            - (n1 * tphi1 / r1)
                * (d * d / 2.0
                    - (5.0 + 3.0 * t1 + 10.0 * c1 - 4.0 * c1 * c1 - 9.0 * ep2) * d.powi(4) / 24.0
                    + (61.0 + 90.0 * t1 + 298.0 * c1 + 45.0 * t1 * t1
                        - 252.0 * ep2
                        - 3.0 * c1 * c1)
                        * d.powi(6)
                        / 720.0);
        let lon = self.central_meridian.to_radians()
            + (d - (1.0 + 2.0 * t1 + c1) * d.powi(3) / 6.0
                + (5.0 - 2.0 * c1 + 28.0 * t1 - 3.0 * c1 * c1 + 8.0 * ep2 + 24.0 * t1 * t1)
                    * d.powi(5)
                    / 120.0)
                / cphi1;
        GeoPoint::from_lat_lon(lat.to_degrees(), lon.to_degrees())
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

/// The conformal-latitude function `t = tan(π/4 − φ/2) / [(1 − e sinφ)/(1 + e
/// sinφ)]^{e/2}` (Snyder 15-9) — the isometric-latitude building block of the
/// conformal conic / stereographic projections.
fn conformal_t(phi: f64, e: f64) -> f64 {
    let es = e * phi.sin();
    (PI / 4.0 - phi / 2.0).tan() / ((1.0 - es) / (1.0 + es)).powf(e / 2.0)
}

/// `m = cosφ / √(1 − e² sin²φ)` (Snyder 14-15) — the radius of the parallel
/// scaled to the unit sphere; the prime-vertical-radius factor.
fn m_factor(phi: f64, e: f64) -> f64 {
    let s = phi.sin();
    phi.cos() / (1.0 - e * e * s * s).sqrt()
}

/// Invert the conformal latitude `t` back to geodetic latitude by Snyder's
/// fixed-point iteration (7-9); converges in a few rounds for any terrestrial
/// latitude. Its `e → 0` (spherical) limit is the **exact** closed form
/// `φ = π/2 − 2·arctan(t)`, proof-backed by `riemann`
/// `Hausdorff/Stereographic.conformalLat_inverse`.
fn lat_from_conformal_t(t: f64, e: f64) -> f64 {
    let mut phi = std::f64::consts::FRAC_PI_2 - 2.0 * t.atan();
    for _ in 0..15 {
        let es = e * phi.sin();
        let next =
            std::f64::consts::FRAC_PI_2 - 2.0 * (t * ((1.0 - es) / (1.0 + es)).powf(e / 2.0)).atan();
        if (next - phi).abs() < 1e-14 {
            return next;
        }
        phi = next;
    }
    phi
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

/// Meridian arc length from the equator to latitude `phi` on the ellipsoid.
fn meridian_arc(phi: f64, a: f64, e2: f64) -> f64 {
    a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2.powi(3) / 256.0) * phi
        - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2.powi(3) / 1024.0) * (2.0 * phi).sin()
        + (15.0 * e2 * e2 / 256.0 + 45.0 * e2.powi(3) / 1024.0) * (4.0 * phi).sin()
        - (35.0 * e2.powi(3) / 3072.0) * (6.0 * phi).sin())
}

