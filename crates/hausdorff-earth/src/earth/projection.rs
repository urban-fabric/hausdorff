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
use hausdorff_kernels::projection::{
    KrugerSeries, transverse_mercator_forward, transverse_mercator_inverse, web_mercator_forward,
    web_mercator_inverse,
};
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
        let (x, y) = web_mercator_forward(p.lat, p.lon, self.radius);
        Coord::new(x, y)
    }

    fn inverse(&self, c: Coord) -> GeoPoint {
        // φ = arcsin(tanh(y/R)) — the **provably exact** inverse (`riemann`
        // `Hausdorff/Mercator.mercator_inverse`: forward∘inverse = id), since
        // `tanh(ψ(φ)) = sinφ` (`Rhumb.tanh_isoLat`).
        let (lat, lon) = web_mercator_inverse(c.x, c.y, self.radius);
        GeoPoint::from_lat_lon(lat, lon)
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

impl TransverseMercator {
    /// The Karney–Krüger series constants for this projection's ellipsoid.
    fn kruger(&self) -> KrugerSeries {
        KrugerSeries::for_ellipsoid(self.ellipsoid.a, self.ellipsoid.f)
    }
}

impl Projection for TransverseMercator {
    fn project(&self, p: GeoPoint) -> Coord {
        let (e, n) = transverse_mercator_forward(
            p.lat,
            p.lon,
            self.central_meridian,
            self.lat_origin,
            self.k0,
            self.false_easting,
            self.false_northing,
            &self.kruger(),
        );
        Coord::new(e, n)
    }

    fn inverse(&self, coord: Coord) -> GeoPoint {
        let (lat, lon) = transverse_mercator_inverse(
            coord.x,
            coord.y,
            self.central_meridian,
            self.lat_origin,
            self.k0,
            self.false_easting,
            self.false_northing,
            &self.kruger(),
        );
        GeoPoint::from_lat_lon(lat, lon)
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
/// scaled to the unit sphere; the prime-vertical-radius factor. A thin `f64`
/// wrapper over the generic [`hausdorff_kernels::projection::m_factor`]: the
/// conformal/equal-area auxiliary-latitude kernels (`taupf`/`tauf`,
/// `conformal_t`/`authalic_q`, the projection families) all now live in
/// [`hausdorff_kernels::projection`], generic over any
/// [`euclid::numeric::LazyReal`], so they run on `f64`, `f32`, and tensors alike.
fn m_factor(phi: f64, e: f64) -> f64 {
    hausdorff_kernels::projection::m_factor(phi, e)
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
