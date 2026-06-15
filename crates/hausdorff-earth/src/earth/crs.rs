//! `crs` — coordinate reference systems, EPSG-keyed, and transforms between them.
//!
//! A [`Crs`] is a coordinate reference system identified (where standard) by an
//! EPSG code. Three kinds are modelled, matching the EPSG `coord_ref_sys_kind`:
//! **geographic 2D** (lat/lon on a datum, e.g. EPSG:4326), **geocentric** (ECEF
//! Cartesian, EPSG:4978), and **projected** (a [`Projection`] of a geographic
//! CRS, e.g. EPSG:3857 Web Mercator or the UTM zones).
//!
//! [`transform`] converts a coordinate from one CRS to another by pivoting
//! through **WGS84 geocentric (ECEF)** — the one frame every CRS can reach:
//! native → datum geodetic → datum ECEF → (Helmert) → WGS84 ECEF, then the
//! reverse into the target. Routing through ECEF (rather than just lat/lon) is
//! what lets a datum shift — a [`Helmert`] transform — sit on the pivot, so
//! cross-datum conversions are exact to the datum parameters.
//!
//! Datum/ellipsoid conversions (geodetic↔ECEF) are the proof-backed exact
//! closed forms of `riemann` `Hausdorff/{Geodetic,Ellipsoid}`; the projection
//! formulas are Snyder series (numerically tested, not Lean-proved — series, like
//! the Karney geodesic, are outside the closed-form proof boundary); the Helmert
//! datum transform is a closed-form ECEF affine map (Lean-backable, cf.
//! `Hausdorff/AffineMap`).

use crate::earth::GeoPoint;
use crate::earth::projection::{
    AlbersEqualArea, LambertConformalConic, LocalTangentPlane, PolarStereographic, Projection,
    TransverseMercator, WebMercator,
};
use hausdorff_base::earth::{EcefPoint, Ellipsoid};
use hausdorff_planar::space::r2::Coord;

mod datum;
pub use datum::{Datum, DatumShift, Helmert};

/// A concrete map projection — the projected-CRS member of the enum, so a [`Crs`]
/// stays `Copy` and match-able (no `dyn Projection`). Phase 2 adds the conic /
/// stereographic / equal-area members.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProjectionKind {
    /// Spherical Web Mercator (EPSG:3857).
    WebMercator(WebMercator),
    /// Transverse Mercator (UTM and national grids).
    TransverseMercator(TransverseMercator),
    /// Lambert Conformal Conic, 2 standard parallels (national/regional grids).
    LambertConformalConic(LambertConformalConic),
    /// Polar Stereographic (the polar regions).
    PolarStereographic(PolarStereographic),
    /// Albers Equal-Area Conic (statistical/choropleth maps).
    AlbersEqualArea(AlbersEqualArea),
    /// Local tangent-plane (topocentric ENU) — a site-local East/North system.
    LocalTangentPlane(LocalTangentPlane),
}

impl ProjectionKind {
    /// Geographic (lat/lon) → projected planar metres.
    #[must_use]
    pub fn project(self, p: GeoPoint) -> Coord {
        match self {
            ProjectionKind::WebMercator(m) => m.project(p),
            ProjectionKind::TransverseMercator(m) => m.project(p),
            ProjectionKind::LambertConformalConic(m) => m.project(p),
            ProjectionKind::PolarStereographic(m) => m.project(p),
            ProjectionKind::AlbersEqualArea(m) => m.project(p),
            ProjectionKind::LocalTangentPlane(m) => m.project(p),
        }
    }

    /// Projected planar metres → geographic (lat/lon).
    #[must_use]
    pub fn inverse(self, c: Coord) -> GeoPoint {
        match self {
            ProjectionKind::WebMercator(m) => m.inverse(c),
            ProjectionKind::TransverseMercator(m) => m.inverse(c),
            ProjectionKind::LambertConformalConic(m) => m.inverse(c),
            ProjectionKind::PolarStereographic(m) => m.inverse(c),
            ProjectionKind::AlbersEqualArea(m) => m.inverse(c),
            ProjectionKind::LocalTangentPlane(m) => m.inverse(c),
        }
    }
}

/// A coordinate value tagged by which kind of CRS it belongs to — geographic
/// (lat/lon [`GeoPoint`]), projected planar metres ([`Coord`]), or geocentric
/// ECEF ([`EcefPoint`]). [`transform`] consumes and produces these.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CrsCoord {
    /// Lat/lon on a geographic CRS.
    Geographic(GeoPoint),
    /// Easting/northing on a projected CRS.
    Projected(Coord),
    /// X/Y/Z on a geocentric CRS.
    Geocentric(EcefPoint),
}

impl CrsCoord {
    fn geographic(self) -> GeoPoint {
        match self {
            CrsCoord::Geographic(p) => p,
            _ => panic!("expected a Geographic coordinate for a geographic CRS"),
        }
    }
    fn projected(self) -> Coord {
        match self {
            CrsCoord::Projected(c) => c,
            _ => panic!("expected a Projected coordinate for a projected CRS"),
        }
    }
    fn geocentric(self) -> EcefPoint {
        match self {
            CrsCoord::Geocentric(p) => p,
            _ => panic!("expected a Geocentric coordinate for a geocentric CRS"),
        }
    }
}

/// A **coordinate reference system**, keyed (where standard) by EPSG code.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Crs {
    /// Geographic 2D — latitude/longitude on a datum (EPSG:4326 is WGS84).
    Geographic { datum: Datum, epsg: Option<u32> },
    /// Geocentric — ECEF Cartesian on a datum (EPSG:4978 is WGS84).
    Geocentric { datum: Datum, epsg: Option<u32> },
    /// Projected — a map projection of a geographic CRS (EPSG:3857, UTM, …).
    Projected {
        datum: Datum,
        projection: ProjectionKind,
        epsg: Option<u32>,
    },
}

impl Crs {
    /// WGS84 geographic (EPSG:4326).
    #[must_use]
    pub fn wgs84() -> Crs {
        Crs::Geographic {
            datum: Datum::WGS84,
            epsg: Some(4326),
        }
    }

    /// A **local tangent-plane** projected CRS anchored at `origin` (WGS84) —
    /// local East/North metres for site-level work. The origin is user-chosen, so
    /// it is not an EPSG-standard CRS (no code). `transform(Crs::wgs84(), this, p)`
    /// gives `p`'s local East/North.
    #[must_use]
    pub fn local_tangent_plane(origin: GeoPoint) -> Crs {
        Crs::Projected {
            datum: Datum::WGS84,
            projection: ProjectionKind::LocalTangentPlane(LocalTangentPlane::new(
                Ellipsoid::WGS84,
                origin,
            )),
            epsg: None,
        }
    }

    /// The CRS for a standard EPSG `code`, or `None` if unrecognised.
    ///
    /// Recognised: `4326` (WGS84 geographic), `4978` (WGS84 geocentric/ECEF),
    /// `3857`/`900913` (Web Mercator), the WGS84 UTM zones `32601..=32660`
    /// (north) / `32701..=32760` (south); `5070` (NAD83 Conus Albers), `2154`
    /// (RGF93 Lambert-93), `3995` (Arctic) / `3031` (Antarctic) Polar
    /// Stereographic; the Helmert-shifted geographic datums `4277` (OSGB36),
    /// `4230` (ED50), `4267` (NAD27); the British National Grid `27700`; and the
    /// **NTv2 grid-shifted** geographic datums `4314` (DHDN), `4275` (NTF),
    /// `4272` (NZGD49). (NAD83/RGF93 are treated as a zero shift from WGS84.)
    #[must_use]
    pub fn from_epsg(code: u32) -> Option<Crs> {
        match code {
            4326 => Some(Crs::wgs84()),
            4978 => Some(Crs::Geocentric {
                datum: Datum::WGS84,
                epsg: Some(4978),
            }),
            3857 | 900_913 => Some(Crs::Projected {
                datum: Datum::WGS84,
                projection: ProjectionKind::WebMercator(WebMercator::EPSG_3857),
                epsg: Some(3857),
            }),
            32_601..=32_660 => Some(Self::utm((code - 32_600) as u8, true, code)),
            32_701..=32_760 => Some(Self::utm((code - 32_700) as u8, false, code)),
            5070 => Some(Crs::Projected {
                datum: Datum::GRS80,
                projection: ProjectionKind::AlbersEqualArea(AlbersEqualArea {
                    ellipsoid: Ellipsoid::GRS80,
                    lat1: 29.5,
                    lat2: 45.5,
                    lat0: 23.0,
                    lon0: -96.0,
                    false_easting: 0.0,
                    false_northing: 0.0,
                }),
                epsg: Some(5070),
            }),
            2154 => Some(Crs::Projected {
                datum: Datum::GRS80,
                projection: ProjectionKind::LambertConformalConic(LambertConformalConic {
                    ellipsoid: Ellipsoid::GRS80,
                    lat1: 44.0,
                    lat2: 49.0,
                    lat0: 46.5,
                    lon0: 3.0,
                    false_easting: 700_000.0,
                    false_northing: 6_600_000.0,
                }),
                epsg: Some(2154),
            }),
            3995 => Some(Self::polar_stereographic(71.0, true, 3995)),
            3031 => Some(Self::polar_stereographic(-71.0, false, 3031)),
            // Datum-shifted geographic CRSes (non-WGS84 datums).
            4277 => Some(Crs::Geographic {
                datum: Datum::osgb36(),
                epsg: Some(4277),
            }),
            4230 => Some(Crs::Geographic {
                datum: Datum::ed50(),
                epsg: Some(4230),
            }),
            4267 => Some(Crs::Geographic {
                datum: Datum::nad27(),
                epsg: Some(4267),
            }),
            // NTv2 grid-shifted geographic datums (~cm accuracy, embedded grids).
            4314 => Some(Crs::Geographic {
                datum: Datum::dhdn(), // German DHDN — BETA2007
                epsg: Some(4314),
            }),
            4275 => Some(Crs::Geographic {
                datum: Datum::ntf(), // French NTF — ntf_r93
                epsg: Some(4275),
            }),
            4272 => Some(Crs::Geographic {
                datum: Datum::nzgd49(), // New Zealand NZGD49 — nzgd2kgrid0005
                epsg: Some(4272),
            }),
            // OSGB36 / British National Grid (Transverse Mercator on Airy 1830).
            27_700 => Some(Crs::Projected {
                datum: Datum::osgb36(),
                projection: ProjectionKind::TransverseMercator(TransverseMercator {
                    ellipsoid: Datum::osgb36().ellipsoid, // Airy 1830
                    central_meridian: -2.0,
                    lat_origin: 49.0,
                    k0: 0.999_601_271_7,
                    false_easting: 400_000.0,
                    false_northing: -100_000.0,
                }),
                epsg: Some(27_700),
            }),
            _ => None,
        }
    }

    fn utm(zone: u8, north: bool, epsg: u32) -> Crs {
        Crs::Projected {
            datum: Datum::WGS84,
            projection: ProjectionKind::TransverseMercator(TransverseMercator::utm(zone, north)),
            epsg: Some(epsg),
        }
    }

    fn polar_stereographic(lat_ts: f64, north: bool, epsg: u32) -> Crs {
        Crs::Projected {
            datum: Datum::WGS84,
            projection: ProjectionKind::PolarStereographic(PolarStereographic {
                ellipsoid: Ellipsoid::WGS84,
                lat_ts,
                lon0: 0.0,
                north,
                false_easting: 0.0,
                false_northing: 0.0,
            }),
            epsg: Some(epsg),
        }
    }

    /// The EPSG code, if this CRS carries one.
    #[must_use]
    pub fn epsg(self) -> Option<u32> {
        match self {
            Crs::Geographic { epsg, .. }
            | Crs::Geocentric { epsg, .. }
            | Crs::Projected { epsg, .. } => epsg,
        }
    }

    /// The datum this CRS is referenced to.
    #[must_use]
    pub fn datum(self) -> Datum {
        match self {
            Crs::Geographic { datum, .. }
            | Crs::Geocentric { datum, .. }
            | Crs::Projected { datum, .. } => datum,
        }
    }

    /// This CRS's native coordinate → **WGS84 geocentric (ECEF)**, the universal
    /// transform pivot.
    fn to_wgs84_ecef(self, c: CrsCoord) -> EcefPoint {
        match self {
            Crs::Geographic { datum, .. } => {
                datum.to_wgs84_ecef(geodetic_to_ecef(datum, c.geographic()))
            }
            Crs::Geocentric { datum, .. } => datum.to_wgs84_ecef(c.geocentric()),
            Crs::Projected {
                datum, projection, ..
            } => {
                let geo = projection.inverse(c.projected());
                datum.to_wgs84_ecef(geodetic_to_ecef(datum, geo))
            }
        }
    }

    /// **WGS84 geocentric (ECEF)** → this CRS's native coordinate (inverse pivot).
    fn decode_wgs84_ecef(self, e: EcefPoint) -> CrsCoord {
        match self {
            Crs::Geographic { datum, .. } => {
                CrsCoord::Geographic(ecef_to_geodetic(datum, datum.decode_wgs84_ecef(e)))
            }
            Crs::Geocentric { datum, .. } => {
                CrsCoord::Geocentric(datum.decode_wgs84_ecef(e))
            }
            Crs::Projected {
                datum, projection, ..
            } => {
                let geo = ecef_to_geodetic(datum, datum.decode_wgs84_ecef(e));
                CrsCoord::Projected(projection.project(geo))
            }
        }
    }
}

/// Transform a coordinate from one CRS to another, pivoting through WGS84 ECEF.
/// `coord` must be `from`'s native [`CrsCoord`] variant (geographic for a
/// geographic CRS, etc.); the result is `to`'s native variant.
#[must_use]
pub fn transform(from: Crs, to: Crs, coord: CrsCoord) -> CrsCoord {
    to.decode_wgs84_ecef(from.to_wgs84_ecef(coord))
}

/// Geodetic → ECEF on a datum's ellipsoid, at ellipsoidal height 0 (geographic
/// 2D carries no height).
fn geodetic_to_ecef(datum: Datum, p: GeoPoint) -> EcefPoint {
    datum.ellipsoid.geodetic_to_ecef(p, 0.0)
}

/// ECEF → geodetic on a datum's ellipsoid (height discarded — 2D).
fn ecef_to_geodetic(datum: Datum, e: EcefPoint) -> GeoPoint {
    datum.ellipsoid.ecef_to_geodetic(e).0
}

#[cfg(test)]
mod tests;

