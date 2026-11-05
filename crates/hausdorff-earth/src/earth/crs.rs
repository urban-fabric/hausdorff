//! `crs` — coordinate reference systems, EPSG-keyed, and transforms between them.
//!
//! A [`Crs`] is a coordinate reference system identified (where standard) by an
//! EPSG code. Three kinds are modelled, matching the EPSG `coord_ref_sys_kind`:
//! **geographic** (lat/lon on a datum — 2-D EPSG:4326, or 3-D EPSG:4979 with
//! ellipsoidal height), **geocentric** (ECEF Cartesian, EPSG:4978), and
//! **projected** (a [`Projection`] of a geographic CRS, e.g. EPSG:3857 Web
//! Mercator or the UTM zones). A [`CompoundCrs`] pairs a horizontal CRS with a
//! [`VerticalCrs`] (height).
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
use hausdorff_base::earth::{BBox, EcefPoint, Ellipsoid};
use hausdorff_planar::space::r2::Coord;

mod datum;
pub use datum::{Datum, DatumShift, Helmert};

mod method;
mod proj_string;
mod projjson;
mod wkt;

mod projection_kind;
pub use projection_kind::ProjectionKind;

mod coordinate;
pub use coordinate::{CrsCoord, GeoDim};

mod registry;

/// A **coordinate reference system**, keyed (where standard) by EPSG code.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Crs {
    /// Geographic — latitude/longitude on a datum, 2-D (EPSG:4326 is WGS84) or
    /// 3-D with ellipsoidal height (`dim = D3`, e.g. EPSG:4979).
    Geographic {
        datum: Datum,
        epsg: Option<u32>,
        dim: GeoDim,
    },
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
    /// WGS84 geographic 2-D (EPSG:4326).
    #[must_use]
    pub fn wgs84() -> Crs {
        Crs::Geographic {
            datum: Datum::WGS84,
            epsg: Some(4326),
            dim: GeoDim::D2,
        }
    }

    /// WGS84 geographic 3-D (EPSG:4979) — lat/lon + ellipsoidal height. Unlike
    /// [`wgs84`](Self::wgs84), a transform into this CRS retains the height.
    #[must_use]
    pub fn wgs84_3d() -> Crs {
        Crs::Geographic {
            datum: Datum::WGS84,
            epsg: Some(4979),
            dim: GeoDim::D3,
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

    /// This CRS's native coordinate → **WGS84 geocentric (ECEF)**, the universal
    /// transform pivot. Geographic coordinates carry their ellipsoidal height
    /// into ECEF; projected ones are planar (height 0).
    fn to_wgs84_ecef(self, c: CrsCoord) -> EcefPoint {
        match self {
            Crs::Geographic { datum, .. } => {
                let (geo, height) = c.as_geographic();
                datum.to_wgs84_ecef(geodetic_to_ecef(datum, geo, height))
            }
            Crs::Geocentric { datum, .. } => datum.to_wgs84_ecef(c.geocentric()),
            Crs::Projected {
                datum, projection, ..
            } => {
                let geo = projection.inverse(c.projected());
                datum.to_wgs84_ecef(geodetic_to_ecef(datum, geo, 0.0))
            }
        }
    }

    /// **WGS84 geocentric (ECEF)** → this CRS's native coordinate (inverse pivot).
    /// A 2-D geographic CRS drops the recovered height; a 3-D one (`dim = D3`)
    /// keeps it.
    fn decode_wgs84_ecef(self, e: EcefPoint) -> CrsCoord {
        match self {
            Crs::Geographic { datum, dim, .. } => {
                let (geo, height) = ecef_to_geodetic(datum, datum.decode_wgs84_ecef(e));
                let height = match dim {
                    GeoDim::D2 => 0.0,
                    GeoDim::D3 => height,
                };
                CrsCoord::Geographic { point: geo, height }
            }
            Crs::Geocentric { datum, .. } => CrsCoord::Geocentric(datum.decode_wgs84_ecef(e)),
            Crs::Projected {
                datum, projection, ..
            } => {
                let (geo, _height) = ecef_to_geodetic(datum, datum.decode_wgs84_ecef(e));
                CrsCoord::Projected(projection.project(geo))
            }
        }
    }
}

/// Transform a coordinate from one CRS to another, pivoting through WGS84 ECEF.
/// `coord` must be `from`'s native [`CrsCoord`] variant (geographic for a
/// geographic CRS, etc.); the result is `to`'s native variant.
///
/// **Transforming to the same CRS returns the coordinate unchanged**, bit for
/// bit. Without the short circuit the pivot still runs — geodetic → ECEF →
/// geodetic — and while the round trip is accurate it is not *exact*: it moved an
/// ellipsoidal height by a couple of nanometres. Harmless on its own, but a
/// no-op that is not a no-op compounds in a pipeline that normalizes everything
/// to a working CRS and transforms repeatedly, and it makes an identity hard to
/// assert.
#[must_use]
pub fn transform(from: Crs, to: Crs, coord: CrsCoord) -> CrsCoord {
    if from == to {
        return coord;
    }
    to.decode_wgs84_ecef(from.to_wgs84_ecef(coord))
}

mod compound;
pub use compound::{CompoundCrs, VerticalCrs};

/// Geodetic → ECEF on a datum's ellipsoid, at ellipsoidal `height` metres.
fn geodetic_to_ecef(datum: Datum, p: GeoPoint, height: f64) -> EcefPoint {
    datum.ellipsoid.geodetic_to_ecef(p, height)
}

/// ECEF → geodetic on a datum's ellipsoid, returning lat/lon and the recovered
/// ellipsoidal height.
fn ecef_to_geodetic(datum: Datum, e: EcefPoint) -> (GeoPoint, f64) {
    datum.ellipsoid.ecef_to_geodetic(e)
}

#[cfg(test)]
mod tests;
