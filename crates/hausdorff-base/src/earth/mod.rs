//! `earth` — the **Earth-specific** layer: WGS84 lat/lon primitives
//! ([`GeoPoint`]) and great-circle geodesics in metres. The generic,
//! planet-agnostic spherical geometry lives in [`crate::space::s1`] /
//! `hausdorff-sphere`'s `space::s2` (unit circle/sphere, angles & steradians);
//! `hausdorff-earth`'s `earth::sphere` generalizes the metric to any planet
//! radius. This layer fixes the radius to Earth's.
//!
//! Downstream crates re-export these so there is one point/bbox/bearing type
//! across the workspace (`point → cell` shares this type with the H3 grid in
//! `hausdorff-index`). Unit-typed lengths use the sibling `euclid` units crate
//! (metres etc. live there, not here).
//!
//! The Earth-facing layers that need the spherical/planar/index machinery —
//! `arc`, `projection`, `similarity`, `batch` — live in `hausdorff-earth`,
//! which re-exports this whole module so the assembled `earth` namespace stays
//! one surface.

mod bbox;
mod bearing;
mod geodesic;
mod karney;
mod point;
mod radial_grid;
mod rhumb;

/// Geocentric (ECEF) coordinates and local east/north/up frames.
pub mod ecef;
/// Reference ellipsoids — the datum/planet model behind projections and
/// ellipsoidal geodesics.
pub mod ellipsoid;
/// The ellipsoid latitude system (geodetic / geocentric / reduced) and the
/// surface-of-revolution radii of curvature, proof-backed by `riemann`.
pub mod latitude;

pub use bbox::{BBox, BBoxParseError};
pub use bearing::Bearing;
pub use ecef::{EcefPoint, Enu, EnuFrame};
pub use ellipsoid::{Ellipsoid, GeodesicDirect, GeodesicInverse, GeodesicLine};
pub use latitude::geodetic_normal;
pub use geodesic::{EARTH_RADIUS_M, destination, destination_m, haversine_m, haversine_meters};
// Workspace-internal: the radius-parameterized kernels behind
// `hausdorff-earth`'s `earth::sphere::Sphere` (hence `#[doc(hidden)]`).
#[doc(hidden)]
pub use geodesic::{destination_on_m, haversine_on_m};
pub use point::{FixedGeoPoint, GeoPoint};
pub use radial_grid::{RadialGrid, RingPoint};
pub use rhumb::{
    rhumb_bearing, rhumb_destination, rhumb_destination_m, rhumb_distance_m, rhumb_distance_meters,
    rhumb_point_at_ratio,
};
