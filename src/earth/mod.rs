//! `earth` — the **Earth-specific** layer: WGS84 lat/lon primitives
//! ([`GeoPoint`]) and great-circle geodesics in metres. The generic,
//! planet-agnostic spherical geometry lives in [`crate::space::s1`]/[`crate::space::s2`] (unit
//! circle/sphere, angles & steradians); [`crate::earth::sphere`] generalizes the metric
//! to any planet radius. This layer fixes the radius to Earth's.
//!
//! Lives next to the H3 grid so `point → cell` shares one type
//! ([`crate::index::h3::H3Cell::for_geo_point`]); downstream crates re-export these so
//! there is one point/bbox/bearing type across the workspace. Unit-typed lengths
//! use the sibling `euclid` units crate (metres etc. live there, not here).

pub mod batch;

mod bbox;
mod bearing;
mod geodesic;
mod karney;
mod point;
mod radial_grid;
mod similarity;

/// Reference ellipsoids — the datum/planet model behind projections and
/// ellipsoidal geodesics.
pub mod ellipsoid;
/// Map projections anchored to the Earth (e.g. UTM).
pub mod projection;
/// The planet-flexible metric layer: arc lengths and areas on a sphere of a
/// chosen radius ([`sphere::Sphere`]), generalizing this layer's fixed radius.
pub mod sphere;

pub use batch::LonLatBatch;
pub use bbox::{BBox, BBoxParseError};
pub use bearing::Bearing;
pub use ellipsoid::{Ellipsoid, GeodesicDirect, GeodesicInverse};
pub use geodesic::{EARTH_RADIUS_M, destination, destination_m, haversine_m, haversine_meters};
pub use point::{FixedGeoPoint, GeoPoint};
pub use radial_grid::{RadialGrid, RingPoint};
pub use similarity::{
    directed_hausdorff_m, frechet_distance, frechet_distance_m, hausdorff_distance,
    hausdorff_distance_m,
};
