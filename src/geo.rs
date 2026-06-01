//! Geographic primitives — the workspace's canonical WGS84 types and spherical
//! geodesics, living next to the H3 grid so `point → cell` shares one type
//! ([`crate::h3::H3Cell::for_geo_point`]).
//!
//! `kolmogorov-geo` and hypatia's `hypatia-osm-primitives` re-export these, so
//! there is one point/bbox/bearing type across the workspace rather than
//! per-crate copies. Unit-typed geodesic sampling helpers use the sibling
//! `euclid` units crate, which stays domain-agnostic and lightweight.

mod bbox;
mod bearing;
mod geodesic;
mod point;
mod radial_grid;
mod similarity;

pub use bbox::{BBox, BBoxParseError};
pub use bearing::Bearing;
pub use geodesic::{EARTH_RADIUS_M, destination, destination_m, haversine_m, haversine_meters};
pub use point::{FixedGeoPoint, GeoPoint};
pub use radial_grid::{RadialGrid, RingPoint};
pub use similarity::{
    directed_hausdorff_m, frechet_distance, frechet_distance_m, hausdorff_distance,
    hausdorff_distance_m,
};
