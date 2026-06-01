//! Geographic primitives — the workspace's canonical WGS84 types and spherical
//! geodesics, living next to the H3 grid so `point → cell` shares one type
//! ([`crate::h3::H3Cell::for_geo_point`]).
//!
//! `kolmogorov-geo` and hypatia's `hypatia-osm-primitives` re-export these, so
//! there is one point/bbox/bearing type across the workspace rather than
//! per-crate copies. The types stay free of `euclid`/unit dependencies on
//! purpose — the unit-typed wrappers (`LengthMeters`, …) live in `kolmogorov-geo`
//! — keeping hausdorff light enough for the minimal PBF reader to depend on.

mod bbox;
mod bearing;
mod geodesic;
mod point;

pub use bbox::{BBox, BBoxParseError};
pub use bearing::Bearing;
pub use geodesic::{destination_m, haversine_m, EARTH_RADIUS_M};
pub use point::{FixedGeoPoint, GeoPoint};
