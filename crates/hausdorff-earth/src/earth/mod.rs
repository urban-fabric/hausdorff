//! `earth` — the complete Earth-specific layer: the WGS84 core from
//! `hausdorff-base` (re-exported below: [`GeoPoint`], [`BBox`], haversine /
//! rhumb / Karney geodesy, [`ecef`], [`ellipsoid`], [`sphere`]) plus the
//! layers that integrate it with the spherical, planar, and index machinery
//! ([`arc`], [`projection`], `similarity`, [`batch`]).

pub mod batch;

mod similarity;

// Each of these carries its own module header; a second description here would
// be a copy to keep in sync, and rustdoc resolves the merged doc block in this
// module's scope rather than the submodule's, which breaks their intra-doc
// links.
pub mod arc;
pub mod crs;
pub mod gridshift;
pub mod polygon;
pub mod projection;
pub mod sphere;

pub use hausdorff_base::earth::*;

pub use batch::LonLatBatch;
pub use polygon::{GeoMultiPolygon, GeoPolygon};
pub use similarity::{
    directed_hausdorff_m, frechet_distance, frechet_distance_m, hausdorff_distance,
    hausdorff_distance_m,
};
