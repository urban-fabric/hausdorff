//! `earth` — the complete Earth-specific layer: the WGS84 core from
//! `hausdorff-base` (re-exported below: [`GeoPoint`], [`BBox`], haversine /
//! rhumb / Karney geodesy, [`ecef`], [`ellipsoid`], [`sphere`]) plus the
//! layers that integrate it with the spherical, planar, and index machinery
//! ([`arc`], [`projection`], `similarity`, [`batch`]).

pub mod batch;

mod similarity;

/// Great-circle arc operations: interpolation, cross/along-track distances,
/// arc intersection, spherical point-in-ring, and Clairaut-aware arc bboxes.
pub mod arc;
/// Map projections anchored to the Earth (e.g. UTM).
pub mod projection;
/// Coordinate reference systems (EPSG-keyed) and transforms between them.
pub mod crs;
/// High-accuracy datum transforms via embedded NTv2 (`.gsb`) grid-shift files.
pub mod gridshift;
/// The planet-flexible metric layer: arc lengths and areas on a sphere of a
/// chosen radius ([`sphere::Sphere`]), generalizing this layer's fixed radius.
pub mod sphere;
/// WGS84 polygons measured natively on the curved Earth ([`polygon::GeoPolygon`] /
/// [`polygon::GeoMultiPolygon`]) — exact area / perimeter, no projection.
pub mod polygon;

pub use hausdorff_base::earth::*;

pub use batch::LonLatBatch;
pub use polygon::{GeoMultiPolygon, GeoPolygon};
pub use similarity::{
    directed_hausdorff_m, frechet_distance, frechet_distance_m, hausdorff_distance,
    hausdorff_distance_m,
};
