#![deny(unsafe_code)]

//! Internal geometry library for UrbanFabric.
//!
//! Organized into layers:
//! - [`algebra`] — cross-cutting abstractions (metric spaces, vector spaces),
//!   after `spire`'s `algebra` package.
//! - [`space`] — the unit-less geometry, after Google's S2:
//!   [`space::r2`]/[`space::r3`] (Euclidean) and [`space::s1`]/[`space::s2`]
//!   (unit-circle/unit-sphere). Angles live here, not in `euclid`.
//! - [`earth`] — the Earth-specific layer (WGS84 lat/lon, great-circle metres);
//!   [`earth::sphere`] generalizes the metric to any planet radius and
//!   [`earth::projection`] holds map projections.
//! - [`index`] — discrete spatial indexes ([`index::h3`], [`index::geohash`];
//!   the S2 cell index lives with its geometry in [`space::s2`]).
//! - [`format`](mod@crate::format) — WKB/WKT/GeoJSON codecs.
//! - [`convert`] bridges the planar and spherical sides.

pub mod algebra;
pub mod config;
pub mod convert;
pub mod earth;
pub mod format;
pub mod index;
pub mod space;

/// Internal hierarchical-profiling hooks (the `profiling` feature); no-ops when
/// the feature is off.
mod profile;

/// Internal scalar-first SIMD dispatch (runtime-selected) for batch kernels.
mod simd;

/// Transitional alias: `earth` was previously named `geo`. Downstream crates can
/// keep using `hausdorff::geo::…` while they migrate to `hausdorff::earth::…`.
#[doc(hidden)]
pub use earth as geo;

/// The profiling substrate, re-exported when the `profiling` feature is on so
/// callers can drive collection: set `ProfileMode`, open a `PerfReport` root,
/// and read the resulting span tree. See the `profile` module for the kernel
/// span hooks.
#[cfg(feature = "profiling")]
pub use huygens;

pub use algebra::{InnerProduct, MetricSpace, Norm, VectorSpace};
pub use earth::{FixedGeoPoint, GeoPoint};
pub use space::r2::{
    Area, Bearing, BooleanOps, BoundingRect, Centroid, ChamberlainDuquetteArea, ConcaveHull,
    Contains, ConvexHull, Coord, CoordsIter, Destination, Dimensions, Direction, Distance,
    Euclidean, Geometry, HasDimensions, Haversine, Intersects, Length, Line, LineString,
    MultiLineString, MultiPoint, MultiPolygon, Operation, Orient, Perimeter, Point, Polygon, Rect,
    Simplify, Triangle, Winding, WindingOrder,
};
