//! Internal geometry library for UrbanFabric.
//!
//! A façade over the `hausdorff-*` workspace crates, preserving the original
//! single-crate module surface. Organized into layers:
//! - [`algebra`] — cross-cutting abstractions (metric spaces, vector spaces),
//!   after `spire`'s `algebra` package (in `hausdorff-planar`).
//! - [`space`] — the unit-less geometry, after Google's S2:
//!   [`space::r2`]/[`space::r3`] (Euclidean) and [`space::s1`]/[`space::s2`]
//!   (unit-circle/unit-sphere). Angles live here, not in `euclid`.
//!   (`r2` in `hausdorff-planar`, `s2` in `hausdorff-sphere`, `s1`/`r3` in
//!   `hausdorff-base`.)
//! - [`earth`] — the Earth-specific layer (WGS84 lat/lon, great-circle metres);
//!   [`earth::sphere`] generalizes the metric to any planet radius and
//!   [`earth::projection`] holds map projections. (Core in `hausdorff-base`,
//!   the integrated layers in `hausdorff-earth`.)
//! - [`index`] — discrete spatial indexes ([`index::h3`], [`index::geohash`];
//!   the S2 cell index lives with its geometry in [`space::s2`]) — in
//!   `hausdorff-index`.
//! - [`format`](mod@crate::format) — WKB/WKT/GeoJSON codecs, in
//!   `hausdorff-format`.
//! - [`convert`] bridges the planar and spherical sides.

pub use hausdorff_earth::{convert, earth};
pub use hausdorff_format::format;
pub use hausdorff_index::index;
pub use hausdorff_planar::{algebra, config};

pub mod space {
    //! `space` — the unit-less geometry spaces, organized after Google's S2
    //! into dimensional pieces: [`r2`]/[`r3`] for Euclidean geometry and
    //! [`s1`]/[`s2`] for unit-circle/unit-sphere spherical geometry. All carry
    //! no physical unit (metres enter only in [`crate::earth`]); angles live
    //! in [`s1`], not `euclid`.

    pub use hausdorff_base::space::{r3, s1};
    pub use hausdorff_planar::space::r2;
    pub use hausdorff_sphere::space::s2;
}

/// Transitional alias: `earth` was previously named `geo`. Downstream crates can
/// keep using `hausdorff::geo::…` while they migrate to `hausdorff::earth::…`.
#[doc(hidden)]
pub use self::earth as geo;

/// The profiling substrate, re-exported when the `profiling` feature is on so
/// callers can drive collection: set `ProfileMode`, open a `PerfReport` root,
/// and read the resulting span tree. The kernel span hooks live in the
/// instrumented subcrates (the feature forwards to them).
#[cfg(feature = "profiling")]
pub use huygens;

// The `coord!` / `point!` / `line_string!` / `polygon!` construction macros,
// `#[macro_export]`ed by `hausdorff-planar`, re-exported here so call sites
// keep using `hausdorff::coord! { … }`.
pub use hausdorff_planar::{coord, line_string, point, polygon};

pub use algebra::{InnerProduct, MetricSpace, Norm, VectorSpace};
pub use earth::{FixedGeoPoint, GeoPoint};
pub use space::r2::{
    Area, Bearing, BooleanOps, BoundingRect, Centroid, ChamberlainDuquetteArea, ConcaveHull,
    Contains, ConvexHull, Coord, CoordsIter, Densify, Destination, Dimensions, Direction, Distance,
    Euclidean, Geometry, HasDimensions, Haversine, InterpolatePoint, Intersects, Length, Line,
    LineString, MultiLineString, MultiPoint, MultiPolygon, Operation, Orient, Perimeter, Point,
    Polygon, Rect, Rhumb, Simplify, Triangle, Winding, WindingOrder,
};
