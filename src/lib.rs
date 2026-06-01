// `deny` rather than `forbid` so the `bytemuck` POD derives in `earth`/`s1` can
// locally `#![allow(unsafe_code)]`; the H3 implementation stays unsafe-free.
#![deny(unsafe_code)]

//! Internal geometry library for UrbanFabric.
//!
//! Organized after Google's S2 into dimensional pieces: [`r2`]/[`r3`] for
//! Euclidean geometry, [`s1`]/[`s2`] for unit-circle/unit-sphere spherical
//! geometry — all unit-less (angles live here, not in `euclid`). [`earth`] is
//! the Earth-specific layer (WGS84 lat/lon + great-circle metres); [`sphere`]
//! generalizes the metric to any planet radius; [`convert`] bridges the planar
//! and spherical sides. Cell indexes live in [`h3`] and [`s2`].

pub mod compat;
pub mod convert;
pub mod earth;
pub mod geohash;
pub mod h3;
pub mod projection;
pub mod r2;
pub mod r2_format;
pub mod r3;
pub mod s1;
pub mod s2;
pub mod sphere;

pub use earth::{FixedGeoPoint, GeoPoint};
pub use r2::{
    Area, Bearing, BoundingRect, Centroid, ChamberlainDuquetteArea, ConcaveHull, Contains,
    ConvexHull, Coord, CoordsIter, Destination, Dimensions, Direction, Distance, Euclidean,
    Geometry, HasDimensions, Haversine, Intersects, Length, Line, LineString, MultiLineString,
    MultiPoint, MultiPolygon, Orient, Point, Polygon, Rect, Simplify, Triangle, Winding,
    WindingOrder,
};
