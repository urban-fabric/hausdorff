// `deny` rather than `forbid` so the `bytemuck` POD derives in `geo` can locally
// `#![allow(unsafe_code)]`; the H3 implementation stays unsafe-free under deny.
#![deny(unsafe_code)]

//! Internal geospatial primitives and index implementations for UrbanFabric.

pub mod convert;
pub mod geo;
pub mod h3;
pub mod planar;
pub mod r3;
pub mod s1;
pub mod s2;
pub mod sphere;

pub use geo::{FixedGeoPoint, GeoPoint};
pub use planar::{
    Area, Bearing, BoundingRect, Centroid, ChamberlainDuquetteArea, Contains, ConvexHull, Coord,
    CoordsIter, Destination, Dimensions, Direction, Distance, Euclidean, Geometry, HasDimensions,
    Haversine, Length, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Orient, Point,
    Polygon, Rect, Simplify, Triangle, Winding, WindingOrder,
};
