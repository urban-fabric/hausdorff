//! Planar (Cartesian) geometry — the in-house replacement for `geo`/`geo-types`.
//!
//! `f64`-scalar geometry types ([`Coord`], [`Point`], [`Line`], [`LineString`],
//! [`Polygon`], [`MultiPolygon`], …) and the algorithms over them, mirroring
//! `geo`'s shapes and trait names so consumers migrate by swapping the import
//! path. This is the planar sibling of the spherical WGS84 types in
//! [`crate::earth`]; move between the two with [`crate::convert`].
//!
//! Ported from georust `geo` / `geo-types` (dual MIT/Apache-2.0); see the crate
//! `NOTICE`.

#[macro_use]
mod macros;

mod coord;
mod geometry;
mod geometry_algorithms;
mod line_string;
mod polygon;
mod primitives;
mod util;

pub mod area;
pub mod batch;
pub mod bool_ops;
pub mod bounding_rect;
pub mod centroid;
pub mod concave_hull;
pub mod contains;
pub mod convex_hull;
pub mod coords_iter;
pub mod densify;
pub mod dimensions;
pub mod geodesic_area;
pub mod index;
pub mod intersects;
pub mod join;
pub mod line_measures;
pub mod map_coords;
pub mod min_distance;
pub mod prepared;
pub mod relate;
pub mod set_distance;
pub mod simplify;
pub mod winding;

pub use coord::{Coord, Point};
pub use geometry::Geometry;
pub use line_string::{LineString, MultiLineString};
pub use polygon::{MultiPoint, MultiPolygon, Polygon};
pub use primitives::{Line, Rect, Triangle};

pub use area::Area;
pub use batch::PointBatch;
pub use bool_ops::{BooleanOps, Operation};
pub use bounding_rect::BoundingRect;
pub use centroid::Centroid;
pub use concave_hull::ConcaveHull;
pub use contains::Contains;
pub use convex_hull::ConvexHull;
pub use coords_iter::CoordsIter;
pub use densify::Densify;
pub use dimensions::{Dimensions, HasDimensions};
pub use geodesic_area::{ChamberlainDuquetteArea, GeodesicArea};
pub use index::StrTree;
pub use intersects::Intersects;
pub use join::{GeometryIndex, nested_loop_join, spatial_join, spatial_join_candidates};
pub use line_measures::{
    Bearing, Destination, Distance, Euclidean, Geodesic, Haversine, Length, Perimeter,
};
pub use map_coords::MapCoords;
pub use min_distance::{MinDistanceWorkspace, euclidean_distance, euclidean_within};
pub use prepared::PreparedGeometry;
pub use relate::SpatialPredicate;
pub use set_distance::{directed_hausdorff_distance, frechet_distance, hausdorff_distance};
pub use simplify::Simplify;
pub use winding::{Direction, Orient, Winding, WindingOrder};
