//! Planar (Cartesian) geometry — the in-house replacement for `geo`/`geo-types`.
//!
//! `f64`-scalar geometry types ([`Coord`], [`Point`], [`Line`], [`LineString`],
//! [`Polygon`], [`MultiPolygon`], …) and the algorithms over them, mirroring
//! `geo`'s shapes and trait names so consumers migrate by swapping the import
//! path. This is the planar sibling of the spherical WGS84 types in
//! [`hausdorff_base::earth`]; move between the two with `hausdorff-earth`'s `convert`.
//!
//! Ported from georust `geo` / `geo-types` (dual MIT/Apache-2.0); see the crate
//! `NOTICE`.

#[macro_use]
mod macros;

mod convex_polygon;
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
pub mod buffer;
pub mod centroid;
pub mod concave_hull;
pub mod contains;
pub mod convex_hull;
pub mod coords_iter;
pub mod densify;
pub mod dimensions;
pub mod exact;
pub mod geodesic_area;
pub mod index;
pub mod intersects;
pub mod join;
pub mod line_measures;
pub mod make_valid;
pub mod map_coords;
pub mod min_distance;
pub mod prepared;
pub mod relate;
pub mod set_distance;
pub mod simplify;
pub mod winding;

pub mod affine;
pub mod closest_point;
pub mod extremes;
pub mod interior_point;
pub mod isometry;
pub mod line_locate_interpolate;
pub mod remove_repeated_points;
pub mod validation;

pub use convex_polygon::ConvexPolygon;
pub use coord::{Coord, Point};
pub use geometry::Geometry;
pub use line_string::{LineString, MultiLineString};
pub use polygon::{MultiPoint, MultiPolygon, Polygon};
pub use primitives::{ConfidenceEllipse, Ellipse, Line, Rect, Triangle};

pub use area::Area;
pub use batch::PointBatch;
pub use bool_ops::{BooleanOps, Operation, UnaryUnionOperand, unary_union};
pub use bounding_rect::BoundingRect;
pub use buffer::{BufferStyle, buffer};
pub use centroid::Centroid;
pub use concave_hull::ConcaveHull;
pub use contains::Contains;
pub use convex_hull::ConvexHull;
pub use coords_iter::CoordsIter;
pub use densify::Densify;
pub use dimensions::{Dimensions, HasDimensions};
pub use exact::{
    ExactCoord, ExactPredicateError, ExactResult, ExactSegmentIntersection, ExactTrianglePosition,
    InCirclePosition, OrientationSign,
};
pub use geodesic_area::{ChamberlainDuquetteArea, GeodesicArea};
pub use index::StrTree;
pub use intersects::Intersects;
pub use join::{GeometryIndex, nested_loop_join, spatial_join, spatial_join_candidates};
pub use line_measures::{
    Bearing, Destination, Distance, Euclidean, Geodesic, Haversine, InterpolatePoint, Length,
    Perimeter, Rhumb,
};
pub use make_valid::{MakeValid, MakeValidError};
pub use map_coords::MapCoords;
pub use min_distance::{MinDistanceWorkspace, euclidean_distance, euclidean_within};
pub use prepared::PreparedGeometry;
pub use relate::SpatialPredicate;
pub use set_distance::{
    directed_hausdorff_distance, directed_hausdorff_witness, frechet_distance, hausdorff_distance,
};
pub use simplify::Simplify;
pub use winding::{Direction, Orient, Winding, WindingOrder};

pub use affine::{AffineOps, AffineTransform, Rotate, Scale, Skew, Translate};
pub use closest_point::{Closest, ClosestPoint};
pub use extremes::{Extreme, Extremes, Outcome};
pub use interior_point::InteriorPoint;
pub use isometry::Isometry;
pub use line_locate_interpolate::{LineInterpolatePoint, LineLocatePoint};
pub use remove_repeated_points::RemoveRepeatedPoints;
pub use validation::Validation;
