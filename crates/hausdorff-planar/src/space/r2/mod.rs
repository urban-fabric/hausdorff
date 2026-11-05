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
//!
//! ## One name, one path
//!
//! **Every public item is re-exported here, and every module holding one is
//! private.** `Area` is `r2::Area` and nothing else; there is no
//! `r2::area::Area`. The file layout is therefore an internal matter and can be
//! regrouped without touching a single caller.
//!
//! This is the rule the workspace manifest gives for deleting the old `hausdorff`
//! umbrella crate — that it "gave most types two public paths … so call sites
//! drifted between them and 'where does this live' stopped having one answer" —
//! applied one level down, because this module had re-created it. 38 modules were
//! `pub mod` *and* re-exported, so 107 names had two paths; meanwhile 26 public
//! items had only the long one, which is what made the long path look load-bearing
//! and invited the drift. Both halves are closed by the same rule.
//!
//! A module that genuinely wants to be a namespace — `validation::InvalidPolygon`
//! rather than a bare `InvalidPolygon` — is a reasonable alternative, but it has to
//! be *either* a namespace *or* behind the facade. It cannot be both, and mixing
//! the two per-module is how this drifted in the first place.
//!
//! [`validation`] is the one module that takes that option, for reasons given on
//! it. Everything else is behind the facade.
//!
//! ## Groups
//!
//! The 44 implementation modules are grouped by the question they answer, because
//! one flat list of 47 sorts well and reads badly — it gave no hint where a new
//! module belongs or which existing one already does the job.
//!
//! | group | what its modules do | n |
//! |---|---|---|
//! | `model` | the geometry types and the traits that walk their coordinates | 9 |
//! | `measure` | geometry in, number out | 9 |
//! | `predicate` | geometry in, verdict out | 6 |
//! | `construct` | geometry in, new geometry out | 10 |
//! | `transform` | the same geometry, moved or resampled | 7 |
//! | `query` | finding geometries among many | 3 |
//!
//! Each group is a private module with its own header explaining where its
//! boundary falls and which members look misfiled and why. The remaining three
//! modules stay at this level: `macros` and `util` because they are below the
//! taxonomy rather than in it, and [`validation`] because its path is public API.
//!
//! **The grouping is a layout decision, not an API one.** Because the facade above
//! is the whole public surface, a module can move between groups without any
//! caller noticing — which is what made this reorganisation safe to do, and is the
//! point of doing the one-path rule first. Nothing outside this file knows that
//! `Area` is implemented in `measure/area.rs`.
//!
//! One consequence worth stating: a group never reaches sideways into another
//! group's internals. Five internal helpers, from three modules, genuinely need to
//! cross a boundary, and they are re-exported `pub(crate)` below so there is a
//! single place that records it. If a sixth appears it belongs there too, rather
//! than as a direct path into a sibling group — the list is short enough to be
//! worth keeping short.

#[macro_use]
mod macros;

mod util;

// Six groups, declared in the order the module docs list them. Each is private:
// the facade below is the whole public surface, which is what lets a module move
// between groups without any caller noticing.
mod construct;
mod measure;
mod model;
mod predicate;
mod query;
mod transform;

// The one public namespace, kept at this level because its path *is* API — see
// its own header for why it is not behind the facade. A grouped module's path is
// internal and may move; this one may not.
pub mod validation;

// Internal helpers that cross a group boundary, routed through the facade for the
// same reason the public items are: a group can be reorganised without its
// siblings having to learn where its files went. `pub(crate)`, so none of this
// widens the public surface — these are the only sideways references in the crate,
// and keeping them visible in one place is what stops the next one being added
// silently.
pub(crate) use construct::convex_hull::convex_hull_coords;
pub(crate) use measure::min_distance::{seg_box, segment_gap_sq, segments};
pub(crate) use predicate::intersects::segments_intersect;

pub use model::convex_polygon::ConvexPolygon;
pub use model::coord::{Coord, Point};
pub use model::coords_iter::CoordsIter;
pub use model::dimensions::{Dimensions, HasDimensions};
pub use model::geometry::Geometry;
pub use model::line_string::{LineString, MultiLineString};
pub use model::polygon::{MultiPoint, MultiPolygon, Polygon};
pub use model::primitives::{ConfidenceEllipse, Ellipse, Line, Rect, Triangle};

pub use measure::area::Area;
pub use measure::bounding_rect::BoundingRect;
pub use measure::bounding_shapes::{
    AsGeometry, BoundingShapes, Circle, OrientedRect, hull_area, oriented_envelope_area,
};
pub use measure::centroid::Centroid;
pub use measure::extremes::{Extreme, Extremes, Outcome};
pub use measure::geodesic_area::{ChamberlainDuquetteArea, GeodesicArea};
pub use measure::line_measures::{
    Bearing, Destination, Distance, Euclidean, Geodesic, Haversine, InterpolatePoint, Length,
    Perimeter, Rhumb,
};
pub use measure::min_distance::{
    MinDistanceWorkspace, euclidean_distance, euclidean_within, point_segment_distance,
    segment_segment_distance,
};
pub use measure::set_distance::{
    directed_hausdorff_distance, directed_hausdorff_witness, frechet_distance, hausdorff_distance,
};

pub use predicate::contains::Contains;
// The `exact` predicates carry a different contract from their `f64` siblings —
// exact arithmetic, with overflow reported rather than absorbed — so their names
// stay verbose enough to say so at the call site.
pub use predicate::exact::{
    ExactCoord, ExactPredicateError, ExactResult, ExactSegmentIntersection, ExactTrianglePosition,
    InCirclePosition, OrientationSign, in_circle_determinant, in_circle_position, orient,
    orientation_determinant, point_in_triangle, point_on_segment, segment_intersection,
    segments_properly_cross, triangle_contains_point,
};
pub use predicate::intersects::Intersects;
pub use predicate::prepared::PreparedGeometry;
pub use predicate::relate::{IntersectionMatrix, SpatialPredicate, geometry_dimension, relate};
pub use predicate::winding::{Direction, Orient, Winding, WindingOrder};

pub use construct::bool_ops::{BooleanOps, Operation, UnaryUnionOperand, unary_union};
pub use construct::buffer::{BufferStyle, buffer};
pub use construct::closest_point::{Closest, ClosestPoint};
pub use construct::concave_hull::ConcaveHull;
pub use construct::convex_hull::ConvexHull;
pub use construct::interior_point::InteriorPoint;
pub use construct::line_locate_interpolate::{LineInterpolatePoint, LineLocatePoint};
pub use construct::make_valid::{MakeValid, MakeValidError};
pub use construct::polygonize::{Polygonize, Polygonized};
pub use construct::triangulate::{
    Triangulate, Triangulation, Voronoi, delaunay, delaunay_triangles, nearest_site, voronoi,
    voronoi_polygons,
};

pub use transform::affine::{AffineOps, AffineTransform, Rotate, Scale, Skew, Translate};
pub use transform::densify::Densify;
pub use transform::isometry::Isometry;
pub use transform::map_coords::MapCoords;
pub use transform::remove_repeated_points::RemoveRepeatedPoints;
pub use transform::simplify::Simplify;
pub use transform::simplify_topology::SimplifyTopologyPreserving;

pub use query::batch::PointBatch;
pub use query::index::StrTree;
pub use query::join::{GeometryIndex, nested_loop_join, spatial_join, spatial_join_candidates};
