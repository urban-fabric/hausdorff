//! In-house S2 cell-indexing implementation — a from-scratch port of the S2
//! Geometry cell hierarchy (Google's S2), giving the workspace a second
//! hierarchical index system alongside `hausdorff-index`'s `index::h3`.
//!
//! Covers the indexing core: the 64-bit [`S2CellId`] codec + validity, the
//! cube-to-sphere quadratic projection (lat/lon ↔ cell), `parent` / `children`
//! / child positions, `range_min`/`range_max` containment, edge/vertex/all
//! neighbors, cell centers, Hilbert-order traversal, and hex tokens. Cells
//! share [`hausdorff_base::earth::GeoPoint`] at the lat/lon boundary, exactly like the H3
//! side, so a `point → cell` mapping uses one coordinate type across both
//! systems.
//!
//! Ported from the Go reference implementation (`github.com/golang/geo/s2`).
//! Beyond cells, this module also provides generic unit-sphere geometry built on
//! the [`S2Point`] unit vector: robust [`measures`] (triangle area, angles),
//! orientation [`predicates`], and the [`Loop`] region (spherical polygon) — all
//! computed from vectors, with no projection.

pub mod arc;
mod cap;
mod cellid;
mod cellunion;
mod convex_loop;
mod coords;
mod cover;
pub mod edges;
pub mod great_circle;
mod loops;
pub mod measures;
pub mod metric;
mod point;
mod polygon;
mod polyline;
pub mod predicates;
mod rect;
mod region;
mod relate;
pub mod simplex;
pub mod triangle;
mod vector;

pub use cap::Cap;
pub use cellid::S2CellId;
pub use cellunion::CellUnion;
pub use convex_loop::ConvexLoop;
pub use cover::{
    bound_r2_geometry, bound_r2_geometry_arc, cover_r2_geometry, cover_r2_geometry_at_level,
};
pub use edges::{Crossing, crossing_sign, intersection};
pub use great_circle::GreatCircle;
pub use loops::Loop;
pub use measures::{
    angle, girard_area, point_area, polygon_signed_area_vos, signed_area, signed_area_vos,
    turn_angle,
};
pub use metric::{AVG_EDGE, LengthMetric, MIN_WIDTH, SphereMetric};
pub use point::S2Point;
pub use polygon::Polygon;
pub use polyline::Polyline;
pub use predicates::{robust_sign, sign};
pub use rect::Rect;
pub use region::{Region, RegionCoverer};
pub use relate::SpatialPredicate;
pub use simplex::SphereTriangle;
pub use triangle::SphericalTriangle;
