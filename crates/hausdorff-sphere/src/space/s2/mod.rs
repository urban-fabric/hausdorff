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
//! the [`S2Point`] unit vector: robust [`point_area`] / [`angle`] measures,
//! orientation [`robust_sign`], and the [`Loop`] region (spherical polygon) — all
//! computed from vectors, with no projection.

//! ## Layout
//!
//! Twenty-three modules grouped into four, in dependency order: `primitive` (the
//! unit vector and the robust predicates on it), `cell` (the cube-sphere
//! projection and the cell hierarchy), `shape` (the spherical regions), and `ops`
//! (booleans and DE-9IM, which compose several shapes). `vector` and [`arc`] stay
//! at this level — the first because every group uses it, the second because its
//! path is public API.
//!
//! The groups are private, so **the facade below is the whole public surface** and
//! a module can move between groups without any caller noticing. Each group's
//! header says where its boundary falls.
//!
//! ## Crossing a group boundary
//!
//! Six internal items — and only six — are `pub(crate)` rather than group-private,
//! because another group genuinely needs them: `S2Point::vector` and `S2Point::from_vector`
//! (reached from every group — sixteen and nine modules respectively),
//! `coords::face_uv_to_xyz` and `S2CellId::face_uv_bounds` (both for `shape::rect`,
//! which inverts the projection), and `cover::polygon_region` / `cover::valid_lonlat`
//! (for `ops::relate`). Each is annotated where it is declared, saying which group
//! reaches for it.
//!
//! Unlike `hausdorff-planar`'s `r2`, these are *not* collected into a
//! `pub(crate) use` block here, because three of the six are methods and cannot be
//! re-exported — a list holding half the crossings would imply a completeness it
//! does not have. Everything else stays `pub(super)`, which inside a group means
//! group-private, and at this level still means visible throughout `s2`.
//!
//! Grouping tightened one item rather than loosening it: `Loop::convex_winding`
//! was `pub(crate)`, and both callers turned out to be in `shape`, so it is now
//! `pub(super)`.

pub mod arc;

mod vector;

mod cell;
mod ops;
mod primitive;
mod shape;

pub use cell::cellid::S2CellId;
pub use cell::cellunion::CellUnion;
pub use cell::cover::{
    bound_r2_geometry, bound_r2_geometry_arc, cover_r2_geometry, cover_r2_geometry_at_level,
};
pub use cell::metric::{AVG_EDGE, LengthMetric, MIN_WIDTH, SphereMetric};
pub use ops::bool_ops::{
    Operation, SphericalBooleanOps, boolean_op, symmetric_difference, total_area,
};
pub use ops::relate::SpatialPredicate;
pub use primitive::edges::{Crossing, crossing_sign, intersection};
pub use primitive::great_circle::GreatCircle;
pub use primitive::measures::{
    angle, girard_area, point_area, polygon_signed_area_vos, signed_area, signed_area_vos,
    turn_angle,
};
pub use primitive::point::S2Point;
pub use primitive::predicates::{robust_sign, sign};
pub use shape::cap::Cap;
pub use shape::convex_loop::ConvexLoop;
pub use shape::loops::Loop;
pub use shape::polygon::Polygon;
pub use shape::polyline::Polyline;
pub use shape::rect::Rect;
pub use shape::region::{Region, RegionCoverer};
pub use shape::simplex::SphereTriangle;
pub use shape::triangle::SphericalTriangle;
