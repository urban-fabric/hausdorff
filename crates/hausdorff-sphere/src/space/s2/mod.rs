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
mod coords;
mod cover;
pub mod edges;
mod loops;
pub mod measures;
mod point;
mod polygon;
mod polyline;
pub mod predicates;
mod rect;
mod region;
mod vector;

pub use cap::Cap;
pub use cellid::S2CellId;
pub use cellunion::CellUnion;
pub use cover::{cover_r2_geometry, cover_r2_geometry_at_level};
pub use edges::{Crossing, crossing_sign, intersection};
pub use loops::Loop;
pub use measures::{angle, girard_area, point_area, signed_area, turn_angle};
pub use point::S2Point;
pub use polygon::Polygon;
pub use polyline::Polyline;
pub use predicates::{robust_sign, sign};
pub use rect::Rect;
pub use region::{Region, RegionCoverer};
