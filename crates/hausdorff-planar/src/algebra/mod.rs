//! Cross-cutting algebraic abstractions shared across the geometry spaces,
//! organized after `spire`'s `algebra` package.
//!
//! - [`metric`] — [`MetricSpace`]s and the set/curve distances (Hausdorff,
//!   Fréchet) defined generically over them.
//! - [`vector`] — the [`VectorSpace`] / [`InnerProduct`] / [`Norm`] hierarchy
//!   that the coordinate types implement.
//! - [`orientation`] — the signed [`Orientation`] predicate (planar signed area
//!   / ℝ³ triple product) and generic orientation-based containment.
//! - [`simplex`] — the [`Simplex`] containment trait (point-in-triangle R2 /
//!   point-in-tetrahedron R3) via a barycentric orientation partition.
//! - [`interpolate`] — the [`Interpolate`] trait (`lerp` with endpoint laws):
//!   straight-line in R2/R3, great-circle SLERP on S2.
//! - [`region`] — the [`PointRegion`] / [`Prunable`] bounding-volume traits
//!   (point containment + index pruning) over every space's `O(1)` region
//!   primitive (interval, sector, rect, ellipse, cap), with generic spatial
//!   filters written once.
//!
//! The concrete spaces in [`crate::space::r2`]/[`hausdorff_base::space::r3`]/`hausdorff-sphere`'s `space::s2` and the
//! Earth layer build on these so generic algorithms are written once. The trait
//! laws are machine-checked in the `riemann` Lean development (`Geometry.*`).

pub mod interpolate;
pub mod metric;
pub mod orientation;
pub mod region;
pub mod simplex;
pub mod vector;

pub use interpolate::{Interpolate, sample};
pub use metric::{
    FrechetWorkspace, MetricSpace, directed_hausdorff, frechet_distance, hausdorff_distance,
};
pub use orientation::{Orientation, in_triangle, is_collinear, left_of_all};
pub use region::{
    PointRegion, Prunable, any_contained, count_contained, disjoint, filter_contained,
};
pub use simplex::{Barycentric, Simplex, Tetrahedron};
pub use vector::{InnerProduct, Norm, VectorSpace};
