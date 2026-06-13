//! Cross-cutting algebraic abstractions shared across the geometry spaces,
//! organized after `spire`'s `algebra` package.
//!
//! - [`metric`] — [`MetricSpace`]s and the set/curve distances (Hausdorff,
//!   Fréchet) defined generically over them.
//! - [`vector`] — the [`VectorSpace`] / [`InnerProduct`] / [`Norm`] hierarchy
//!   that the coordinate types implement.
//! - [`orientation`] — the signed [`Orientation`] predicate (planar signed area
//!   / ℝ³ triple product) and generic orientation-based containment.
//!
//! The concrete spaces in [`crate::space::r2`]/[`hausdorff_base::space::r3`]/`hausdorff-sphere`'s `space::s2` and the
//! Earth layer build on these so generic algorithms are written once. The trait
//! laws are machine-checked in the `riemann` Lean development (`Geometry.*`).

pub mod metric;
pub mod orientation;
pub mod vector;

pub use metric::{
    FrechetWorkspace, MetricSpace, directed_hausdorff, frechet_distance, hausdorff_distance,
};
pub use orientation::{Orientation, in_triangle};
pub use vector::{InnerProduct, Norm, VectorSpace};
