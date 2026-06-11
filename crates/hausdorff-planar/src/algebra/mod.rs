//! Cross-cutting algebraic abstractions shared across the geometry spaces,
//! organized after `spire`'s `algebra` package.
//!
//! - [`metric`] — [`MetricSpace`]s and the set/curve distances (Hausdorff,
//!   Fréchet) defined generically over them.
//! - [`vector`] — the [`VectorSpace`] / [`InnerProduct`] / [`Norm`] hierarchy
//!   that the coordinate types implement.
//!
//! The concrete spaces in [`crate::space::r2`]/[`hausdorff_base::space::r3`]/`hausdorff-sphere`'s `space::s2` and the
//! Earth layer build on these so generic algorithms are written once.

pub mod metric;
pub mod vector;

pub use metric::{
    FrechetWorkspace, MetricSpace, directed_hausdorff, frechet_distance, hausdorff_distance,
};
pub use vector::{InnerProduct, Norm, VectorSpace};
