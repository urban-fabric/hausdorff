//! Cross-cutting algebraic abstractions shared across the geometry spaces,
//! organized after `spire`'s `algebra` package.
//!
//! The abstract trait surface and the generic algorithms over it live upstream
//! in [`hausdorff_geometry_traits`] — [`Orientation`] (+ [`in_triangle`],
//! [`left_of_all`], [`is_collinear`]), the [`VectorSpace`] / [`InnerProduct`] /
//! [`Norm`] hierarchy, [`Simplex`] / [`Barycentric`], [`Interpolate`] (+
//! [`sample`]), and [`PointRegion`] / [`Prunable`] (+ [`filter_contained`],
//! [`count_contained`], [`any_contained`], [`disjoint`]) — so the `r1` / `s1` /
//! `r3` substrates in `hausdorff-base` implement them too. They are re-exported
//! here unchanged, so the canonical `hausdorff_planar::algebra::*` paths every
//! consumer uses are stable.
//!
//! This module holds the parts that belong to `hausdorff-planar`:
//! - the **planar instances** of those traits — [`crate::space::r2::Coord`] /
//!   `Point` (orientation, vector space, interpolation), `Rect` / `Ellipse` /
//!   `ConfidenceEllipse` (regions), `Triangle` (simplex) — plus the R3
//!   [`Tetrahedron`] simplex;
//! - [`metric`] — the [`MetricSpace`]s and the set/curve distances (Hausdorff,
//!   Fréchet) defined generically over them. It stays here because it is tied to
//!   this crate's [`Distance`](crate::space::r2::Distance) trait.
//!
//! The concrete spaces in [`crate::space::r2`] / `hausdorff_base::space::{r1,s1,r3}`
//! / `hausdorff-sphere`'s `space::s2` and the Earth layer build on these so
//! generic algorithms are written once. The trait laws are machine-checked in the
//! `riemann` Lean development (`Geometry.*` / `Hausdorff.*`).

mod interpolate;
pub mod metric;
mod orientation;
mod region;
mod simplex;
mod vector;

pub use hausdorff_geometry_traits::{
    Barycentric, InnerProduct, Interpolate, Norm, Orientation, PointRegion, Prunable, Simplex,
    VectorSpace, any_contained, count_contained, disjoint, filter_contained, in_triangle,
    is_collinear, left_of_all, sample,
};
pub use metric::{
    FrechetWorkspace, HausdorffWitness, MetricSpace, directed_hausdorff,
    directed_hausdorff_witness, frechet_distance, hausdorff_distance,
};
pub use simplex::Tetrahedron;
