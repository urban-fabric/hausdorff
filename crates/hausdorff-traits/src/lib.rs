//! Abstract geometry traits shared across the hausdorff geometry spaces
//! (R1/R2/R3/S1/S2/ellipsoid), and the generic algorithms written once over
//! them.
//!
//! The trait *definitions* and space-agnostic free functions live here,
//! upstream of every concrete space, so each space implements them for its own
//! point / region types beside that type's definition — including the `r1` /
//! `s1` / `r3` substrates in [`hausdorff-base`], which would otherwise have no
//! way to reach an abstraction defined in the downstream `hausdorff-planar`.
//!
//! - [`orientation`] — the signed [`Orientation`] predicate + generic
//!   containment ([`in_triangle`], [`left_of_all`], [`is_collinear`]).
//! - [`vector`] — the [`VectorSpace`] / [`InnerProduct`] / [`Norm`] hierarchy.
//! - [`simplex`] — the [`Simplex`] / [`Barycentric`] point-in-simplex trait.
//! - [`interpolate`] — the [`Interpolate`] segment trait + generic [`sample`].
//! - [`region`] — the [`PointRegion`] / [`Prunable`] bounding-volume traits +
//!   the generic spatial filters ([`filter_contained`], [`count_contained`],
//!   [`any_contained`], [`disjoint`]).
//! - [`interval`] — the [`ClosedInterval`] algebra (`contains` / `intersection` /
//!   `union` / `expanded` / hull), shared by `r1::Interval` and `s2::LatInterval`.
//!
//! Each trait's laws are machine-checked in the `riemann` Lean development
//! (`Geometry.*` / `Hausdorff.*`); the per-space impls cite the theorems. The
//! metric / Hausdorff / Fréchet layer is *not* here — it is tied to
//! `hausdorff-planar`'s `Distance` trait and lives in that crate's `algebra`.
//!
//! [`hausdorff-base`]: https://docs.rs/hausdorff-base

pub mod interpolate;
pub mod interval;
pub mod orientation;
pub mod region;
pub mod simplex;
pub mod vector;

pub use interpolate::{Interpolate, sample, sample_at};
pub use interval::ClosedInterval;
pub use orientation::{
    Orientation, contains_convex, convex_winding, in_triangle, is_collinear, is_convex_ring,
    left_of_all,
};
pub use region::{
    PointRegion, Prunable, any_contained, count_contained, disjoint, filter_contained,
};
pub use simplex::{Barycentric, Simplex};
pub use vector::{InnerProduct, Norm, VectorSpace};
