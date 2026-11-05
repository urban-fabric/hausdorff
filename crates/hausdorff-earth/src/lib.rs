//! The integrated Earth layer — the pieces of the `earth` namespace that need
//! the spherical, planar, or index machinery, sitting above `hausdorff-base`'s
//! self-contained earth core:
//!
//! - [`earth::arc`] — great-circle arcs over the S2 unit-sphere kernels.
//! - [`earth::projection`] — map projections (planar `r2` output).
//! - `earth::similarity` — great-circle Hausdorff/Fréchet measures (the
//!   generic `algebra` kernels under the `r2` `Haversine` metric).
//! - [`earth::batch`] — columnar lon/lat batches with bulk kernels (distance
//!   fans, whole-batch H3 encoding, bboxes).
//! - [`convert`] — the planar ↔ geographic bridge.
//!
//! [`earth`] re-exports the whole `hausdorff-base` earth core, so a consumer
//! that needs the integrated layers gets one coherent `earth` module rather than
//! having to import the core from one crate and its extensions from another. A
//! consumer that needs *only* the core should depend on `hausdorff-base`
//! directly and skip this crate — and the planar, spherical, and index layers
//! it pulls in — entirely.

pub mod convert;
pub mod earth;

#[cfg(test)]
mod tests;
