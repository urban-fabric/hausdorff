//! Backend-neutral geospatial formula kernels.
//!
//! This crate owns formulas that can be expressed over `euclid`'s lazy numeric
//! traits. It deliberately avoids Hausdorff's eager data records, storage
//! encodings, spatial indexes, and backend crates so the same formulas can be
//! lowered by scalar, tensor/XLA, or dataframe/expression backends.

pub mod earth;
pub mod interval;
pub mod planar;
pub mod projection;
pub mod spherical1;
