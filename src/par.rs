//! Internal: the one feature-gated seam between sequential and parallel batch
//! execution, after nalgebra's `par_iter` pattern.
//!
//! Batch kernels express their independent per-element work as a closure over an
//! index and call [`map_indexed`]. With the `parallel` feature enabled it runs on
//! rayon's work-stealing pool; otherwise it is a plain sequential map. The scalar
//! path is the reference, and both produce identical results, so the kernels are
//! written once and the choice is a build-time flag.

/// Map `f` over `0..n`, collecting into a `Vec` — on rayon when the `parallel`
/// feature is on, sequentially otherwise.
#[cfg(feature = "parallel")]
pub(crate) fn map_indexed<T, F>(n: usize, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Send + Sync,
{
    use rayon::prelude::*;
    (0..n).into_par_iter().map(f).collect()
}

/// Map `f` over `0..n`, collecting into a `Vec` (sequential reference path).
#[cfg(not(feature = "parallel"))]
pub(crate) fn map_indexed<T, F>(n: usize, f: F) -> Vec<T>
where
    F: Fn(usize) -> T,
{
    (0..n).map(f).collect()
}
