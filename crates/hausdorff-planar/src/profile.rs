//! Inert profiling hooks.
//!
//! Kernel call sites still bind [`span()`](crate::profile::span) guards so the
//! instrumentation shape remains visible in code, but hausdorff no longer
//! depends on a profiling substrate.

/// Inert profiling span. Bind it to a named local so it lives to the end of the
/// scope (`let _span = profile::span("name")`).
#[inline(always)]
#[must_use]
pub(crate) fn span(_name: &'static str) -> Inert {
    Inert
}

/// Zero-sized inert span guard.
pub(crate) struct Inert;
