//! Optional in-process profiling hooks (the `profiling` feature).
//!
//! [`span()`](crate::profile::span) opens a hierarchical huygens span for the
//! enclosing scope; bind it
//! to a local (`let _span = profile::span("…")`) and it closes on drop, nesting
//! under any ambient [`PerfReport`](huygens::PerfReport) root. With the
//! `profiling` feature **off** it is an inert zero-sized guard the optimizer
//! removes entirely, so the kernels carry no profiling cost in normal builds.
//!
//! To collect, enable `profiling`, set the runtime mode
//! (`hausdorff::huygens::set_mode(ProfileMode::Heavy)`), and open a root
//! [`PerfReport::span`](huygens::PerfReport::span) around the work; the kernel
//! spans then assemble into its tree.

/// Open a profiling span for the enclosing scope. Bind it to a named local so it
/// lives to the end of the scope (`let _span = profile::span("name")`); it
/// closes on drop. No-op (an inert ZST) without the `profiling` feature.
#[cfg(feature = "profiling")]
#[inline]
#[must_use]
pub(crate) fn span(name: &'static str) -> huygens::span::SpanGuard<'static> {
    huygens::span::span(name)
}

/// Inert profiling span — the `profiling` feature is disabled, so this is a
/// zero-sized no-op the optimizer elides.
#[cfg(not(feature = "profiling"))]
#[inline(always)]
#[must_use]
pub(crate) fn span(_name: &'static str) -> Inert {
    Inert
}

/// Zero-sized inert span guard used when the `profiling` feature is disabled.
#[cfg(not(feature = "profiling"))]
pub(crate) struct Inert;
