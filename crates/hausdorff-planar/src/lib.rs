//! Planar (Cartesian) geometry and the cross-cutting abstractions over it —
//! the in-house replacement for `geo`/`geo-types`.
//!
//! - [`space::r2`] — the planar `f64` geometry types and algorithms (ported
//!   from georust `geo`; see the repo `NOTICE`).
//! - [`algebra`] — metric-/vector-space abstractions, after `spire`'s
//!   `algebra` package. Lives with `r2` because the two are mutually coupled:
//!   the blanket `MetricSpace` impl is over `r2`'s `Distance`, and `r2`'s
//!   set-distances lift the generic `algebra` kernels.
//! - [`config`] — the `hausdorff.toml` machine-tunable kernel parameters.
//!
//! The geographic metrics (`Haversine`, `Rhumb`, `Geodesic`) delegate to the
//! validated kernels in `hausdorff-base`'s `earth` layer.

pub mod algebra;
pub mod config;
pub mod space;

/// Internal hierarchical-profiling hooks (the `profiling` feature); no-ops when
/// the feature is off.
mod profile;

/// Internal scalar-first SIMD dispatch (runtime-selected) for batch kernels.
mod simd;

/// The profiling substrate, re-exported when the `profiling` feature is on so
/// callers can drive collection: set `ProfileMode`, open a `PerfReport` root,
/// and read the resulting span tree. See the `profile` module for the kernel
/// span hooks.
#[cfg(feature = "profiling")]
pub use huygens;
