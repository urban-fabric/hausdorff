//! Ellipsoidal geodesics — Karney's algorithms for the direct and inverse
//! geodesic problems on an ellipsoid of revolution, accurate to round-off
//! (~15 nm on Earth) for any pair of points, including near-antipodal ones
//! where Vincenty's method fails to converge.
//!
//! Ported from `geographiclib-rs` 0.2.7 (MIT), itself a port of Charles F. F.
//! Karney's GeographicLib (see the crate `NOTICE`); the algorithms are from
//! Karney, *Algorithms for geodesics*, J. Geodesy 87, 43–55 (2013). The port
//! drops the capability masks, arc-mode positioning, and the reduced-length /
//! geodesic-scale outputs — [`Karney::inverse`] and [`Karney::direct`] always
//! solve for distance and azimuths, and [`Karney::ring_perimeter_area`] is the
//! `PolygonArea` accumulator (the C4 area series). `sincosd` uses
//! GeographicLib's `fmod`-based quadrant reduction (the Python/JS form) rather
//! than `remquo`, which hausdorff has no dependency for.
//!
//! The public surface is on [`Ellipsoid`](crate::earth::Ellipsoid)
//! (`distance` / `inverse` / `destination` / `direct`); this module is the
//! kernel. Differential-tested against `geographiclib-rs` (see
//! `tests/geodesic.rs`).
//!
//! Internally the kernel is split into the scalar/angular [`math`] primitives,
//! the Fourier [`series`] coefficient machinery, the [`solver`] (the [`Karney`]
//! per-ellipsoid state and inverse/direct drivers), and the [`line`](mod@line)
//! ([`LineKernel`]) direct-problem state.

mod line;
mod math;
mod series;
mod solver;

#[cfg(test)]
mod tests;

pub(crate) use line::LineKernel;
pub(crate) use solver::{Karney, wgs84};
