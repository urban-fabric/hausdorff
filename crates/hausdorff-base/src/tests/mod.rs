//! Test suites for the geospatial foundations, organized by the surface under
//! test.
//!
//! These live here, inside the crate that owns the code, which matters most for
//! this crate: it is the *bottom* of the stack, so a defect in a bearing or an
//! interval should be reproducible without compiling the planar, spherical, or
//! Earth layers that sit above it.
//!
//! | suite | question it answers |
//! |---|---|
//! | [`earth`] | do the WGS84 core types and the spherical geodesics behave? |
//! | [`r1`] | does `r1::Interval` match golang/geo's `r1.Interval` semantics? |
//! | [`geodesy_extremes`] | antipodes, poles, degenerate pairs, and invalid input |
//!
//! Two neighbouring suites deliberately live elsewhere, because they need
//! crates this one sits below:
//!
//! - the **ellipsoidal (Karney) geodesic** goldens and the `r2` metric-trait
//!   conformance are in `hausdorff-planar`, which owns the `Geodesic`/`Haversine`
//!   /`Rhumb` metric implementations those tests drive;
//! - the **great-circle similarity** measures (Hausdorff/Fréchet) are in
//!   `hausdorff-earth`, which owns them.
//!
//! Tests that need private series or crossover constants sit beside their
//! implementation instead — `earth/ellipsoid/clairaut_tests.rs`,
//! `earth/karney/tests.rs`, `earth/latitude/tests.rs`, `space/r3/tests.rs`, and
//! `space/s1/tests.rs`. Anything reachable through the public surface belongs
//! here, where it keeps passing when the algorithm behind it is replaced.

mod earth;
mod geodesy_extremes;
mod r1;
