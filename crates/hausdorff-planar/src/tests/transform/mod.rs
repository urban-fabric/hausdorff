//! The same geometry, moved or resampled — and the invariant that it is still the
//! same geometry.
//!
//! What these suites check is mostly **what did not change**: an affine map must
//! preserve area up to its determinant, densifying must not move the line it adds
//! vertices to, and dropping repeated points must leave the shape identical. That
//! makes the oracle the input itself, which is as independent as an oracle gets.
//!
//! [`simplify_topology`] is the one that can genuinely break its subject, since
//! Douglas–Peucker may introduce a self-intersection; its suite is therefore about
//! validity rather than preservation.

mod affine;
mod remove_repeated_points;
mod simplify_topology;
