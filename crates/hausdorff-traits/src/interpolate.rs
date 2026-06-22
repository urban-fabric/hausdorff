//! Interpolation along a segment — after the `riemann` Lean `Interp` abstraction
//! (`Geometry/Interpolate`).
//!
//! Every space has a `lerp p q t` taking `0 ↦ p` and `1 ↦ q` (the endpoint
//! laws). The planar / 3-D straight-line interpolation (riemann `planarInterp` /
//! `chordInterp`; impl'd for `Coord` in `hausdorff-planar` and `Vector` in
//! `hausdorff-base`) and the spherical great-circle SLERP (`Hausdorff/Slerp`,
//! `slerp_zero` / `slerp_angle`; impl'd for `S2Point` in `hausdorff-sphere`) are
//! instances of [`euclid::Interpolate`], so path sampling, midpoints, and
//! subdivision are written once and run in every space.

pub use euclid::{Interpolate, sample, sample_at};
