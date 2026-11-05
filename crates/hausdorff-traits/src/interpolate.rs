//! Interpolation along a segment — after the `riemann` Lean `Interp` abstraction
//! (`Geometry/Interpolate`).
//!
//! Every space has a `lerp p q t` taking `0 ↦ p` and `1 ↦ q` (the endpoint
//! laws). The planar / 3-D straight-line interpolation (riemann `planarInterp` /
//! `chordInterp`; impl'd for `Coord` in `hausdorff-planar` and `Vector` in
//! `hausdorff-base`) and the spherical great-circle SLERP (`Hausdorff/Slerp`,
//! `slerp_zero` / `slerp_angle`; impl'd for `S2Point` in `hausdorff-sphere`) are
//! instances of [`euclid_numeric::Interpolate`], so path sampling, midpoints, and
//! subdivision are written once and run in every space.

use euclid_numeric::{Natural, Real};

pub use euclid_numeric::Interpolate;

/// Interpolate at a natural-number sample index.
///
/// `segments` is the denominator of the interpolation parameter and should be nonzero.
#[must_use]
pub fn sample_at<P, N>(start: P, end: P, index: N, segments: N) -> P
where
    P: Interpolate,
    N: Natural<Real = P::Scalar>,
{
    start.lerp(end, index.to_real() / segments.to_real())
}

/// Return `segments + 1` samples including both endpoints.
///
/// When `segments == 0`, this returns only `start`.
#[must_use]
pub fn sample<P>(start: P, end: P, segments: u64) -> Vec<P>
where
    P: Interpolate,
{
    if segments == 0 {
        return vec![start];
    }

    let mut values = Vec::with_capacity(segments as usize + 1);
    for i in 0..=segments {
        let t = P::Scalar::scalar(i as f64 / segments as f64);
        values.push(start.clone().lerp(end.clone(), t));
    }
    values
}
