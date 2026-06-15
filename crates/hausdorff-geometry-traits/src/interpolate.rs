//! Interpolation along a segment — after the `riemann` Lean `Interp` abstraction
//! (`Geometry/Interpolate`).
//!
//! Every space has a `lerp p q t` taking `0 ↦ p` and `1 ↦ q` (the endpoint
//! laws). The planar / 3-D straight-line interpolation (riemann `planarInterp` /
//! `chordInterp`; impl'd for `Coord` in `hausdorff-planar` and `Vector` in
//! `hausdorff-base`) and the spherical great-circle SLERP (`Hausdorff/Slerp`,
//! `slerp_zero` / `slerp_angle`; impl'd for `S2Point` in `hausdorff-sphere`) are
//! instances of one [`Interpolate`] trait, so path sampling, midpoints, and
//! subdivision are written once and run in every space.

/// A space whose points interpolate along a segment by a fraction `t`.
pub trait Interpolate: Copy {
    /// Interpolate from `self` toward `q` by fraction `t`. Obeys the endpoint
    /// laws `lerp(p, q, 0) = p` and `lerp(p, q, 1) = q` (riemann `Interp`).
    #[must_use]
    fn lerp(self, q: Self, t: f64) -> Self;

    /// The midpoint `lerp(self, q, 0.5)` — the straight-line midpoint in R2/R3,
    /// the great-circle midpoint on `S2`.
    #[inline]
    #[must_use]
    fn midpoint(self, q: Self) -> Self {
        self.lerp(q, 0.5)
    }
}

/// `segments + 1` points evenly spaced along the segment from `p` to `q`
/// (inclusive of both ends) — generic over any [`Interpolate`] space, so it is
/// straight-line in R2/R3 and great-circle on `S2`. `segments == 0` yields
/// `[p]`.
#[must_use]
pub fn sample<P: Interpolate>(p: P, q: P, segments: usize) -> Vec<P> {
    if segments == 0 {
        return vec![p];
    }
    (0..=segments)
        .map(|i| p.lerp(q, i as f64 / segments as f64))
        .collect()
}
