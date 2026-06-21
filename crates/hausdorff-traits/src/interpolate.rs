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

use euclid::numeric::{LazyNatural, LazyReal};

/// A space whose points interpolate along a segment by a fraction `t`.
pub trait Interpolate: Clone {
    /// The scalar field the interpolation fraction is drawn from.
    type Scalar: LazyReal;

    /// Interpolate from `self` toward `q` by fraction `t`. Obeys the endpoint
    /// laws `lerp(p, q, 0) = p` and `lerp(p, q, 1) = q` (riemann `Interp`).
    #[must_use]
    fn lerp(self, q: Self, t: Self::Scalar) -> Self;

    /// The midpoint `lerp(self, q, 0.5)` — the straight-line midpoint in R2/R3,
    /// the great-circle midpoint on `S2`.
    #[inline]
    #[must_use]
    fn midpoint(self, q: Self) -> Self {
        self.lerp(q, Self::Scalar::half())
    }
}

/// The `index / segments` point along the segment from `p` to `q`.
///
/// This is the fully lazy sampling primitive: the count/index values are
/// [`LazyNatural`] values and are converted into the interpolation scalar with
/// [`LazyNatural::to_real`] without inspecting them as host integers. Callers must
/// avoid `segments == 0`; the host collection adapter [`sample`] handles that
/// structural case.
#[must_use]
pub fn sample_at<P, N>(p: P, q: P, index: N, segments: N) -> P
where
    P: Interpolate,
    N: LazyNatural<Real = P::Scalar>,
{
    p.lerp(q, index.to_real() / segments.to_real())
}

/// `segments + 1` points evenly spaced along the segment from `p` to `q`
/// (inclusive of both ends) — generic over any [`Interpolate`] space, so it is
/// straight-line in R2/R3 and great-circle on `S2`.
///
/// This materializes a host [`Vec`], so the segment count is structural even
/// though it uses the concrete [`LazyNatural`] carrier `u64`. Use [`sample_at`]
/// when the index/count values themselves must stay lazy.
#[must_use]
pub fn sample<P: Interpolate>(p: P, q: P, segments: u64) -> Vec<P> {
    if segments == 0 {
        return vec![p];
    }
    (0..=segments)
        .map(|i| {
            p.clone()
                .lerp(q.clone(), P::Scalar::scalar(i as f64 / segments as f64))
        })
        .collect()
}
