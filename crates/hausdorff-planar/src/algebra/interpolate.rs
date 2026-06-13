//! Interpolation along a segment — after the `riemann` Lean `Interp` abstraction
//! (`Geometry/Interpolate`).
//!
//! Every space has a `lerp p q t` taking `0 ↦ p` and `1 ↦ q` (the endpoint
//! laws). The planar / 3-D straight-line interpolation (riemann `planarInterp` /
//! `chordInterp`) and the spherical great-circle SLERP (`Hausdorff/Slerp`,
//! `slerp_zero` / `slerp_angle`; impl'd for `S2Point` in `hausdorff-sphere`) are
//! instances of one [`Interpolate`] trait, so path sampling, midpoints, and
//! subdivision are written once and run in every space.

use crate::space::r2::Coord;
use hausdorff_base::space::r3::Vector;

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

impl Interpolate for Coord {
    /// Planar straight-line `lerp` (riemann `planarInterp`).
    #[inline]
    fn lerp(self, q: Self, t: f64) -> Self {
        Coord::new(self.x + t * (q.x - self.x), self.y + t * (q.y - self.y))
    }
}

impl Interpolate for Vector {
    /// The R3 chord `lerp` (riemann `chordInterp`).
    #[inline]
    fn lerp(self, q: Self, t: f64) -> Self {
        self.add(q.sub(self).mul(t))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coord_lerp_endpoints_and_midpoint() {
        let p = Coord::new(0.0, 0.0);
        let q = Coord::new(2.0, 4.0);
        assert_eq!(p.lerp(q, 0.0), p); // lerp_zero
        assert_eq!(p.lerp(q, 1.0), q); // lerp_one
        assert_eq!(p.midpoint(q), Coord::new(1.0, 2.0));
    }

    #[test]
    fn vector_lerp_endpoints_and_midpoint() {
        let p = Vector::new(1.0, 2.0, 3.0);
        let q = Vector::new(4.0, 6.0, 8.0);
        assert_eq!(p.lerp(q, 0.0), p);
        assert_eq!(p.lerp(q, 1.0), q);
        assert_eq!(p.midpoint(q), Vector::new(2.5, 4.0, 5.5));
    }

    #[test]
    fn sample_spans_inclusive_and_evenly() {
        let p = Coord::new(0.0, 0.0);
        let q = Coord::new(4.0, 0.0);
        let pts = sample(p, q, 4);
        assert_eq!(pts.len(), 5);
        assert_eq!(pts[0], p);
        assert_eq!(pts[4], q);
        assert_eq!(pts[2], Coord::new(2.0, 0.0));
        assert_eq!(sample(p, q, 0), vec![p]);
    }
}
