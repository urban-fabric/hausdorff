//! Planar [`Interpolate`] instance — straight-line `lerp` on [`Coord`] (riemann
//! `planarInterp`). The trait, the endpoint laws, and the generic `sample` /
//! `midpoint` live upstream in [`hausdorff_geometry_traits`].

use hausdorff_geometry_traits::Interpolate;

use crate::space::r2::Coord;

impl Interpolate for Coord {
    type Scalar = f64;

    /// Planar straight-line `lerp` (riemann `planarInterp`).
    #[inline]
    fn lerp(self, q: Self, t: f64) -> Self {
        Coord::new(self.x + t * (q.x - self.x), self.y + t * (q.y - self.y))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_geometry_traits::sample;

    #[test]
    fn coord_lerp_endpoints_and_midpoint() {
        let p = Coord::new(0.0, 0.0);
        let q = Coord::new(2.0, 4.0);
        assert_eq!(p.lerp(q, 0.0), p); // lerp_zero
        assert_eq!(p.lerp(q, 1.0), q); // lerp_one
        assert_eq!(p.midpoint(q), Coord::new(1.0, 2.0));
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
