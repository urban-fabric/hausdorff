//! Orientation predicate — the signed turn / volume test shared across the
//! geometry spaces, after the `riemann` Lean `OrientGeom` abstraction.
//!
//! [`crate::space::r2::Coord`] uses the planar signed area
//! `(b−a) × (c−a)`; [`hausdorff_base::space::r3::Vector`] (the ℝ³ substrate of
//! `space::s2`) uses the scalar triple product `(a × b) · c`. Both are bundled
//! here as [`Orientation`] so generic orientation-based algorithms (containment,
//! hulls, in/out tests) are written once.
//!
//! The laws the predicate obeys — antisymmetry under a swap, invariance under a
//! cyclic rotation, affinity along a segment — are machine-checked in the
//! `riemann` development (`triDet_swap`/`triDet_cyclic`/`triDet_lerp` for the
//! plane, `triple_*` for ℝ³/S²; unified as `Geometry.OrientGeom`). Any space
//! whose `orient` satisfies them inherits the proved containment/convexity
//! results.

use crate::space::r2::Coord;
use hausdorff_base::space::r3::Vector;

/// A space carrying a signed orientation predicate on triples.
pub trait Orientation: Copy {
    /// Signed orientation of `(a, b, c)`: `> 0` a positive turn
    /// (counter-clockwise in R2, right-handed in R3), `< 0` negative,
    /// `0` degenerate (collinear in R2, coplanar in R3).
    fn orient(a: Self, b: Self, c: Self) -> f64;

    /// `true` iff `a, b, c` make a strictly positive turn.
    #[inline]
    #[must_use]
    fn turns_positive(a: Self, b: Self, c: Self) -> bool {
        Self::orient(a, b, c) > 0.0
    }

    /// `true` iff `a, b, c` are degenerate (collinear / coplanar).
    #[inline]
    #[must_use]
    fn degenerate(a: Self, b: Self, c: Self) -> bool {
        Self::orient(a, b, c) == 0.0
    }
}

impl Orientation for Coord {
    /// Planar signed area `(b−a) × (c−a)` (twice the triangle area).
    #[inline]
    fn orient(a: Self, b: Self, c: Self) -> f64 {
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    }
}

impl Orientation for Vector {
    /// Scalar triple product `(a × b) · c` — signed volume of the
    /// parallelepiped; its sign is the right-handed orientation.
    #[inline]
    fn orient(a: Self, b: Self, c: Self) -> f64 {
        a.cross(b).dot(c)
    }
}

/// **Strict point-in-triangle**, generic over any [`Orientation`] space.
///
/// `p` is strictly inside triangle `a, b, c` iff it lies on the inner side of
/// all three directed edges — the same-sign orientation test. The barycentric
/// soundness of this certificate is proved in `riemann`
/// (`point_in_triangle_edges`; the spherical mirror `sphereInTriangle`).
#[must_use]
pub fn in_triangle<O: Orientation>(a: O, b: O, c: O, p: O) -> bool {
    if O::orient(a, b, c) > 0.0 {
        O::orient(a, b, p) > 0.0 && O::orient(b, c, p) > 0.0 && O::orient(c, a, p) > 0.0
    } else {
        O::orient(a, b, p) < 0.0 && O::orient(b, c, p) < 0.0 && O::orient(c, a, p) < 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coord_orient_sign() {
        let a = Coord::new(0.0, 0.0);
        let b = Coord::new(1.0, 0.0);
        let c = Coord::new(0.0, 1.0);
        assert!(Coord::turns_positive(a, b, c)); // CCW
        assert!(Coord::turns_positive(a, c, b) == false); // CW
        assert!(Coord::degenerate(a, b, Coord::new(2.0, 0.0))); // collinear
    }

    #[test]
    fn point_in_triangle_works() {
        let a = Coord::new(0.0, 0.0);
        let b = Coord::new(4.0, 0.0);
        let c = Coord::new(0.0, 4.0);
        assert!(in_triangle(a, b, c, Coord::new(1.0, 1.0)));
        assert!(!in_triangle(a, b, c, Coord::new(3.0, 3.0)));
        // orientation-independent: reversed (CW) winding gives the same verdict
        assert!(in_triangle(a, c, b, Coord::new(1.0, 1.0)));
    }

    #[test]
    fn vector_orient_triple_product() {
        let x = Vector::new(1.0, 0.0, 0.0);
        let y = Vector::new(0.0, 1.0, 0.0);
        let z = Vector::new(0.0, 0.0, 1.0);
        assert!(Vector::orient(x, y, z) > 0.0); // right-handed
        assert!(Vector::degenerate(x, y, Vector::new(1.0, 1.0, 0.0))); // coplanar
    }
}
