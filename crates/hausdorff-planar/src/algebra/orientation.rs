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

use crate::space::r2::{Coord, Point};
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

impl Orientation for Point {
    /// A [`Point`] is a [`Coord`] with geometry semantics, so its orientation is
    /// the same planar signed area `(b−a) × (c−a)`; delegates to `Coord`.
    #[inline]
    fn orient(a: Self, b: Self, c: Self) -> f64 {
        Coord::orient(a.0, b.0, c.0)
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

/// **Collinearity** of three planar points — `true` iff `a`, `b`, `c` lie on one
/// line, i.e. the orientation determinant `cross(a, b, c)` is zero within `tol`.
/// The exact `tol = 0` test is [`Orientation::degenerate`]; this is the tolerant
/// convenience for inexact (floating-point) coordinates.
///
/// Proof-backed by `riemann` `Hausdorff/LineGeometry`: `collinear_interp`
/// (points interpolated along `a–b` have `cross = 0`) and its converse
/// `interp_locate_of_collinear` — together `cross a b c = 0 ⇔ c` lies on the
/// line `a–b`.
#[must_use]
pub fn is_collinear<O: Orientation>(a: O, b: O, c: O, tol: f64) -> bool {
    O::orient(a, b, c).abs() <= tol
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

/// **Strict point-in-convex-polygon**, generic over any [`Orientation`] space —
/// the n-gon generalization of [`in_triangle`].
///
/// `p` is strictly inside the convex polygon `vertices` (given in
/// counter-clockwise order, `≥ 3` vertices) iff it lies strictly to the **left**
/// of every directed edge `vᵢ → vᵢ₊₁`. The "left of every edge" region is convex
/// — proved once generically in `riemann` (`OrientGeom.leftOfAll_convex`,
/// instantiated planar `leftOfAllStrict_convex` and spherical
/// `sphereLeftOfAll_convex`) — so this is sound for planar polygons and spherical
/// convex polygons alike. Fewer than 3 vertices bound no region (`false`).
#[must_use]
pub fn left_of_all<O: Orientation>(vertices: &[O], p: O) -> bool {
    let n = vertices.len();
    n >= 3 && (0..n).all(|i| O::orient(vertices[i], vertices[(i + 1) % n], p) > 0.0)
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
        assert!(!Coord::turns_positive(a, c, b)); // CW
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
    fn point_orient_matches_coord() {
        use crate::space::r2::Point;
        // `Point` flows through the generic predicates identically to `Coord`.
        let a = Point::new(0.0, 0.0);
        let b = Point::new(4.0, 0.0);
        let c = Point::new(0.0, 4.0);
        assert!(Point::turns_positive(a, b, c)); // CCW
        assert!(Point::degenerate(a, b, Point::new(2.0, 0.0))); // collinear
        assert!(in_triangle(a, b, c, Point::new(1.0, 1.0)));
        assert!(!in_triangle(a, b, c, Point::new(3.0, 3.0)));
        // Same verdict as the underlying `Coord`.
        assert_eq!(Point::orient(a, b, c), Coord::orient(a.0, b.0, c.0));
    }

    #[test]
    fn vector_orient_triple_product() {
        let x = Vector::new(1.0, 0.0, 0.0);
        let y = Vector::new(0.0, 1.0, 0.0);
        let z = Vector::new(0.0, 0.0, 1.0);
        assert!(Vector::orient(x, y, z) > 0.0); // right-handed
        assert!(Vector::degenerate(x, y, Vector::new(1.0, 1.0, 0.0))); // coplanar
    }

    #[test]
    fn is_collinear_detects_the_line() {
        let a = Coord::new(0.0, 0.0);
        let b = Coord::new(2.0, 2.0);
        // c on the line a–b → collinear.
        assert!(is_collinear(a, b, Coord::new(5.0, 5.0), 1e-9));
        // c off the line → not collinear.
        assert!(!is_collinear(a, b, Coord::new(5.0, 4.0), 1e-9));
        // An interpolated point (collinear_interp) is collinear.
        assert!(is_collinear(a, b, Coord::new(0.5, 0.5), 1e-9));
        // The tolerance admits a tiny perpendicular offset.
        assert!(is_collinear(a, b, Coord::new(1.0, 1.0 + 1e-12), 1e-9));
        // Works through the `Point` orientation too.
        assert!(is_collinear(
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(3.0, 0.0),
            1e-9
        ));
    }

    #[test]
    fn point_in_convex_polygon_works() {
        // A CCW unit square.
        let square = [
            Coord::new(0.0, 0.0),
            Coord::new(2.0, 0.0),
            Coord::new(2.0, 2.0),
            Coord::new(0.0, 2.0),
        ];
        assert!(left_of_all(&square, Coord::new(1.0, 1.0))); // centre, inside
        assert!(!left_of_all(&square, Coord::new(3.0, 1.0))); // outside
        assert!(!left_of_all(&square, Coord::new(1.0, 0.0))); // on an edge → not strict
        // Agrees with in_triangle on the first three (CCW) vertices.
        assert!(left_of_all(&square[..3], Coord::new(1.5, 0.25)));
        // Degenerate: fewer than three vertices bound no region.
        assert!(!left_of_all(&square[..2], Coord::new(1.0, 1.0)));
    }
}
