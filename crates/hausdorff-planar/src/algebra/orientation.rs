//! Planar [`Orientation`] instances — [`Coord`] and the [`Point`] newtype both
//! use the signed area `(b−a) × (c−a)` (twice the triangle area). The trait and
//! the generic orientation algorithms (`in_triangle`, `left_of_all`,
//! `is_collinear`) live upstream in [`hausdorff_traits`]; the laws are
//! machine-checked in `riemann` (`Geometry.OrientGeom`).

use hausdorff_traits::Orientation;

use crate::space::r2::{Coord, Point};

impl Orientation for Coord {
    type Scalar = f64;

    /// Planar signed area `(b−a) × (c−a)` (twice the triangle area).
    #[inline]
    fn orient(a: Self, b: Self, c: Self) -> f64 {
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    }
}

impl Orientation for Point {
    type Scalar = f64;

    /// A [`Point`] is a [`Coord`] with geometry semantics, so its orientation is
    /// the same planar signed area `(b−a) × (c−a)`; delegates to `Coord`.
    #[inline]
    fn orient(a: Self, b: Self, c: Self) -> f64 {
        Coord::orient(a.0, b.0, c.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_traits::{in_triangle, is_collinear, left_of_all};

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
