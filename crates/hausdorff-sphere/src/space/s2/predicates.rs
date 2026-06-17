//! `s2` orientation predicates on the unit sphere, ported from golang/geo's
//! `predicates.go`. Computed from vector triple products — no projection.

use super::S2Point;

/// True iff `a, b, c` wind counter-clockwise, i.e. `c` lies to the left of the
/// directed great-circle edge `a → b`.
///
/// Computed as the signed volume of the parallelepiped, `(c × a) · b > 0` — the
/// order golang/geo's `Sign` uses so that `sign(a,b,c)` and `sign(c,b,a)` are
/// never both true.
#[must_use]
pub fn sign(a: S2Point, b: S2Point, c: S2Point) -> bool {
    c.vector().cross(a.vector()).dot(b.vector()) > 0.0
}

/// Orientation of `a, b, c`: `+1` counter-clockwise, `−1` clockwise, `0`
/// degenerate (two points equal, or exactly collinear in floating point).
///
/// This is the fast triage form; unlike golang/geo's full `RobustSign` it does
/// not fall back to exact arithmetic for points that are collinear to within
/// floating-point error (it returns `0` instead). Rotation-invariant and
/// sign-reversing under argument swaps for all non-degenerate inputs.
#[must_use]
pub fn robust_sign(a: S2Point, b: S2Point, c: S2Point) -> i32 {
    let det = c.vector().cross(a.vector()).dot(b.vector());
    if det > 0.0 {
        1
    } else if det < 0.0 {
        -1
    } else {
        0
    }
}

/// `S2Point` implements the shared [`Orientation`](hausdorff_planar::algebra::Orientation)
/// predicate as the scalar triple product `(a × b) · c` of the position
/// vectors — the same right-handed convention as the ℝ³ `Vector` instance, so
/// the generic orientation algorithms (`in_triangle`, …) apply to spherical
/// triangles. (Its laws are the `riemann` `triple_*` theorems.)
impl hausdorff_planar::algebra::Orientation for S2Point {
    type Scalar = f64;

    #[inline]
    fn orient(a: Self, b: Self, c: Self) -> f64 {
        a.vector().cross(b.vector()).dot(c.vector())
    }
}

#[cfg(test)]
mod orientation_trait_tests {
    use super::S2Point;
    use hausdorff_base::earth::GeoPoint;
    use hausdorff_planar::algebra::{Orientation, in_triangle};

    fn pt(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    #[test]
    fn s2point_orientation_and_containment() {
        // Three mid-latitude points spread around, with a near-pole point inside.
        let a = pt(10.0, 0.0);
        let b = pt(10.0, 120.0);
        let c = pt(10.0, 240.0);
        let inside = pt(80.0, 0.0);
        assert!(S2Point::orient(a, b, c).abs() > 0.0);
        assert!(in_triangle(a, b, c, inside));
    }
}
