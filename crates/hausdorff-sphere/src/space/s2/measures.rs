//! `s2` spherical measures — robust triangle area and interior angles on the
//! unit sphere, ported from golang/geo's `point_measures.go`.
//!
//! Areas are in **steradians** (multiply by R² for m²; see
//! [`hausdorff_base::earth::EARTH_RADIUS_M`]). Everything is computed from the unit-vector
//! [`S2Point`]s directly, so there is no projection and no pole/antimeridian
//! special-casing.

use hausdorff_base::space::r3::Vector;
use hausdorff_base::space::s1::Angle;

use super::S2Point;

#[inline]
fn v(p: S2Point) -> Vector {
    p.vector()
}

/// A vector perpendicular to both `a` and `b`, numerically robust when the two
/// points are nearly equal (`(a + b) × (b − a)`, golang/geo's `PointCross`).
fn point_cross(a: S2Point, b: S2Point) -> Vector {
    v(a).add(v(b)).cross(v(b).sub(v(a)))
}

/// Area of the spherical triangle `abc` in **steradians**, via l'Huilier's
/// theorem — accurate for triangles of every size. Always non-negative; points
/// must be unit length and no two antipodal.
#[must_use]
pub fn point_area(a: S2Point, b: S2Point, c: S2Point) -> f64 {
    let sa = v(b).angle(v(c)).radians();
    let sb = v(c).angle(v(a)).radians();
    let sc = v(a).angle(v(b)).radians();
    let s = 0.5 * (sa + sb + sc);

    if s >= 3e-4 {
        let dmin = s - sa.max(sb).max(sc);
        if dmin < 1e-2 * s.powi(5) {
            let area = girard_area(a, b, c);
            if dmin < s * 0.1 * (area + 5e-15) {
                return area;
            }
        }
    }

    let tan_prod =
        (0.5 * s).tan() * (0.5 * (s - sa)).tan() * (0.5 * (s - sb)).tan() * (0.5 * (s - sc)).tan();
    4.0 * tan_prod.max(0.0).sqrt().atan()
}

/// Area of the spherical triangle `abc` via Girard's formula (spherical excess).
/// Faster than [`point_area`] but less accurate for tiny triangles.
#[must_use]
pub fn girard_area(a: S2Point, b: S2Point, c: S2Point) -> f64 {
    let ab = point_cross(a, b);
    let bc = point_cross(b, c);
    let ac = point_cross(a, c);
    let area = ab.angle(ac).radians() - ab.angle(bc).radians() + bc.angle(ac).radians();
    area.max(0.0)
}

/// Signed triangle area: `+point_area` when `abc` is counter-clockwise, `−`
/// otherwise (`0` if the three points are degenerate/collinear).
#[must_use]
pub fn signed_area(a: S2Point, b: S2Point, c: S2Point) -> f64 {
    f64::from(super::predicates::robust_sign(a, b, c)) * point_area(a, b, c)
}

/// Interior angle at vertex `b` of triangle `abc`, in `[0, π]`. Symmetric:
/// `angle(a, b, c) == angle(c, b, a)`.
#[must_use]
pub fn angle(a: S2Point, b: S2Point, c: S2Point) -> Angle {
    point_cross(a, b).angle(point_cross(c, b))
}

/// Exterior ("turn") angle at vertex `b` walking `a → b → c`: positive for a
/// left/counter-clockwise turn, negative for a right/clockwise turn — the
/// geodesic curvature at `b`. `turn_angle(a, b, c) == −turn_angle(c, b, a)`.
/// This is the term summed by the Gauss–Bonnet loop-area formula.
#[must_use]
pub fn turn_angle(a: S2Point, b: S2Point, c: S2Point) -> Angle {
    let magnitude = point_cross(a, b).angle(point_cross(b, c));
    if super::predicates::robust_sign(a, b, c) > 0 {
        magnitude
    } else {
        -magnitude
    }
}
