//! `s2` spherical measures — robust triangle area and interior angles on the
//! unit sphere, ported from golang/geo's `point_measures.go`.
//!
//! Areas are in **steradians** (multiply by R² for m²; see
//! [`hausdorff_base::earth::EARTH_RADIUS_M`]). Everything is computed from the unit-vector
//! [`S2Point`]s directly, so there is no projection and no pole/antimeridian
//! special-casing.

use hausdorff_base::space::r3::Vector;
use hausdorff_base::space::s1::Angle;

use crate::space::s2::S2Point;

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

/// Signed area of the spherical triangle `abc` in **steradians** by the
/// **van Oosterom–Strackee** formula `2·atan2([a,b,c], 1 + a·b + b·c + c·a)`,
/// `[a,b,c] = a·(b×c)` — a one-shot signed excess straight from the vertices
/// (sign = orientation), with no per-side angles or separate sign test. Equal to
/// [`signed_area`].
///
/// Proof-backed by `riemann` Hausdorff/SphereArea: `excessTan` is the rational
/// `tan(E/2) = [a,b,c]/(1+a·b+b·c+c·a)`, `sphericalExcess = 2·arctan` of it, with
/// `sphericalExcess_pos` (sign = the orientation `triple`) and
/// `sphericalExcess_swap` (a transposition negates it).
#[must_use]
pub fn signed_area_vos(a: S2Point, b: S2Point, c: S2Point) -> f64 {
    let (va, vb, vc) = (v(a), v(b), v(c));
    let num = va.dot(vb.cross(vc));
    let den = 1.0 + va.dot(vb) + vb.dot(vc) + vc.dot(va);
    2.0 * num.atan2(den)
}

/// Signed area (steradians) of the spherical **polygon** through the vertex ring
/// `verts` (≥ 3 points), as the sum of the fan-triangle van Oosterom signed areas
/// from `verts[0]` — the elliptic-geometry analogue of the planar shoelace.
/// Positive for a counter-clockwise ring; `0` for fewer than three vertices.
///
/// Each triangle is [`signed_area_vos`] (proof-backed by `riemann`
/// Hausdorff/SphereArea); the polygon area is their additive fan sum.
#[must_use]
pub fn polygon_signed_area_vos(verts: &[S2Point]) -> f64 {
    if verts.len() < 3 {
        return 0.0;
    }
    let v0 = verts[0];
    verts[1..]
        .windows(2)
        .map(|w| signed_area_vos(v0, w[0], w[1]))
        .sum()
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

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::earth::GeoPoint;

    fn p(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    #[test]
    fn vos_matches_signed_area_and_is_signed() {
        let a = p(0.0, 0.0);
        let b = p(0.0, 20.0);
        let c = p(20.0, 0.0);
        // One-shot van Oosterom signed area matches the l'Huilier·sign area.
        assert!((signed_area_vos(a, b, c) - signed_area(a, b, c)).abs() < 1e-9);
        // Reversing the orientation negates it (sphericalExcess_swap).
        assert!((signed_area_vos(a, c, b) + signed_area_vos(a, b, c)).abs() < 1e-12);
        // A degenerate (repeated-vertex) triangle has zero area.
        assert!(signed_area_vos(a, a, c).abs() < 1e-12);
    }

    #[test]
    fn polygon_area_is_the_fan_sum() {
        // A CCW quadrilateral near the equator.
        let q = [p(0.0, 0.0), p(0.0, 20.0), p(20.0, 20.0), p(20.0, 0.0)];
        let area = polygon_signed_area_vos(&q);
        // Positive (counter-clockwise) and equal to the two fan triangles.
        assert!(area > 0.0);
        let fan = signed_area_vos(q[0], q[1], q[2]) + signed_area_vos(q[0], q[2], q[3]);
        assert!((area - fan).abs() < 1e-12);
        // On a triangle it agrees with the single-triangle area.
        let t = [p(0.0, 0.0), p(0.0, 20.0), p(20.0, 0.0)];
        assert!((polygon_signed_area_vos(&t) - signed_area_vos(t[0], t[1], t[2])).abs() < 1e-12);
        // Degenerate input is zero.
        assert_eq!(polygon_signed_area_vos(&q[..2]), 0.0);
    }
}
