//! Spherical-triangle solver on the unit sphere: build a triangle from three
//! vertices or three side lengths, recover its interior angles by the spherical
//! law of cosines, and its area by l'Huilier's theorem. Sides and angles are
//! [`Angle`]s (radians); area is in **steradians** — radius-free, like the rest
//! of `s2`.
//!
//! Everything is computed from the central angles between the unit-vector
//! [`S2Point`]s (the same angle-between [`angle` uses), so
//! there is no projection and no pole/antimeridian special-casing.

use hausdorff_base::space::s1::Angle;

use crate::space::s2::S2Point;

/// A spherical triangle, held as its three **side lengths** (the arc lengths
/// opposite each vertex): `a` is opposite vertex `A`, `b` opposite `B`, `c`
/// opposite `C`. Each side is an [`Angle`] in `[0, π]`.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SphericalTriangle {
    /// The side opposite vertex `A` — the arc length `B–C`.
    pub a: Angle,
    /// The side opposite vertex `B` — the arc length `C–A`.
    pub b: Angle,
    /// The side opposite vertex `C` — the arc length `A–B`.
    pub c: Angle,
}

impl SphericalTriangle {
    /// Build the triangle from its three vertices, taking each side as the
    /// **central angle** (geodesic arc length) between the opposite pair of
    /// vertex directions — `a = ∠(B, C)`, `b = ∠(C, A)`, `c = ∠(A, B)` — the
    /// same vector angle [`point_area`
    /// measures the sides with.
    #[must_use]
    pub fn from_vertices(a: S2Point, b: S2Point, c: S2Point) -> Self {
        Self {
            a: b.vector().angle(c.vector()),
            b: c.vector().angle(a.vector()),
            c: a.vector().angle(b.vector()),
        }
    }

    /// Build the triangle directly from its three side lengths.
    #[must_use]
    pub fn from_sides(a: Angle, b: Angle, c: Angle) -> Self {
        Self { a, b, c }
    }

    /// The three sides `(a, b, c)` — the arc lengths opposite vertices
    /// `A`, `B`, `C`.
    #[must_use]
    pub fn sides(&self) -> (Angle, Angle, Angle) {
        (self.a, self.b, self.c)
    }

    /// The three **interior angles** `(A, B, C)` at the vertices, by the
    /// **spherical law of cosines** — `cos A = (cos a − cos b cos c)/(sin b sin c)`,
    /// and cyclically for `B` and `C` (the `acos` argument clamped to `[-1, 1]`).
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphericalLawCosines.sphere_law_cosines`
    /// (the vector identity `A·B = (A·C)(B·C) + (tangent C A)·(tangent C B)`,
    /// whose tangent factors have lengths `sin a`, `sin b` by `tangent_normSq`,
    /// giving `cos c = cos a cos b + sin a sin b cos C`), with the right-angle
    /// case `spherical_pythagoras`.
    #[must_use]
    pub fn interior_angles(&self) -> (Angle, Angle, Angle) {
        let (ca, cb, cc) = (self.a.cos(), self.b.cos(), self.c.cos());
        let (sa, sb, sc) = (self.a.sin(), self.b.sin(), self.c.sin());
        let solve = |num: f64, den: f64| Angle::from_radians((num / den).clamp(-1.0, 1.0).acos());
        (
            solve(ca - cb * cc, sb * sc),
            solve(cb - cc * ca, sc * sa),
            solve(cc - ca * cb, sa * sb),
        )
    }

    /// The **area** of the triangle in **steradians**, by **l'Huilier's
    /// theorem** from the three sides:
    /// `E = 4·atan√(tan(s/2)·tan((s−a)/2)·tan((s−b)/2)·tan((s−c)/2))`, where
    /// `s = (a + b + c)/2` is the semiperimeter. Always non-negative — the same
    /// l'Huilier formula [`point_area`
    /// uses, evaluated from the stored sides.
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphericalLawSines.gram_triple_unit`
    /// (the Gram determinant fixing the side/angle relation) and the l'Huilier
    /// bounds in `Hausdorff/SphereAreaBounds` (`pointArea_nonneg`,
    /// `pointArea_lt_two_pi`): the area is non-negative and below `2π`.
    #[must_use]
    pub fn area(&self) -> f64 {
        let (sa, sb, sc) = (self.a.radians(), self.b.radians(), self.c.radians());
        let s = 0.5 * (sa + sb + sc);
        let tan_prod = (0.5 * s).tan()
            * (0.5 * (s - sa)).tan()
            * (0.5 * (s - sb)).tan()
            * (0.5 * (s - sc)).tan();
        4.0 * tan_prod.max(0.0).sqrt().atan()
    }

    /// The **spherical excess** `E = (A + B + C) − π`, the amount by which the
    /// triangle's interior angles overshoot a flat triangle's `π` — computed from
    /// the [`interior_angles`](Self::interior_angles).
    ///
    /// By **Girard's theorem** the excess *is* the area in **steradians**, so
    /// `spherical_excess() == area()` (to round-off): a spherical triangle's area
    /// is literally how much its angle sum exceeds `π`. For the octant triangle
    /// `E = 3·(π/2) − π = π/2`, matching its `π/2`-sr area.
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphereAreaBounds` (the area is
    /// non-negative and below `2π` — `pointArea_nonneg` / `pointArea_lt_two_pi`,
    /// so `0 ≤ E < 2π`) and `Hausdorff/SphericalLawSines.gram_triple_unit` (the
    /// Gram determinant fixing the side/angle relation behind the angle sum).
    #[must_use]
    pub fn spherical_excess(&self) -> f64 {
        let (a, b, c) = self.interior_angles();
        a.radians() + b.radians() + c.radians() - std::f64::consts::PI
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::space::r3::Vector;

    fn pt(x: f64, y: f64, z: f64) -> S2Point {
        S2Point::from_direction(Vector::new(x, y, z))
    }

    #[test]
    fn octant_triangle_is_a_trirectangular_eighth_of_the_sphere() {
        // The positive-octant triangle: vertices on the three axes.
        let a = pt(1.0, 0.0, 0.0);
        let b = pt(0.0, 1.0, 0.0);
        let c = pt(0.0, 0.0, 1.0);
        let tri = SphericalTriangle::from_vertices(a, b, c);

        // All three sides are quarter-circles (π/2).
        let (sa, sb, sc) = tri.sides();
        let half_pi = std::f64::consts::FRAC_PI_2;
        for s in [sa, sb, sc] {
            assert!((s.radians() - half_pi).abs() < 1e-12, "{}", s.radians());
        }

        // All three interior angles are right angles (π/2) — trirectangular.
        let (ia, ib, ic) = tri.interior_angles();
        for ang in [ia, ib, ic] {
            assert!((ang.radians() - half_pi).abs() < 1e-12, "{}", ang.radians());
        }

        // Area = π/2 sr — one eighth of the whole sphere (4π / 8 = π/2).
        assert!((tri.area() - half_pi).abs() < 1e-12, "{}", tri.area());
        assert!(
            (tri.area() - 4.0 * std::f64::consts::PI / 8.0).abs() < 1e-12,
            "{}",
            tri.area()
        );

        // Girard's theorem: the spherical excess E = 3·(π/2) − π = π/2 IS the
        // area (steradians).
        assert!(
            (tri.spherical_excess() - half_pi).abs() < 1e-12,
            "{}",
            tri.spherical_excess()
        );
        assert!((tri.spherical_excess() - tri.area()).abs() < 1e-12);
    }

    #[test]
    fn spherical_excess_equals_area_on_a_scalene_triangle() {
        // Girard's theorem holds for an irregular triangle too: excess == area.
        let a = pt(1.0, 0.0, 0.0);
        let b = pt(0.0, 1.0, 0.0);
        let c = pt(0.3, 0.2, 1.0);
        let tri = SphericalTriangle::from_vertices(a, b, c);
        assert!(
            (tri.spherical_excess() - tri.area()).abs() < 1e-9,
            "{} vs {}",
            tri.spherical_excess(),
            tri.area()
        );
    }

    #[test]
    fn from_sides_round_trips_and_matches_measures_area() {
        // A smaller scalene triangle built directly from sides equals one built
        // from the matching vertices.
        let a = pt(1.0, 0.0, 0.0);
        let b = pt(0.0, 1.0, 0.0);
        let c = pt(0.2, 0.3, 1.0);
        let from_v = SphericalTriangle::from_vertices(a, b, c);
        let (sa, sb, sc) = from_v.sides();
        let from_s = SphericalTriangle::from_sides(sa, sb, sc);
        assert_eq!(from_v.sides(), from_s.sides());

        // The l'Huilier area from the sides matches measures::point_area on the
        // vertices.
        let measured = crate::space::s2::point_area(a, b, c);
        assert!(
            (from_v.area() - measured).abs() < 1e-12,
            "{} vs {measured}",
            from_v.area()
        );

        // Spherical excess: A + B + C − π equals the area (steradians).
        let (ia, ib, ic) = from_v.interior_angles();
        let excess = ia.radians() + ib.radians() + ic.radians() - std::f64::consts::PI;
        assert!(
            (excess - from_v.area()).abs() < 1e-9,
            "{excess} vs {}",
            from_v.area()
        );
    }
}
