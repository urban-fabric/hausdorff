//! `r3` — three-dimensional Euclidean geometry: the [`Vector`] type, mirroring
//! golang/geo's `r3` package.
//!
//! `r3::Vector` is the substrate of the spherical world: an `S2Point` (in `hausdorff-sphere`)
//! is a *unit* `r3::Vector`, and spherical measures/predicates are computed from
//! vector dot/cross products rather than lat/lon, which is what makes them free
//! of pole and antimeridian singularities.

use crate::space::s1::Angle;

/// A vector in ℝ³.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// The three axes of ℝ³.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Vector {
    /// The vector `(x, y, z)`.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Dot (scalar) product.
    #[inline]
    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross (vector) product.
    #[inline]
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Vector sum. (Inherent, named after golang/geo's `r3.Vector`; the
    /// `VectorSpace` (in `hausdorff-planar`'s `algebra`) impl routes through it.)
    #[inline]
    #[must_use]
    #[allow(clippy::should_implement_trait)]
    pub fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }

    /// Vector difference.
    #[inline]
    #[must_use]
    #[allow(clippy::should_implement_trait)]
    pub fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }

    /// Scalar multiple.
    #[inline]
    #[must_use]
    #[allow(clippy::should_implement_trait)]
    pub fn mul(self, m: f64) -> Self {
        Self::new(self.x * m, self.y * m, self.z * m)
    }

    /// Squared Euclidean norm (`self · self`).
    #[inline]
    #[must_use]
    pub fn norm2(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm(self) -> f64 {
        self.norm2().sqrt()
    }

    /// A unit vector in the same direction, or the zero vector if `self` is zero
    /// (matching `r3.Vector.Normalize`).
    #[must_use]
    pub fn normalize(self) -> Self {
        let n2 = self.norm2();
        if n2 == 0.0 {
            return Self::new(0.0, 0.0, 0.0);
        }
        self.mul(1.0 / n2.sqrt())
    }

    /// Componentwise absolute value.
    #[must_use]
    pub fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// True iff the vector is (near) unit length.
    #[must_use]
    pub fn is_unit(self) -> bool {
        (self.norm2() - 1.0).abs() <= 1e-14
    }

    /// The angle between `self` and `o`, robust for both near-parallel and
    /// near-antiparallel vectors (`atan2(|a×b|, a·b)`).
    #[must_use]
    pub fn angle(self, o: Self) -> Angle {
        Angle::from_radians(self.cross(o).norm().atan2(self.dot(o)))
    }

    /// The axis of the largest-magnitude component. Ties resolve toward `Z` then
    /// `Y` (strict `>`), exactly as `r3.Vector.LargestComponent`.
    #[must_use]
    pub fn largest_component(self) -> Axis {
        let (tx, ty, tz) = (self.x.abs(), self.y.abs(), self.z.abs());
        if tx > ty {
            if tx > tz { Axis::X } else { Axis::Z }
        } else if ty > tz {
            Axis::Y
        } else {
            Axis::Z
        }
    }

    /// The component of `self` **orthogonal to** `normal` — the projection onto
    /// the plane through the origin perpendicular to `normal` (the tangent plane
    /// to the unit sphere at `normal`): `self − normal·(self·normal)`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/CrossTrack.foot_perp_normal` (the
    /// result is ⟂ `normal`, when `normal` is unit) and
    /// `Hausdorff/SphericalLawCosines.tangent_normSq` (its squared length is
    /// `‖self‖² − (self·normal)²`). It is the `tangent` of the Lean development:
    /// the foot of `self`'s perpendicular onto the great circle with pole
    /// `normal`, and the arc-tangent whose length is the sine of the side.
    #[inline]
    #[must_use]
    pub fn project_onto_tangent(self, normal: Self) -> Self {
        self.sub(normal.mul(self.dot(normal)))
    }

    /// The **cosine of the angle** between `self` and `other` —
    /// `(self·other) / (‖self‖·‖other‖)`. By the Cauchy–Schwarz inequality the
    /// result always lies in `[-1, 1]`, so callers may pass it straight to
    /// `acos` **without clamping**. If either vector is the zero vector the angle
    /// is undefined; `0.0` is returned.
    ///
    /// Proof-backed by `riemann` `Hausdorff/HaversineBounds.dot3_sq_le`
    /// (Cauchy–Schwarz, `(p·q)² ≤ (p·p)(q·q)`, so the normalized cosine is in
    /// `[-1, 1]`) and `unitVec_dot_le_one` (the unit-vector special case
    /// `|p·q| ≤ 1`).
    #[inline]
    #[must_use]
    pub fn cos_angle(self, other: Self) -> f64 {
        let denom = self.norm() * other.norm();
        if denom == 0.0 {
            return 0.0;
        }
        // Cauchy–Schwarz guarantees the ratio is in [-1, 1]; the explicit clamp
        // only guards against f64 rounding nudging it a hair outside.
        (self.dot(other) / denom).clamp(-1.0, 1.0)
    }

    /// The **Gram determinant** of the three vectors — the determinant of their
    /// `3×3` matrix of pairwise dot products `Gram(a, b, c)ᵢⱼ = vᵢ·vⱼ`. It equals
    /// the squared volume of the parallelepiped they span, and so satisfies the
    /// identity `triple_product(a, b, c)² == gram_det(a, b, c)`, where the triple
    /// product is `a·(b×c)` (the ℝ³ `Orientation` / scalar triple product of
    /// `(a, b, c)`). For unit vectors it reduces to the spherical form
    /// `1 − cos²α_ab − cos²α_bc − cos²α_ca + 2·cos α_ab·cos α_bc·cos α_ca`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphericalLawSines.gram_triple` (the
    /// unconditional identity `(A·(B×C))² = det Gram(A, B, C)`) and its
    /// unit-vector specialization `gram_triple_unit`.
    #[must_use]
    pub fn gram_det(a: Self, b: Self, c: Self) -> f64 {
        let (aa, bb, cc) = (a.dot(a), b.dot(b), c.dot(c));
        let (ab, bc, ca) = (a.dot(b), b.dot(c), c.dot(a));
        // det of the symmetric Gram matrix [[aa, ab, ca], [ab, bb, bc], [ca, bc, cc]].
        aa * (bb * cc - bc * bc) - ab * (ab * cc - bc * ca) + ca * (ab * bc - bb * ca)
    }

    /// A unit vector orthogonal to `self`, chosen deterministically (built from
    /// the smallest-magnitude axis for numerical stability), as `r3.Vector.Ortho`.
    #[must_use]
    pub fn ortho(self) -> Self {
        let (tx, ty, tz) = (self.x.abs(), self.y.abs(), self.z.abs());
        let other = if tx < ty && tx < tz {
            Self::new(1.0, 0.0, 0.0)
        } else if ty < tz {
            Self::new(0.0, 1.0, 0.0)
        } else {
            Self::new(0.0, 0.0, 1.0)
        };
        self.cross(other).normalize()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn project_onto_tangent_is_perpendicular_with_pythagorean_length() {
        let v = Vector::new(2.0, -3.0, 5.0);
        let n = Vector::new(1.0, 1.0, 1.0).normalize();
        let proj = v.project_onto_tangent(n);
        // The projection is orthogonal to the normal (foot_perp_normal).
        assert!(proj.dot(n).abs() < 1e-12, "{}", proj.dot(n));
        // Its squared length is ‖v‖² − (v·n)² (tangent_normSq).
        let expected = v.dot(v) - v.dot(n).powi(2);
        assert!((proj.norm2() - expected).abs() < 1e-12, "{}", proj.norm2());
    }

    #[test]
    fn project_onto_tangent_matches_inline_subtraction() {
        // The two duplicated inline copies in s2/arc.rs equal this primitive.
        let vp = Vector::new(0.3, 0.7, -0.4);
        let n = Vector::new(-0.2, 0.5, 0.84).normalize();
        let inline = vp.sub(n.mul(vp.dot(n)));
        assert_eq!(vp.project_onto_tangent(n), inline);
        // The latitude-extreme vertex ẑ − n·n_z is ẑ.project_onto_tangent(n).
        let z = Vector::new(0.0, 0.0, 1.0);
        let vertex = Vector::new(-n.x * n.z, -n.y * n.z, 1.0 - n.z * n.z);
        let proj = z.project_onto_tangent(n);
        assert!((proj.x - vertex.x).abs() < 1e-15);
        assert!((proj.y - vertex.y).abs() < 1e-15);
        assert!((proj.z - vertex.z).abs() < 1e-15);
    }

    #[test]
    fn cos_angle_known_configurations_and_range() {
        let x = Vector::new(1.0, 0.0, 0.0);
        let y = Vector::new(0.0, 1.0, 0.0);
        // Orthogonal → 0.
        assert!(x.cos_angle(y).abs() < 1e-15);
        // Same direction (any positive scale) → 1.
        assert!((x.cos_angle(x.mul(7.0)) - 1.0).abs() < 1e-15);
        // Opposite → −1.
        assert!((x.cos_angle(x.mul(-3.0)) + 1.0).abs() < 1e-15);
        // 45° between (1,0,0) and (1,1,0) → 1/√2.
        let d = Vector::new(1.0, 1.0, 0.0);
        assert!((x.cos_angle(d) - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-15);
        // Zero vector is graceful.
        assert_eq!(x.cos_angle(Vector::new(0.0, 0.0, 0.0)), 0.0);
        // Always in [-1, 1] (Cauchy–Schwarz), for assorted vectors.
        let vecs = [
            Vector::new(3.0, -1.0, 2.0),
            Vector::new(-4.0, 5.0, 0.5),
            Vector::new(0.1, 0.2, -9.0),
            Vector::new(-2.0, -2.0, -2.0),
        ];
        for &u in &vecs {
            for &v in &vecs {
                let c = u.cos_angle(v);
                assert!((-1.0..=1.0).contains(&c), "{c}");
                // It is the cosine of the angle the robust `angle` method reports.
                assert!((c - u.angle(v).radians().cos()).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn gram_det_equals_triple_product_squared() {
        let triples = [
            (
                Vector::new(1.0, 0.0, 0.0),
                Vector::new(0.0, 1.0, 0.0),
                Vector::new(0.0, 0.0, 1.0),
            ),
            (
                Vector::new(2.0, -1.0, 3.0),
                Vector::new(0.5, 4.0, -2.0),
                Vector::new(-3.0, 1.0, 1.0),
            ),
            (
                Vector::new(1.0, 1.0, 0.0),
                Vector::new(0.0, 1.0, 1.0),
                Vector::new(1.0, 0.0, 1.0),
            ),
            // A coplanar (degenerate) triple → triple product 0 → gram_det 0.
            (
                Vector::new(1.0, 2.0, 3.0),
                Vector::new(2.0, 4.0, 6.0),
                Vector::new(-1.0, 0.0, 1.0),
            ),
        ];
        for (a, b, c) in triples {
            let triple = a.dot(b.cross(c));
            let gram = Vector::gram_det(a, b, c);
            assert!(
                (triple * triple - gram).abs() < 1e-9,
                "triple²={} vs gram={gram}",
                triple * triple
            );
        }
    }

    #[test]
    fn gram_det_unit_vectors_match_spherical_cosine_form() {
        // For unit vectors: gram_det = 1 − Σcos² + 2∏cos.
        let unit = [
            Vector::new(0.3, -0.5, 0.8).normalize(),
            Vector::new(-0.7, 0.2, 0.4).normalize(),
            Vector::new(0.1, 0.9, -0.2).normalize(),
        ];
        let (a, b, c) = (unit[0], unit[1], unit[2]);
        let (ab, bc, ca) = (a.dot(b), b.dot(c), c.dot(a));
        let spherical = 1.0 - ab * ab - bc * bc - ca * ca + 2.0 * ab * bc * ca;
        assert!((Vector::gram_det(a, b, c) - spherical).abs() < 1e-12);
    }
}
