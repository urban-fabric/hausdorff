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
}
