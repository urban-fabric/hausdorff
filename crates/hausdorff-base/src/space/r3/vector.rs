use super::*;

/// A vector in ℝ³, generic over the lazy-real scalar `T` (defaulting to `f64`).
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Vector<T = f64> {
    pub x: T,
    pub y: T,
    pub z: T,
}

/// The three axes of ℝ³.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl<T> Vector<T> {
    /// The vector `(x, y, z)`.
    #[must_use]
    pub const fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }
}

impl<T> Vector<T>
where
    T: Real,
{
    /// Dot (scalar) product.
    #[inline]
    #[must_use]
    pub fn dot(self, o: Self) -> T {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross (vector) product.
    ///
    /// Proof-backed by `riemann` `Hausdorff/CrossProduct`: it is orthogonal to
    /// both inputs (`cross3_dot_left`/`_right`), anticommutative
    /// (`cross3_anticomm`), and satisfies the Lagrange identity
    /// `‖a×b‖² = ‖a‖²‖b‖² − (a·b)²` (`cross3_normSq`) — the basis of
    /// [`angle`](Self::angle).
    #[inline]
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y.clone() * o.z.clone() - self.z.clone() * o.y.clone(),
            self.z.clone() * o.x.clone() - self.x.clone() * o.z,
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
    pub fn mul(self, m: T) -> Self {
        Self::new(self.x * m.clone(), self.y * m.clone(), self.z * m)
    }

    /// Squared Euclidean norm (`self · self`).
    #[inline]
    #[must_use]
    pub fn norm2(self) -> T {
        self.clone().dot(self)
    }

    /// Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm(self) -> T {
        self.norm2().sqrt()
    }

    /// A unit vector in the same direction, or the zero vector if `self` is zero
    /// (matching `r3.Vector.Normalize`).
    #[must_use]
    pub fn normalize(self) -> Self {
        let n2 = self.clone().norm2();
        let is_zero = n2.eq(&T::zero());
        let safe_n2 = T::select(&is_zero, T::one(), n2);
        let unit = self.mul(safe_n2.sqrt().recip());
        Self::select(&is_zero, Self::new(T::zero(), T::zero(), T::zero()), unit)
    }

    /// Componentwise absolute value.
    #[must_use]
    pub fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// True iff the vector is (near) unit length.
    #[must_use]
    pub fn is_unit(self) -> <T as Compare>::Mask {
        self.is_unit_with_tolerance(T::scalar(1e-14))
    }

    /// True iff the vector is unit length within a caller-supplied squared-norm
    /// tolerance.
    #[must_use]
    pub fn is_unit_with_tolerance(self, tolerance: T) -> <T as Compare>::Mask {
        (self.norm2() - T::one()).abs().le(&tolerance)
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
        let scale = self.clone().dot(normal.clone());
        self.sub(normal.mul(scale))
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
    pub fn cos_angle(self, other: Self) -> T {
        let denom = self.clone().norm() * other.clone().norm();
        let is_zero = denom.eq(&T::zero());
        let safe_denom = T::select(&is_zero, T::one(), denom);
        // Cauchy–Schwarz guarantees the ratio is in [-1, 1]; the explicit clamp
        // only guards against f64 rounding nudging it a hair outside.
        let cosine = (self.dot(other) / safe_denom).clamp(-T::one(), T::one());
        T::select(&is_zero, T::zero(), cosine)
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
    pub fn gram_det(a: Self, b: Self, c: Self) -> T {
        let (aa, bb, cc) = (
            a.clone().dot(a.clone()),
            b.clone().dot(b.clone()),
            c.clone().dot(c.clone()),
        );
        let (ab, bc, ca) = (
            a.clone().dot(b.clone()),
            b.clone().dot(c.clone()),
            c.clone().dot(a.clone()),
        );
        // det of the symmetric Gram matrix [[aa, ab, ca], [ab, bb, bc], [ca, bc, cc]].
        aa * (bb.clone() * cc.clone() - bc.clone() * bc.clone())
            - ab.clone() * (ab.clone() * cc - bc.clone() * ca.clone())
            + ca.clone() * (ab * bc - bb * ca)
    }

    #[inline]
    #[must_use]
    pub(crate) fn select(mask: &<T as Compare>::Mask, on_true: Self, on_false: Self) -> Self {
        Self::new(
            T::select(mask, on_true.x, on_false.x),
            T::select(mask, on_true.y, on_false.y),
            T::select(mask, on_true.z, on_false.z),
        )
    }

    /// A unit vector orthogonal to `self`, chosen deterministically (built from
    /// the smallest-magnitude axis for numerical stability), as `r3.Vector.Ortho`.
    #[must_use]
    pub fn ortho(self) -> Self {
        let tx = self.x.clone().abs();
        let ty = self.y.clone().abs();
        let tz = self.z.clone().abs();
        let choose_x = tx.clone().lt(&ty).and(&tx.lt(&tz));
        let choose_y = choose_x.clone().not().and(&ty.lt(&tz));
        let x_axis = Self::new(T::one(), T::zero(), T::zero());
        let y_axis = Self::new(T::zero(), T::one(), T::zero());
        let z_axis = Self::new(T::zero(), T::zero(), T::one());
        let other = Self::select(&choose_x, x_axis, Self::select(&choose_y, y_axis, z_axis));
        self.cross(other).normalize()
    }
}

impl<T> Vector<T>
where
    T: Real + PartialOrd,
{
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
}

impl Vector<f64> {
    /// The angle between `self` and `o`, robust for both near-parallel and
    /// near-antiparallel vectors (`atan2(|a×b|, a·b)`).
    ///
    /// Proof-backed by `riemann` `Hausdorff/CrossProduct.cross3_normSq_unit`:
    /// for unit vectors `‖a×b‖² = 1 − (a·b)² = sin²c`, so `|a×b|`/`a·b` are the
    /// `sin`/`cos` of the angle and `atan2` recovers it exactly.
    ///
    /// Returns an [`Angle`], which is `f64`-valued, so this is defined only for
    /// the `f64` instantiation of [`Vector`].
    #[must_use]
    pub fn angle(self, o: Self) -> Angle {
        Angle::from_radians(self.cross(o).norm().atan2(self.dot(o)))
    }
}
