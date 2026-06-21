//! `r3` — three-dimensional Euclidean geometry: the [`Vector`] type, mirroring
//! golang/geo's `r3` package.
//!
//! `r3::Vector` is the substrate of the spherical world: an `S2Point` (in `hausdorff-sphere`)
//! is a *unit* `r3::Vector`, and spherical measures/predicates are computed from
//! vector dot/cross products rather than lat/lon, which is what makes them free
//! of pole and antimeridian singularities.

use euclid::numeric::LazyReal;

use crate::space::s1::Angle;

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
    T: LazyReal + Copy + PartialOrd,
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
    pub fn mul(self, m: T) -> Self {
        Self::new(self.x * m, self.y * m, self.z * m)
    }

    /// Squared Euclidean norm (`self · self`).
    #[inline]
    #[must_use]
    pub fn norm2(self) -> T {
        self.dot(self)
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
        let n2 = self.norm2();
        if n2 == T::zero() {
            return Self::new(T::zero(), T::zero(), T::zero());
        }
        self.mul(T::one() / n2.sqrt())
    }

    /// Componentwise absolute value.
    #[must_use]
    pub fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// True iff the vector is (near) unit length.
    #[must_use]
    pub fn is_unit(self) -> bool {
        (self.norm2() - T::one()).abs() <= T::scalar(1e-14)
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
    pub fn cos_angle(self, other: Self) -> T {
        let denom = self.norm() * other.norm();
        if denom == T::zero() {
            return T::zero();
        }
        // Cauchy–Schwarz guarantees the ratio is in [-1, 1]; the explicit clamp
        // only guards against f64 rounding nudging it a hair outside.
        (self.dot(other) / denom).clamp(-T::one(), T::one())
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
            Self::new(T::one(), T::zero(), T::zero())
        } else if ty < tz {
            Self::new(T::zero(), T::one(), T::zero())
        } else {
            Self::new(T::zero(), T::zero(), T::one())
        };
        self.cross(other).normalize()
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

/// An **axis-aligned 3-D bounding box**, stored by its `min` (lowest x, y, z) and
/// `max` (highest) corners — the ℝ³ member of the Euclidean box ladder alongside
/// [`r1::Interval`](crate::space::r1::Interval) (1-D) and `r2::Rect` (2-D, in
/// `hausdorff-planar`), and the box an R-tree over `Vector` points prunes with.
///
/// Membership and the overlap / containment pruning tests are per-axis interval
/// arithmetic, proof-backed by `riemann` `ExactArith/AABB` (`AABB3`):
/// [`overlaps`](Self::overlaps) by `overlaps_of_mem3` (a shared point forces an
/// overlap, so a non-overlapping pair is safe to prune — `prune_sound3`),
/// [`contains_box`](Self::contains_box) by `contains_of_subset3` (containment is
/// monotone) with `subset_trans3` carrying the leaf⊆…⊆root invariant.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Aabb<T = f64> {
    min: Vector<T>,
    max: Vector<T>,
}

impl<T> Aabb<T>
where
    T: LazyReal + Copy + PartialOrd,
{
    /// The box spanning two opposite corners, given in any order (normalized to
    /// per-axis `min`/`max`).
    #[must_use]
    pub fn new(c1: Vector<T>, c2: Vector<T>) -> Self {
        Self {
            min: Vector::new(c1.x.min(c2.x), c1.y.min(c2.y), c1.z.min(c2.z)),
            max: Vector::new(c1.x.max(c2.x), c1.y.max(c2.y), c1.z.max(c2.z)),
        }
    }

    /// The box from already-ordered corners (`min ≤ max` per axis assumed).
    #[must_use]
    pub fn from_min_max(min: Vector<T>, max: Vector<T>) -> Self {
        Self { min, max }
    }

    /// The degenerate box containing only the point `p`.
    #[must_use]
    pub fn from_point(p: Vector<T>) -> Self {
        Self { min: p, max: p }
    }

    /// The lower corner (minimum x, y, z).
    #[inline]
    #[must_use]
    pub fn min(self) -> Vector<T> {
        self.min
    }

    /// The upper corner (maximum x, y, z).
    #[inline]
    #[must_use]
    pub fn max(self) -> Vector<T> {
        self.max
    }

    /// The centre point.
    #[inline]
    #[must_use]
    pub fn center(self) -> Vector<T> {
        self.min.add(self.max).mul(T::half())
    }

    /// True iff `p` lies in the box (edges inclusive) — the per-axis range test
    /// (`riemann` `AABB3.contains`).
    #[inline]
    #[must_use]
    pub fn contains(self, p: Vector<T>) -> bool {
        self.min.x <= p.x
            && p.x <= self.max.x
            && self.min.y <= p.y
            && p.y <= self.max.y
            && self.min.z <= p.z
            && p.z <= self.max.z
    }

    /// True iff this box and `other` overlap (share any point, faces inclusive) —
    /// the per-axis interval-overlap test. A non-overlapping pair shares no point,
    /// so a spatial index may prune it (`riemann` `AABB3.overlaps_of_mem3` /
    /// `prune_sound3`).
    #[inline]
    #[must_use]
    pub fn overlaps(self, other: Aabb<T>) -> bool {
        self.min.x <= other.max.x
            && other.min.x <= self.max.x
            && self.min.y <= other.max.y
            && other.min.y <= self.max.y
            && self.min.z <= other.max.z
            && other.min.z <= self.max.z
    }

    /// True iff this box fully contains `other` (`riemann` `AABB3.subset`, with
    /// `contains_of_subset3` ⟹ every point of `other` is in `self`).
    #[inline]
    #[must_use]
    pub fn contains_box(self, other: Aabb<T>) -> bool {
        self.min.x <= other.min.x
            && other.max.x <= self.max.x
            && self.min.y <= other.min.y
            && other.max.y <= self.max.y
            && self.min.z <= other.min.z
            && other.max.z <= self.max.z
    }

    /// The smallest box containing both `self` and `other`.
    #[must_use]
    pub fn union(self, other: Aabb<T>) -> Aabb<T> {
        Aabb {
            min: Vector::new(
                self.min.x.min(other.min.x),
                self.min.y.min(other.min.y),
                self.min.z.min(other.min.z),
            ),
            max: Vector::new(
                self.max.x.max(other.max.x),
                self.max.y.max(other.max.y),
                self.max.z.max(other.max.z),
            ),
        }
    }

    /// The box grown by `margin` on every side (a negative `margin` shrinks it).
    #[must_use]
    pub fn expanded(self, margin: T) -> Aabb<T> {
        let d = Vector::new(margin, margin, margin);
        Aabb {
            min: self.min.sub(d),
            max: self.max.add(d),
        }
    }
}

/// A **plane** in ℝ³, stored as a `normal` and an `offset` so the plane is the
/// set `{ p : normal·p = offset }`. The 3-D analogue of an R2 line; the splitting
/// surface behind half-space tests, clipping, and BSP/k-d partitioning.
///
/// The signed [`eval`](Self::eval) `normal·p − offset` is zero on the plane and
/// its **sign is the side**. For a plane through three points `a, b, c`
/// ([`from_points`](Self::from_points)), `eval(d)` equals the tetrahedron
/// orientation `orient4(a,b,c,d) = ((b−a)×(c−a))·(d−a)` exactly, so the side-test
/// is the proof-backed `riemann` `ExactArith/Tetrahedron.orient4`; and because
/// `orient4` is affine along a segment in its last vertex
/// (`orient4_lerp_last`), opposite signs at two points certify the segment
/// crosses the plane (the 3-D straddle test).
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Plane<T = f64> {
    normal: Vector<T>,
    offset: T,
}

impl<T> Plane<T>
where
    T: LazyReal + Copy + PartialOrd,
{
    /// The plane with the given `normal` passing through `point`
    /// (`offset = normal·point`).
    #[must_use]
    pub fn from_point_normal(point: Vector<T>, normal: Vector<T>) -> Self {
        Self {
            normal,
            offset: normal.dot(point),
        }
    }

    /// The plane through three points, oriented by the right-hand normal
    /// `(b − a) × (c − a)`. Degenerate (zero normal) iff the points are collinear.
    #[must_use]
    pub fn from_points(a: Vector<T>, b: Vector<T>, c: Vector<T>) -> Self {
        let normal = b.sub(a).cross(c.sub(a));
        Self {
            normal,
            offset: normal.dot(a),
        }
    }

    /// The (unnormalized) plane normal.
    #[inline]
    #[must_use]
    pub fn normal(self) -> Vector<T> {
        self.normal
    }

    /// The plane offset `d` in `normal·p = d`.
    #[inline]
    #[must_use]
    pub fn offset(self) -> T {
        self.offset
    }

    /// The signed value `normal·p − offset`: zero on the plane, and its **sign is
    /// the side** of `p`. For a [`from_points`](Self::from_points) plane this is
    /// exactly `orient4(a,b,c,p)` (`riemann` `Tetrahedron.orient4`).
    #[inline]
    #[must_use]
    pub fn eval(self, p: Vector<T>) -> T {
        self.normal.dot(p) - self.offset
    }

    /// The true **perpendicular signed distance** from `p` to the plane —
    /// [`eval`](Self::eval) divided by `‖normal‖` (positive on the normal's side).
    /// `NaN` for a degenerate (zero-normal) plane.
    #[inline]
    #[must_use]
    pub fn signed_distance(self, p: Vector<T>) -> T {
        self.eval(p) / self.normal.norm()
    }

    /// Which side of the plane `p` lies on: `+1` on the normal's side, `-1` on the
    /// far side, `0` exactly on the plane — the sign of [`eval`](Self::eval).
    #[inline]
    #[must_use]
    pub fn side_of(self, p: Vector<T>) -> i32 {
        let e = self.eval(p);
        if e > T::zero() {
            1
        } else if e < T::zero() {
            -1
        } else {
            0
        }
    }

    /// True iff `p` lies on the plane within an absolute `eval` tolerance `eps`.
    #[inline]
    #[must_use]
    pub fn contains(self, p: Vector<T>, eps: T) -> bool {
        self.eval(p).abs() <= eps
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_is_generic_over_the_real_scalar_f32_parity() {
        // The ℝ³ `Vector` is no longer `f64`-locked: the same operations run at
        // `f32` and agree with the `f64` reference within single-precision
        // tolerance. This is the proof the geometry is generic over the
        // `LazyReal` scalar trait.
        let a64 = Vector::<f64>::new(1.0, 2.0, 3.0);
        let b64 = Vector::<f64>::new(-4.0, 0.5, 6.0);
        let c64 = Vector::<f64>::new(0.0, -1.0, 2.0);
        let a32 = Vector::<f32>::new(1.0, 2.0, 3.0);
        let b32 = Vector::<f32>::new(-4.0, 0.5, 6.0);
        let c32 = Vector::<f32>::new(0.0, -1.0, 2.0);

        // dot product
        assert!((f64::from(a32.dot(b32)) - a64.dot(b64)).abs() < 1e-5);
        // norm
        assert!((f64::from(a32.norm()) - a64.norm()).abs() < 1e-5);
        // cross product, componentwise
        let x32 = a32.cross(b32);
        let x64 = a64.cross(b64);
        assert!((f64::from(x32.x) - x64.x).abs() < 1e-4);
        assert!((f64::from(x32.y) - x64.y).abs() < 1e-4);
        assert!((f64::from(x32.z) - x64.z).abs() < 1e-4);
        // normalized vector is unit length at f32 too
        let u32 = a32.normalize();
        assert!((f64::from(u32.norm()) - 1.0).abs() < 1e-6);
        // the scalar triple product / orientation predicate
        let o32 = a32.cross(b32).dot(c32);
        let o64 = a64.cross(b64).dot(c64);
        assert!((f64::from(o32) - o64).abs() < 1e-3);
    }

    #[test]
    fn project_onto_tangent_is_perpendicular_with_pythagorean_length() {
        let v = Vector::<f64>::new(2.0, -3.0, 5.0);
        let n = Vector::<f64>::new(1.0, 1.0, 1.0).normalize();
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
        let vp = Vector::<f64>::new(0.3, 0.7, -0.4);
        let n = Vector::<f64>::new(-0.2, 0.5, 0.84).normalize();
        let inline = vp.sub(n.mul(vp.dot(n)));
        assert_eq!(vp.project_onto_tangent(n), inline);
        // The latitude-extreme vertex ẑ − n·n_z is ẑ.project_onto_tangent(n).
        let z = Vector::<f64>::new(0.0, 0.0, 1.0);
        let vertex = Vector::<f64>::new(-n.x * n.z, -n.y * n.z, 1.0 - n.z * n.z);
        let proj = z.project_onto_tangent(n);
        assert!((proj.x - vertex.x).abs() < 1e-15);
        assert!((proj.y - vertex.y).abs() < 1e-15);
        assert!((proj.z - vertex.z).abs() < 1e-15);
    }

    #[test]
    fn cos_angle_known_configurations_and_range() {
        let x = Vector::<f64>::new(1.0, 0.0, 0.0);
        let y = Vector::<f64>::new(0.0, 1.0, 0.0);
        // Orthogonal → 0.
        assert!(x.cos_angle(y).abs() < 1e-15);
        // Same direction (any positive scale) → 1.
        assert!((x.cos_angle(x.mul(7.0)) - 1.0).abs() < 1e-15);
        // Opposite → −1.
        assert!((x.cos_angle(x.mul(-3.0)) + 1.0).abs() < 1e-15);
        // 45° between (1,0,0) and (1,1,0) → 1/√2.
        let d = Vector::<f64>::new(1.0, 1.0, 0.0);
        assert!((x.cos_angle(d) - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-15);
        // Zero vector is graceful.
        assert_eq!(x.cos_angle(Vector::<f64>::new(0.0, 0.0, 0.0)), 0.0);
        // Always in [-1, 1] (Cauchy–Schwarz), for assorted vectors.
        let vecs = [
            Vector::<f64>::new(3.0, -1.0, 2.0),
            Vector::<f64>::new(-4.0, 5.0, 0.5),
            Vector::<f64>::new(0.1, 0.2, -9.0),
            Vector::<f64>::new(-2.0, -2.0, -2.0),
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
                Vector::<f64>::new(1.0, 0.0, 0.0),
                Vector::<f64>::new(0.0, 1.0, 0.0),
                Vector::<f64>::new(0.0, 0.0, 1.0),
            ),
            (
                Vector::<f64>::new(2.0, -1.0, 3.0),
                Vector::<f64>::new(0.5, 4.0, -2.0),
                Vector::<f64>::new(-3.0, 1.0, 1.0),
            ),
            (
                Vector::<f64>::new(1.0, 1.0, 0.0),
                Vector::<f64>::new(0.0, 1.0, 1.0),
                Vector::<f64>::new(1.0, 0.0, 1.0),
            ),
            // A coplanar (degenerate) triple → triple product 0 → gram_det 0.
            (
                Vector::<f64>::new(1.0, 2.0, 3.0),
                Vector::<f64>::new(2.0, 4.0, 6.0),
                Vector::<f64>::new(-1.0, 0.0, 1.0),
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
            Vector::<f64>::new(0.3, -0.5, 0.8).normalize(),
            Vector::<f64>::new(-0.7, 0.2, 0.4).normalize(),
            Vector::<f64>::new(0.1, 0.9, -0.2).normalize(),
        ];
        let (a, b, c) = (unit[0], unit[1], unit[2]);
        let (ab, bc, ca) = (a.dot(b), b.dot(c), c.dot(a));
        let spherical = 1.0 - ab * ab - bc * bc - ca * ca + 2.0 * ab * bc * ca;
        assert!((Vector::gram_det(a, b, c) - spherical).abs() < 1e-12);
    }

    #[test]
    fn aabb_contains_overlaps_and_nesting() {
        // Corners accepted in any order and normalized.
        let b = Aabb::<f64>::new(
            Vector::<f64>::new(4.0, 0.0, 6.0),
            Vector::<f64>::new(0.0, 4.0, 2.0),
        );
        assert_eq!(b.min(), Vector::<f64>::new(0.0, 0.0, 2.0));
        assert_eq!(b.max(), Vector::<f64>::new(4.0, 4.0, 6.0));
        assert_eq!(b.center(), Vector::<f64>::new(2.0, 2.0, 4.0));
        // Membership (AABB3.contains): inside, on a face, and outside (per axis).
        assert!(b.contains(Vector::<f64>::new(2.0, 2.0, 4.0)));
        assert!(b.contains(Vector::<f64>::new(0.0, 4.0, 6.0))); // corner, inclusive
        assert!(!b.contains(Vector::<f64>::new(2.0, 2.0, 7.0))); // past max z
        // Overlap (overlaps_of_mem3) vs. a disjoint box (prune_sound3).
        let touching = Aabb::<f64>::new(
            Vector::<f64>::new(4.0, 4.0, 6.0),
            Vector::<f64>::new(6.0, 6.0, 8.0),
        );
        let far = Aabb::<f64>::new(
            Vector::<f64>::new(10.0, 10.0, 10.0),
            Vector::<f64>::new(12.0, 12.0, 12.0),
        );
        assert!(b.overlaps(touching)); // shares the corner
        assert!(!b.overlaps(far));
        // Nesting (contains_of_subset3) and its transitive union.
        let inner = Aabb::<f64>::new(
            Vector::<f64>::new(1.0, 1.0, 3.0),
            Vector::<f64>::new(3.0, 3.0, 5.0),
        );
        assert!(b.contains_box(inner));
        assert!(!inner.contains_box(b));
        assert!(b.union(far).contains_box(b));
        assert!(b.union(far).contains_box(far));
        // Expansion grows every face.
        assert!(
            b.expanded(1.0)
                .contains(Vector::<f64>::new(-0.5, -0.5, 1.5))
        );
    }

    #[test]
    fn plane_side_test_equals_orient4() {
        // The z = 0 plane (normal +z) via three points.
        let (a, b, c) = (
            Vector::<f64>::new(0.0, 0.0, 0.0),
            Vector::<f64>::new(1.0, 0.0, 0.0),
            Vector::<f64>::new(0.0, 1.0, 0.0),
        );
        let plane = Plane::from_points(a, b, c);
        assert_eq!(plane.normal(), Vector::<f64>::new(0.0, 0.0, 1.0));
        // Sides: above (+z), below (−z), and on the plane.
        assert_eq!(plane.side_of(Vector::<f64>::new(0.5, 0.5, 2.0)), 1);
        assert_eq!(plane.side_of(Vector::<f64>::new(0.5, 0.5, -2.0)), -1);
        assert_eq!(plane.side_of(Vector::<f64>::new(0.5, 0.5, 0.0)), 0);
        assert!(plane.contains(Vector::<f64>::new(9.0, -3.0, 0.0), 1e-12));
        // Perpendicular signed distance is the bare z (unit normal).
        assert!((plane.signed_distance(Vector::<f64>::new(0.0, 0.0, 3.0)) - 3.0).abs() < 1e-12);

        // `eval(d)` equals orient4(a, b, c, d) = ((b−a)×(c−a))·(d−a) exactly.
        let orient4 =
            |a: Vector, b: Vector, c: Vector, d: Vector| b.sub(a).cross(c.sub(a)).dot(d.sub(a));
        for d in [
            Vector::<f64>::new(0.3, 0.4, 1.7),
            Vector::<f64>::new(-2.0, 1.0, -0.5),
            Vector::<f64>::new(5.0, 5.0, 0.0),
        ] {
            assert!((plane.eval(d) - orient4(a, b, c, d)).abs() < 1e-12);
        }
    }

    #[test]
    fn plane_from_point_normal_and_straddle() {
        let plane = Plane::from_point_normal(
            Vector::<f64>::new(0.0, 0.0, 1.0),
            Vector::<f64>::new(0.0, 0.0, 2.0),
        );
        // Offset = normal·point = 2; the plane is z = 1.
        assert_eq!(plane.offset(), 2.0);
        assert_eq!(plane.side_of(Vector::<f64>::new(0.0, 0.0, 5.0)), 1);
        assert_eq!(plane.side_of(Vector::<f64>::new(0.0, 0.0, 0.0)), -1);
        // Opposite sides ⟹ the segment crosses the plane (orient4_lerp_last):
        // eval is affine along the segment, so it has a zero in (0, 1).
        let p = Vector::<f64>::new(0.0, 0.0, 0.0);
        let q = Vector::<f64>::new(0.0, 0.0, 4.0);
        assert!(plane.side_of(p) * plane.side_of(q) < 0);
        let t = -plane.eval(p) / (plane.eval(q) - plane.eval(p));
        assert!((0.0..=1.0).contains(&t));
        let crossing = p.add(q.sub(p).mul(t));
        assert!(plane.contains(crossing, 1e-12));
    }
}
