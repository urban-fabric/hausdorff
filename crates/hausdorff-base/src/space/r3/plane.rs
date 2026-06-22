use super::*;

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
    T: Real,
{
    /// The plane with the given `normal` passing through `point`
    /// (`offset = normal·point`).
    #[must_use]
    pub fn from_point_normal(point: Vector<T>, normal: Vector<T>) -> Self {
        let offset = normal.clone().dot(point);
        Self { normal, offset }
    }

    /// The plane through three points, oriented by the right-hand normal
    /// `(b − a) × (c − a)`. Degenerate (zero normal) iff the points are collinear.
    #[must_use]
    pub fn from_points(a: Vector<T>, b: Vector<T>, c: Vector<T>) -> Self {
        let normal = b.sub(a.clone()).cross(c.sub(a.clone()));
        let offset = normal.clone().dot(a);
        Self { normal, offset }
    }

    /// The (unnormalized) plane normal.
    #[inline]
    #[must_use]
    pub fn normal(&self) -> Vector<T> {
        self.normal.clone()
    }

    /// The plane offset `d` in `normal·p = d`.
    #[inline]
    #[must_use]
    pub fn offset(&self) -> T {
        self.offset.clone()
    }

    /// The signed value `normal·p − offset`: zero on the plane, and its **sign is
    /// the side** of `p`. For a [`from_points`](Self::from_points) plane this is
    /// exactly `orient4(a,b,c,p)` (`riemann` `Tetrahedron.orient4`).
    #[inline]
    #[must_use]
    pub fn eval(&self, p: Vector<T>) -> T {
        self.normal.clone().dot(p) - self.offset.clone()
    }

    /// The true **perpendicular signed distance** from `p` to the plane —
    /// [`eval`](Self::eval) divided by `‖normal‖` (positive on the normal's side).
    /// `NaN` for a degenerate (zero-normal) plane.
    #[inline]
    #[must_use]
    pub fn signed_distance(&self, p: Vector<T>) -> T {
        self.eval(p) / self.normal.clone().norm()
    }
}

impl<T> Plane<T>
where
    T: Real + PartialOrd,
{
    /// Which side of the plane `p` lies on: `+1` on the normal's side, `-1` on the
    /// far side, `0` exactly on the plane — the sign of [`eval`](Self::eval).
    #[inline]
    #[must_use]
    pub fn side_of(&self, p: Vector<T>) -> i32 {
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
    pub fn contains(&self, p: Vector<T>, eps: T) -> bool {
        self.eval(p).abs() <= eps
    }
}

impl<T> Plane<T>
where
    T: Real,
{
    /// Mask-valued plane membership within an absolute `eval` tolerance `eps`.
    #[inline]
    #[must_use]
    pub fn contains_mask(&self, p: Vector<T>, eps: T) -> <T as Compare>::Mask {
        self.eval(p).abs().le(&eps)
    }
}
