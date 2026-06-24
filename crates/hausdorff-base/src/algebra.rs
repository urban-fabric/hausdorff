//! Implementations of the shared [geometry traits](hausdorff_traits)
//! for this crate's substrate types — the ℝ³ [`Vector`](crate::space::r3::Vector)
//! and the bounding regions on the line / circle / ℝ³
//! ([`Interval`](crate::space::r1::Interval),
//! [`Sector`](crate::space::s1::Sector), [`Aabb`](crate::space::r3::Aabb)).
//!
//! These instances live here, beside their types, rather than in the downstream
//! `hausdorff-planar` `algebra` module: the abstract trait surface is upstream
//! (in `hausdorff-traits`), so `hausdorff-base` can implement it
//! directly and its own consumers reach the abstractions via
//! `hausdorff_base::algebra::*` (which re-exports the trait surface, mirroring
//! `hausdorff_planar::algebra` minus that crate's `Distance`-bound metric layer).
//!
//! Each impl's law is machine-checked in `riemann`; the citations sit beside the
//! tests below and in the trait-definition crate's module docs.

pub use hausdorff_traits::{
    Barycentric, InnerProduct, Interpolate, Norm, Orientation, PointRegion, Prunable, Simplex,
    VectorSpace, any_contained, containment_masks, count_contained, disjoint, filter_contained,
    in_triangle, is_collinear, left_of_all, sample, sample_at,
};

use euclid::numeric::{Compare, Element, Real};

use crate::space::r1::Interval;
use crate::space::r3::{Aabb, Vector};
use crate::space::s1::{Angle, Sector};

// ---- R3: the ℝ³ vector — orientation, vector space, interpolation ----

impl<T> Orientation for Vector<T>
where
    T: Real,
{
    type Scalar = T;

    /// Scalar triple product `(a × b) · c` — signed volume of the
    /// parallelepiped; its sign is the right-handed orientation.
    #[inline]
    fn orient(a: Self, b: Self, c: Self) -> T {
        a.cross(b).dot(c)
    }
}

impl<T> VectorSpace for Vector<T>
where
    T: Real,
{
    type Scalar = T;

    fn zero() -> Self {
        Vector::new(T::zero(), T::zero(), T::zero())
    }
    fn add(self, rhs: Self) -> Self {
        Vector::add(self, rhs)
    }
    fn sub(self, rhs: Self) -> Self {
        Vector::sub(self, rhs)
    }
    fn scale(self, factor: T) -> Self {
        self.mul(factor)
    }
    fn select(mask: &<T as Compare>::Mask, on_true: Self, on_false: Self) -> Self {
        Vector::select(mask, on_true, on_false)
    }
}

impl<T> InnerProduct for Vector<T>
where
    T: Real,
{
    fn dot(self, rhs: Self) -> T {
        Vector::dot(self, rhs)
    }
}

impl<T> Norm for Vector<T> where T: Real {}

impl<T> Interpolate for Vector<T>
where
    T: Real,
{
    type Scalar = T;

    /// The R3 chord `lerp` (riemann `chordInterp`).
    #[inline]
    fn lerp(self, q: Self, t: T) -> Self {
        self.clone().add(q.sub(self).mul(t))
    }

    #[inline]
    fn midpoint(self, q: Self) -> Self {
        self.clone().add(q.sub(self).mul(T::half()))
    }
}

// ---- R1: the real-line interval ----

impl PointRegion for Interval {
    type Point = f64;
    type Boolean = bool;

    #[inline]
    fn contains(&self, p: f64) -> bool {
        Interval::contains(*self, p)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        Interval::is_empty(*self)
    }
}

impl Prunable for Interval {
    #[inline]
    fn intersects(&self, other: &Self) -> bool {
        Interval::intersects(*self, *other)
    }
}

// ---- S1: the angular sector (1-D cap) ----

impl PointRegion for Sector {
    type Point = Angle;
    type Boolean = bool;

    #[inline]
    fn contains(&self, p: Angle) -> bool {
        Sector::contains_angle(*self, p)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        Sector::is_empty(*self)
    }
}

impl Prunable for Sector {
    #[inline]
    fn intersects(&self, other: &Self) -> bool {
        Sector::intersects(*self, *other)
    }
}

// ---- R3: the axis-aligned bounding box ----

impl<T> PointRegion for Aabb<T>
where
    T: Real,
{
    type Point = Vector<T>;
    type Boolean = <T as Compare>::Mask;

    #[inline]
    fn contains(&self, p: Vector<T>) -> Self::Boolean {
        Aabb::contains(self, p)
    }

    #[inline]
    fn is_empty(&self) -> Self::Boolean {
        // A normalized `Aabb` (`min ≤ max` per axis) is never empty.
        <Self::Boolean as Element>::zero()
    }
}

impl<T> Prunable for Aabb<T>
where
    T: Real,
{
    #[inline]
    fn intersects(&self, other: &Self) -> Self::Boolean {
        Aabb::overlaps(self, other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_orient_triple_product() {
        let x = Vector::<f64>::new(1.0, 0.0, 0.0);
        let y = Vector::<f64>::new(0.0, 1.0, 0.0);
        let z = Vector::<f64>::new(0.0, 0.0, 1.0);
        assert!(Vector::orient(x, y, z) > 0.0); // right-handed
        assert!(Vector::degenerate(x, y, Vector::<f64>::new(1.0, 1.0, 0.0))); // coplanar
    }

    #[test]
    fn vector_space_and_norm() {
        let v = Vector::<f64>::new(3.0, 4.0, 0.0);
        assert_eq!(
            <Vector as VectorSpace>::zero(),
            Vector::<f64>::new(0.0, 0.0, 0.0)
        );
        assert!((Norm::norm(v) - 5.0).abs() < 1e-12);
        assert!((Norm::norm_squared(v) - 25.0).abs() < 1e-12);
        // `normalized` is a unit vector in the same direction.
        let u = Norm::normalized(v);
        assert!((Norm::norm(u) - 1.0).abs() < 1e-12);
        // The zero vector normalizes to zero (no division by zero).
        assert_eq!(
            Norm::normalized(<Vector as VectorSpace>::zero()),
            Vector::<f64>::new(0.0, 0.0, 0.0)
        );
    }

    #[test]
    fn vector_lerp_endpoints_and_midpoint() {
        let p = Vector::<f64>::new(1.0, 2.0, 3.0);
        let q = Vector::<f64>::new(4.0, 6.0, 8.0);
        assert_eq!(p.lerp(q, 0.0), p);
        assert_eq!(p.lerp(q, 1.0), q);
        assert_eq!(p.midpoint(q), Vector::<f64>::new(2.5, 4.0, 5.5));
        // The generic `sample` walks the chord inclusively.
        let pts = sample(p, q, 3);
        assert_eq!(pts.len(), 4);
        assert_eq!(pts[0], p);
        assert_eq!(pts[3], q);
    }

    #[test]
    fn interval_region_and_pruning() {
        let a = Interval::new(0.0, 10.0);
        let b = Interval::new(5.0, 15.0);
        let c = Interval::new(20.0, 30.0);
        assert!(a.contains(3.0));
        assert!(!a.contains(11.0));
        assert!(!a.is_empty());
        assert!(Interval::EMPTY.is_empty());
        assert!(!disjoint(&a, &b)); // Prunable::intersects via the trait
        assert!(disjoint(&a, &c)); // prune_sound: no shared point
        // The generic filters work on the 1-D interval.
        let contained: u64 = count_contained(&a, &[3.0, 11.0, 9.0, -1.0]);
        assert_eq!(contained, 2);
        assert!(!any_contained(&Interval::EMPTY, &[0.5]));
    }

    #[test]
    fn sector_region_uses_the_proof_backed_membership() {
        let s = Sector::from_center_radius(Angle::from_degrees(350.0), Angle::from_degrees(30.0));
        assert!(s.contains(Angle::from_degrees(10.0)));
        assert!(!s.contains(Angle::from_degrees(90.0)));
        let near = Sector::from_center_radius(Angle::from_degrees(20.0), Angle::from_degrees(20.0));
        let far = Sector::from_center_radius(Angle::from_degrees(160.0), Angle::from_degrees(20.0));
        assert!(!disjoint(&s, &near)); // Prunable::intersects via the trait
        assert!(disjoint(&s, &far)); // arcs_disjoint_of_far
    }

    #[test]
    fn aabb_region_and_pruning() {
        let a = Aabb::<f64>::new(
            Vector::<f64>::new(0.0, 0.0, 0.0),
            Vector::<f64>::new(4.0, 4.0, 4.0),
        );
        let b = Aabb::<f64>::new(
            Vector::<f64>::new(3.0, 3.0, 3.0),
            Vector::<f64>::new(8.0, 8.0, 8.0),
        );
        let c = Aabb::<f64>::new(
            Vector::<f64>::new(5.0, 5.0, 5.0),
            Vector::<f64>::new(9.0, 9.0, 9.0),
        );
        assert!(a.contains(Vector::<f64>::new(2.0, 2.0, 2.0)));
        assert!(!a.is_empty());
        assert!(!disjoint(&a, &b)); // overlaps_of_mem3
        assert!(disjoint(&a, &c)); // prune_sound3
        // The generic filter works in R3 too.
        let pts = [
            Vector::<f64>::new(1.0, 1.0, 1.0),
            Vector::<f64>::new(9.0, 9.0, 9.0),
        ];
        let contained: u64 = count_contained(&a, &pts);
        assert_eq!(contained, 1);
    }
}
