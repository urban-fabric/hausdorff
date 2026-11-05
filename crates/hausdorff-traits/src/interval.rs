//! [`ClosedInterval`] — the closed-interval algebra shared by `r1::Interval`
//! (the real line) and `s2::LatInterval` (the `Rect` latitude coordinate), both
//! ports of golang/geo's `r1.Interval`.
//!
//! A type provides the constructor [`from_bounds`](ClosedInterval::from_bounds),
//! the [`lo`](ClosedInterval::lo) / [`hi`](ClosedInterval::hi) accessors, and
//! mask-driven selection; the containment, intersection, union, expansion, and
//! two-point hull come free — written once here, the minimal structure the
//! interval algebra actually needs.

use euclid_numeric::{Boolean as BooleanMask, Compare, Element, Real};

use crate::orientation::Boolean;

#[inline]
fn mask_and<M: BooleanMask>(a: M, b: M) -> M {
    a.and(&b)
}

#[inline]
fn mask_or<M: BooleanMask>(a: M, b: M) -> M {
    a.or(&b)
}

/// A closed interval `[lo, hi]` on a [`Real`] axis, with `lo > hi` the empty
/// interval (canonical empty `[1, 0]`).
///
/// The laws are machine-checked in `riemann` — `Hausdorff/R1Interval` and
/// `Hausdorff/SphereRect` (`LatInterval`):
/// - `contains` / `contains_interval` are sound — a passing containment test
///   carries membership (`mem_of_subset` / `subset_sound`), with the empty
///   interval contained in everything;
/// - `intersection` is exactly set intersection (`mem_inter_iff`);
/// - `union` contains both inputs and is empty-absorbing
///   (`mem_hull_of_mem_left` / `mem_hull_of_mem_right`, `union_mem_*`);
/// - `expanded` by a non-negative margin only grows the interval
///   (`expanded_mem_of_mem`), capturing the whole margin band
///   (`expanded_mem_of_near`);
/// - `from_point_pair` is the smallest interval through both points
///   (`fromPointPair_minimal`, with `fromPointPair_mem_left/right`).
pub trait ClosedInterval: Clone + Sized {
    /// The scalar field the endpoints live on.
    type Scalar: Real;

    /// The interval `[lo, hi]`.
    fn from_bounds(lo: Self::Scalar, hi: Self::Scalar) -> Self;

    /// The low endpoint.
    fn lo(&self) -> Self::Scalar;

    /// The high endpoint.
    fn hi(&self) -> Self::Scalar;

    /// Branchless interval selection driven by a scalar comparison mask.
    fn select(mask: &Boolean<Self::Scalar>, on_true: Self, on_false: Self) -> Self;

    /// The canonical empty interval `[1, 0]`.
    #[must_use]
    fn empty() -> Self {
        Self::from_bounds(Self::Scalar::one(), Self::Scalar::zero())
    }

    /// True iff the interval is empty (`lo > hi`).
    #[must_use]
    fn is_empty(&self) -> Boolean<Self::Scalar> {
        self.lo().gt(&self.hi())
    }

    /// True iff `p ∈ [lo, hi]` (endpoints included).
    #[must_use]
    fn contains(&self, p: Self::Scalar) -> Boolean<Self::Scalar> {
        mask_and(self.lo().le(&p), p.le(&self.hi()))
    }

    /// True iff `other` is wholly contained in `self` (the empty interval is
    /// contained in everything).
    #[must_use]
    fn contains_interval(&self, other: &Self) -> Boolean<Self::Scalar> {
        mask_or(
            other.is_empty(),
            mask_and(self.lo().le(&other.lo()), other.hi().le(&self.hi())),
        )
    }

    /// The intersection of `self` and `other` (empty if they are disjoint).
    #[must_use]
    fn intersection(&self, other: &Self) -> Self {
        Self::from_bounds(self.lo().max(other.lo()), self.hi().min(other.hi()))
    }

    /// The smallest interval containing both `self` and `other` (the other if one
    /// is empty).
    #[must_use]
    fn union(&self, other: &Self) -> Self {
        let self_empty = self.is_empty();
        let other_empty = other.is_empty();
        let hull = Self::from_bounds(self.lo().min(other.lo()), self.hi().max(other.hi()));
        let when_self_nonempty = Self::select(&other_empty, self.clone(), hull);
        Self::select(&self_empty, other.clone(), when_self_nonempty)
    }

    /// `self` grown (negative `margin` shrinks) by `margin` on each side; the
    /// empty interval stays empty.
    #[must_use]
    fn expanded(&self, margin: Self::Scalar) -> Self {
        let expanded = Self::from_bounds(self.lo() - margin.clone(), self.hi() + margin);
        Self::select(&self.is_empty(), self.clone(), expanded)
    }

    /// The smallest interval containing both points.
    #[must_use]
    fn from_point_pair(a: Self::Scalar, b: Self::Scalar) -> Self {
        Self::from_bounds(a.clone().min(b.clone()), a.max(b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Debug, PartialEq)]
    struct I {
        lo: f64,
        hi: f64,
    }
    impl ClosedInterval for I {
        type Scalar = f64;
        fn from_bounds(lo: f64, hi: f64) -> Self {
            Self { lo, hi }
        }
        fn lo(&self) -> f64 {
            self.lo
        }
        fn hi(&self) -> f64 {
            self.hi
        }
        fn select(mask: &bool, on_true: Self, on_false: Self) -> Self {
            if *mask { on_true } else { on_false }
        }
    }

    #[test]
    fn interval_algebra_written_once() {
        let a = I::from_bounds(0.0, 2.0);
        let b = I::from_bounds(1.0, 3.0);
        // intersection is set intersection
        let i = a.intersection(&b);
        assert_eq!(i, I::from_bounds(1.0, 2.0));
        assert!(i.contains(1.5) && a.contains(1.5) && b.contains(1.5));
        // union contains both, empty-absorbing
        let u = a.union(&b);
        assert_eq!(u, I::from_bounds(0.0, 3.0));
        assert_eq!(a.union(&I::empty()), a);
        assert!(I::empty().is_empty() && !a.is_empty());
        // contains_interval: empty ⊆ everything; soundness
        assert!(a.contains_interval(&I::empty()));
        assert!(u.contains_interval(&a) && u.contains_interval(&b));
        // expanded only grows; from_point_pair is the min/max hull
        assert!(a.expanded(1.0).contains_interval(&a));
        assert_eq!(I::from_point_pair(3.0, -1.0), I::from_bounds(-1.0, 3.0));
    }
}
