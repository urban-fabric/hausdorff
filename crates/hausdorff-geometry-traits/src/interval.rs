//! [`ClosedInterval`] — the closed-interval algebra shared by `r1::Interval`
//! (the real line) and `s2::LatInterval` (the `Rect` latitude coordinate), both
//! ports of golang/geo's `r1.Interval`.
//!
//! A type provides only the constructor [`from_bounds`](ClosedInterval::from_bounds)
//! and the [`lo`](ClosedInterval::lo) / [`hi`](ClosedInterval::hi) accessors; the
//! containment, intersection, union, expansion, and two-point hull come free —
//! written once here, the minimal structure the interval algebra actually needs.

/// A closed interval `[lo, hi]` on an `f64` axis, with `lo > hi` the empty
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
pub trait ClosedInterval: Copy + Sized {
    /// The interval `[lo, hi]`.
    fn from_bounds(lo: f64, hi: f64) -> Self;

    /// The low endpoint.
    fn lo(self) -> f64;

    /// The high endpoint.
    fn hi(self) -> f64;

    /// The canonical empty interval `[1, 0]`.
    #[must_use]
    fn empty() -> Self {
        Self::from_bounds(1.0, 0.0)
    }

    /// True iff the interval is empty (`lo > hi`).
    #[must_use]
    fn is_empty(self) -> bool {
        self.lo() > self.hi()
    }

    /// True iff `p ∈ [lo, hi]` (endpoints included).
    #[must_use]
    fn contains(self, p: f64) -> bool {
        self.lo() <= p && p <= self.hi()
    }

    /// True iff `other` is wholly contained in `self` (the empty interval is
    /// contained in everything).
    #[must_use]
    fn contains_interval(self, other: Self) -> bool {
        other.is_empty() || (self.lo() <= other.lo() && other.hi() <= self.hi())
    }

    /// The intersection of `self` and `other` (empty if they are disjoint).
    #[must_use]
    fn intersection(self, other: Self) -> Self {
        Self::from_bounds(self.lo().max(other.lo()), self.hi().min(other.hi()))
    }

    /// The smallest interval containing both `self` and `other` (the other if one
    /// is empty).
    #[must_use]
    fn union(self, other: Self) -> Self {
        if self.is_empty() {
            other
        } else if other.is_empty() {
            self
        } else {
            Self::from_bounds(self.lo().min(other.lo()), self.hi().max(other.hi()))
        }
    }

    /// `self` grown (negative `margin` shrinks) by `margin` on each side; the
    /// empty interval stays empty.
    #[must_use]
    fn expanded(self, margin: f64) -> Self {
        if self.is_empty() {
            self
        } else {
            Self::from_bounds(self.lo() - margin, self.hi() + margin)
        }
    }

    /// The smallest interval containing both points.
    #[must_use]
    fn from_point_pair(a: f64, b: f64) -> Self {
        Self::from_bounds(a.min(b), a.max(b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Copy, Clone, Debug, PartialEq)]
    struct I {
        lo: f64,
        hi: f64,
    }
    impl ClosedInterval for I {
        fn from_bounds(lo: f64, hi: f64) -> Self {
            Self { lo, hi }
        }
        fn lo(self) -> f64 {
            self.lo
        }
        fn hi(self) -> f64 {
            self.hi
        }
    }

    #[test]
    fn interval_algebra_written_once() {
        let a = I::from_bounds(0.0, 2.0);
        let b = I::from_bounds(1.0, 3.0);
        // intersection is set intersection
        let i = a.intersection(b);
        assert_eq!(i, I::from_bounds(1.0, 2.0));
        assert!(i.contains(1.5) && a.contains(1.5) && b.contains(1.5));
        // union contains both, empty-absorbing
        let u = a.union(b);
        assert_eq!(u, I::from_bounds(0.0, 3.0));
        assert_eq!(a.union(I::empty()), a);
        assert!(I::empty().is_empty() && !a.is_empty());
        // contains_interval: empty ⊆ everything; soundness
        assert!(a.contains_interval(I::empty()));
        assert!(u.contains_interval(a) && u.contains_interval(b));
        // expanded only grows; from_point_pair is the min/max hull
        assert!(a.expanded(1.0).contains_interval(a));
        assert_eq!(I::from_point_pair(3.0, -1.0), I::from_bounds(-1.0, 3.0));
    }
}
