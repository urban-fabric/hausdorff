//! Spatial regions — the bounding-volume abstraction shared by every space's
//! `O(1)` containment primitive, and the pruning predicate spatial indexes ride.
//!
//! Each space carries a constant-size region with constant-time point
//! membership: [`r1::Interval`](hausdorff_base::space::r1::Interval) on the line,
//! [`s1::Sector`](hausdorff_base::space::s1::Sector) on the circle,
//! [`Rect`](crate::space::r2::Rect) / [`Ellipse`](crate::space::r2::Ellipse) /
//! [`ConfidenceEllipse`](crate::space::r2::ConfidenceEllipse) in the plane, and
//! `Cap` on the sphere (in `hausdorff-sphere`). They differ only in their point
//! type, so [`PointRegion`] factors out membership once and generic spatial
//! filters ([`filter_contained`], [`count_contained`], [`any_contained`]) are
//! written against the trait rather than per space.
//!
//! [`Prunable`] adds the second indexing operation — a pairwise
//! [`intersects`](Prunable::intersects) test with **no false negatives**, so a
//! `false` is always safe to prune on. It is exact for the interval / metric-ball
//! family ([`Interval`], [`Sector`], [`Rect`], [`Aabb`](hausdorff_base::space::r3::Aabb),
//! and `Cap`); regions whose pairwise intersection has no closed form (the
//! ellipses) implement only [`PointRegion`].
//!
//! The laws behind each impl are machine-checked in `riemann`:
//! - membership — `Hausdorff/Ellipse.inEllipse`, `SphereCap.inCap_iff_dot`,
//!   `Circular.inArc_iff_dot`, and interval/box order facts;
//! - pruning — `ExactArith/AABB.{overlaps_of_mem, prune_sound}` (the R-tree
//!   soundness), `SphereCap.caps_disjoint_of_far`, `Circular.arcs_disjoint_of_far`
//!   (far-apart centres share no point — exactly `!intersects`).

use hausdorff_base::space::r1::Interval;
use hausdorff_base::space::r3::{Aabb, Vector};
use hausdorff_base::space::s1::{Angle, Sector};

use crate::space::r2::{ConfidenceEllipse, Coord, Ellipse, Rect};

/// A spatial region: a set of points of type [`Point`](Self::Point) with
/// constant-time membership. The cross-space generalization of every
/// `contains_*` primitive (`contains_coord`, `contains_point`, `contains_angle`,
/// …). Each impl's membership test is proof-backed in `riemann` (see the module
/// docs).
pub trait PointRegion {
    /// The point type this region is tested against (`f64` for an interval,
    /// `Angle` for a sector, `Coord` for the planar regions, `S2Point` for a cap).
    type Point;

    /// True iff `p` lies in the region (boundary inclusive).
    fn contains(&self, p: Self::Point) -> bool;

    /// True iff the region contains no points.
    fn is_empty(&self) -> bool;
}

/// A [`PointRegion`] that supports a pairwise **intersection** test — the
/// bounding-volume pruning predicate behind R-tree / cap / sector indexes.
///
/// The test never yields a false negative: if two regions share any point it
/// returns `true`, so a `false` result certifies they are disjoint and the pair
/// may be pruned. For the interval / metric-ball family it is *exact*, backed by
/// the `riemann` disjointness proofs (`AABB.prune_sound`,
/// `caps_disjoint_of_far`, `arcs_disjoint_of_far`).
pub trait Prunable: PointRegion {
    /// True iff `self` and `other` may share a point.
    fn intersects(&self, other: &Self) -> bool;
}

/// The points of `points` that lie in `region` — a generic spatial filter over
/// any [`PointRegion`], written once for every space.
#[must_use]
pub fn filter_contained<R>(region: &R, points: &[R::Point]) -> Vec<R::Point>
where
    R: PointRegion,
    R::Point: Copy,
{
    points
        .iter()
        .copied()
        .filter(|&p| region.contains(p))
        .collect()
}

/// How many of `points` lie in `region`.
#[must_use]
pub fn count_contained<R>(region: &R, points: &[R::Point]) -> usize
where
    R: PointRegion,
    R::Point: Copy,
{
    points.iter().copied().filter(|&p| region.contains(p)).count()
}

/// True iff any of `points` lies in `region` (`false` for an empty region or no
/// points).
#[must_use]
pub fn any_contained<R>(region: &R, points: &[R::Point]) -> bool
where
    R: PointRegion,
    R::Point: Copy,
{
    !region.is_empty() && points.iter().copied().any(|p| region.contains(p))
}

/// True iff `a` and `b` provably share no point, so an index may prune the pair.
/// The certified contrapositive of [`Prunable::intersects`] — the `riemann`
/// `*_disjoint_of_far` / `prune_sound` soundness.
#[must_use]
pub fn disjoint<R: Prunable>(a: &R, b: &R) -> bool {
    !a.intersects(b)
}

// ---- R1: the real-line interval ----

impl PointRegion for Interval {
    type Point = f64;

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

// ---- R2: the axis-aligned rectangle ----

impl PointRegion for Rect {
    type Point = Coord;

    #[inline]
    fn contains(&self, p: Coord) -> bool {
        Rect::contains_coord(*self, p)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        // A `Rect` is normalized to `min ≤ max`, so it is never empty (a
        // degenerate point/segment still contains its locus).
        false
    }
}

impl Prunable for Rect {
    #[inline]
    fn intersects(&self, other: &Self) -> bool {
        Rect::overlaps(*self, *other)
    }
}

// ---- R2: the elliptical disk and the confidence ellipse (membership only) ----

impl PointRegion for Ellipse {
    type Point = Coord;

    #[inline]
    fn contains(&self, p: Coord) -> bool {
        Ellipse::contains_coord(*self, p)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.semi_x() <= 0.0 || self.semi_y() <= 0.0
    }
}

impl PointRegion for ConfidenceEllipse {
    type Point = Coord;

    #[inline]
    fn contains(&self, p: Coord) -> bool {
        ConfidenceEllipse::contains_coord(*self, p)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.semi_major() <= 0.0
    }
}

// ---- R3: the axis-aligned bounding box ----

impl PointRegion for Aabb {
    type Point = Vector;

    #[inline]
    fn contains(&self, p: Vector) -> bool {
        Aabb::contains(*self, p)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        // A normalized `Aabb` (`min ≤ max` per axis) is never empty.
        false
    }
}

impl Prunable for Aabb {
    #[inline]
    fn intersects(&self, other: &Self) -> bool {
        Aabb::overlaps(*self, *other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn rect_region_and_overlap() {
        let a = Rect::new(Coord::new(0.0, 0.0), Coord::new(4.0, 4.0));
        let b = Rect::new(Coord::new(3.0, 3.0), Coord::new(8.0, 8.0));
        let c = Rect::new(Coord::new(5.0, 5.0), Coord::new(9.0, 9.0));
        assert!(a.contains(Coord::new(2.0, 2.0)));
        assert!(!a.is_empty());
        assert!(a.intersects(&b)); // overlaps_of_mem
        assert!(disjoint(&a, &c));
    }

    #[test]
    fn ellipse_regions_membership_only() {
        let e = Ellipse::new(Coord::new(0.0, 0.0), 3.0, 2.0);
        assert!(e.contains(Coord::new(0.0, 0.0)));
        assert!(!e.contains(Coord::new(4.0, 0.0)));
        assert!(!e.is_empty());
        assert!(Ellipse::new(Coord::new(0.0, 0.0), 0.0, 2.0).is_empty());
        let ce = ConfidenceEllipse::from_covariance(Coord::new(0.0, 0.0), 4.0, 0.0, 1.0, 1.0);
        assert!(ce.contains(Coord::new(0.0, 0.0)));
        assert!(!ce.is_empty());
    }

    #[test]
    fn aabb_region_and_pruning() {
        use hausdorff_base::space::r3::{Aabb, Vector};
        let a = Aabb::new(Vector::new(0.0, 0.0, 0.0), Vector::new(4.0, 4.0, 4.0));
        let b = Aabb::new(Vector::new(3.0, 3.0, 3.0), Vector::new(8.0, 8.0, 8.0));
        let c = Aabb::new(Vector::new(5.0, 5.0, 5.0), Vector::new(9.0, 9.0, 9.0));
        assert!(a.contains(Vector::new(2.0, 2.0, 2.0)));
        assert!(!a.is_empty());
        assert!(!disjoint(&a, &b)); // overlaps_of_mem3
        assert!(disjoint(&a, &c)); // prune_sound3
        // The generic filter works in R3 too.
        let pts = [Vector::new(1.0, 1.0, 1.0), Vector::new(9.0, 9.0, 9.0)];
        assert_eq!(count_contained(&a, &pts), 1);
    }

    #[test]
    fn generic_filters_are_written_once() {
        let region = Rect::new(Coord::new(0.0, 0.0), Coord::new(10.0, 10.0));
        let pts = [
            Coord::new(1.0, 1.0),
            Coord::new(5.0, 5.0),
            Coord::new(11.0, 1.0), // outside
            Coord::new(-1.0, 5.0), // outside
        ];
        assert_eq!(count_contained(&region, &pts), 2);
        assert_eq!(filter_contained(&region, &pts).len(), 2);
        assert!(any_contained(&region, &pts));
        // The same generic functions work on a 1-D interval.
        let line = Interval::new(0.0, 1.0);
        assert_eq!(count_contained(&line, &[0.5, 2.0, 0.9, -3.0]), 2);
        assert!(!any_contained(&Interval::EMPTY, &[0.5]));
    }
}
