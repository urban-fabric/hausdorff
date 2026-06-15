//! Planar [`PointRegion`] / [`Prunable`] instances — the axis-aligned [`Rect`]
//! (the metric-ball family member with an exact pairwise overlap test) and the
//! elliptical disks [`Ellipse`] / [`ConfidenceEllipse`] (membership only, since
//! their pairwise intersection has no closed form).
//!
//! The traits and the generic spatial filters (`filter_contained`,
//! `count_contained`, `any_contained`, `disjoint`) live upstream in
//! [`hausdorff_geometry_traits`]; the line / circle / ℝ³ regions (`Interval`,
//! `Sector`, `Aabb`) implement them in `hausdorff-base`, and `Cap` in
//! `hausdorff-sphere`. Each membership/overlap law is machine-checked in
//! `riemann` (cited per impl).

use hausdorff_geometry_traits::{PointRegion, Prunable};

use crate::space::r2::{ConfidenceEllipse, Coord, Ellipse, Rect};

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
    /// Per-axis interval overlap — exact, and a `false` certifies disjointness
    /// (`riemann` `ExactArith/AABB.{overlaps_of_mem, prune_sound}`).
    #[inline]
    fn intersects(&self, other: &Self) -> bool {
        Rect::overlaps(*self, *other)
    }
}

// ---- R2: the elliptical disk and the confidence ellipse (membership only) ----

impl PointRegion for Ellipse {
    type Point = Coord;

    /// Division-free membership `b²(x−cx)² + a²(y−cy)² ≤ a²b²` (`riemann`
    /// `Hausdorff/Ellipse.inEllipse`; the region is convex, `inEllipse_convex`).
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

    /// The Mahalanobis sublevel set — a positive-semidefinite quadratic form, so
    /// convex (`riemann` `Hausdorff/Ellipse.inQuad_convex`).
    #[inline]
    fn contains(&self, p: Coord) -> bool {
        ConfidenceEllipse::contains_coord(*self, p)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.semi_major() <= 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_geometry_traits::{any_contained, count_contained, disjoint, filter_contained};

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
    }
}
