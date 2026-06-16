//! `s2::Cap` — a spherical cap (a disc on the unit sphere): every point within
//! an angular radius of a center. The simplest spherical region — O(1) area and
//! containment, exact at the poles and antimeridian since it is all unit-vector
//! arithmetic. Ported from golang/geo's `s2.Cap`.

use std::f64::consts::PI;

use hausdorff_base::space::s1::Angle;
use hausdorff_planar::algebra::{PointRegion, Prunable};

use super::S2Point;

/// A spherical cap: all points `p` with angular distance `center → p ≤ radius`.
/// A negative radius is the empty cap; a radius of `π` or more is the full
/// sphere.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cap {
    center: S2Point,
    radius: Angle,
}

impl Cap {
    /// The cap centered at `center` with angular `radius`.
    #[must_use]
    pub fn from_center_angle(center: S2Point, radius: Angle) -> Self {
        Self { center, radius }
    }

    /// The degenerate cap containing only `p` (radius 0).
    #[must_use]
    pub fn from_point(p: S2Point) -> Self {
        Self {
            center: p,
            radius: Angle::from_radians(0.0),
        }
    }

    /// The empty cap centered at `center` (contains nothing).
    #[must_use]
    pub fn empty(center: S2Point) -> Self {
        Self {
            center,
            radius: Angle::from_radians(-1.0),
        }
    }

    /// The full cap centered at `center` (the whole sphere).
    #[must_use]
    pub fn full(center: S2Point) -> Self {
        Self {
            center,
            radius: Angle::from_radians(PI),
        }
    }

    /// The center point.
    #[must_use]
    pub fn center(self) -> S2Point {
        self.center
    }

    /// The angular radius.
    #[must_use]
    pub fn radius(self) -> Angle {
        self.radius
    }

    /// True iff the cap contains no points.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.radius.radians() < 0.0
    }

    /// True iff the cap is the whole sphere.
    #[must_use]
    pub fn is_full(self) -> bool {
        self.radius.radians() >= PI
    }

    /// The cap height `1 − cos(radius)` (0 for empty, 2 for full).
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphereCap`: `capHeight ∈ [0, 2]`
    /// (`capHeight_mem`, so the `clamp` is a no-op on the true value), `0` at
    /// `r = 0` (`capHeight_zero`) and `2` at `r = π` (`capHeight_pi`), monotone in
    /// the radius (`capHeight_mono`).
    #[must_use]
    pub fn height(self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            (1.0 - self.radius.radians().cos()).clamp(0.0, 2.0)
        }
    }

    /// Surface area in **steradians**: `2π·height`, in `[0, 4π]`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphereCap`: `capArea = 2π(1 − cos r)`
    /// (Archimedes' hat-box) is non-negative (`capArea_nonneg`), `0` for the empty
    /// cap (`capArea_zero`) and `4π` for the full sphere (`capArea_pi`, the
    /// consistency check), and monotone in the radius (`capArea_mono`).
    #[must_use]
    pub fn area(self) -> f64 {
        2.0 * PI * self.height()
    }

    /// True iff `p` lies in the cap (`center · p ≥ cos(radius)`).
    #[must_use]
    pub fn contains_point(self, p: S2Point) -> bool {
        if self.is_empty() {
            return false;
        }
        self.center.vector().dot(p.vector()) >= self.radius.radians().cos()
    }

    /// The complementary cap: the antipodal center and the complementary radius
    /// `π − radius`, covering exactly the points this cap excludes.
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphereCap.inCap_complement_iff` (via
    /// `sphereDist_antipode`, distance to `-c` is `π −` distance to `c`): the
    /// `(-c, π − r)` cap contains exactly the points at great-circle distance
    /// `≥ r` from `c` — the closed complement of this cap (they share the
    /// boundary `dist = r`, as on S2).
    #[must_use]
    pub fn complement(self) -> Cap {
        let antipode = S2Point::from_vector(self.center.vector().mul(-1.0));
        Cap {
            center: antipode,
            radius: Angle::from_radians(PI - self.radius.radians()),
        }
    }

    /// True iff this cap fully contains `other`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphereCap`: the radius test is sound
    /// by `cap_contains_cap` (`sphereDist c₁ c₂ + r₂ ≤ r₁ ⟹ (c₂, r₂) ⊆ (c₁, r₁)`,
    /// the great-circle triangle inequality — the complement of the
    /// `caps_disjoint_of_far` pruning), and the `is_full` short-circuit by
    /// `cap_full_contains` (a radius-`≥ π` cap is the whole sphere).
    #[must_use]
    pub fn contains_cap(self, other: Cap) -> bool {
        if self.is_full() || other.is_empty() {
            return true;
        }
        self.radius.radians() >= self.center_distance(other) + other.radius.radians()
    }

    /// True iff this cap and `other` share any point.
    #[must_use]
    pub fn intersects(self, other: Cap) -> bool {
        if self.is_empty() || other.is_empty() {
            return false;
        }
        self.center_distance(other) <= self.radius.radians() + other.radius.radians()
    }

    /// The cap grown by `distance` (radius increased, clamped to a full sphere).
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphereCap`: expanding only grows the
    /// cap — every point of the original stays in it (`cap_expanded_contains`) —
    /// and it is a sound `distance`-buffer: any point within `distance` of the cap
    /// lands inside (`cap_expanded_margin`, the great-circle triangle inequality).
    #[must_use]
    pub fn expanded(self, distance: Angle) -> Cap {
        Cap {
            center: self.center,
            radius: Angle::from_radians((self.radius.radians() + distance.radians()).min(PI)),
        }
    }

    fn center_distance(self, other: Cap) -> f64 {
        self.center.vector().angle(other.center.vector()).radians()
    }
}

/// The spherical cap as a cross-space [`PointRegion`] — the S2 instance of the
/// bounding-volume abstraction (alongside `Interval`/`Sector`/`Rect`/`Ellipse`).
/// Membership is the `acos`-free dot test, proof-backed by `riemann`
/// `Hausdorff/SphereCap.inCap_iff_dot`.
impl PointRegion for Cap {
    type Point = S2Point;

    #[inline]
    fn contains(&self, p: S2Point) -> bool {
        Cap::contains_point(*self, p)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        Cap::is_empty(*self)
    }
}

/// Cap–cap pruning: far-apart centres share no point, so the pair prunes —
/// proof-backed by `riemann` `Hausdorff/SphereCap.caps_disjoint_of_far`.
impl Prunable for Cap {
    #[inline]
    fn intersects(&self, other: &Self) -> bool {
        Cap::intersects(*self, *other)
    }
}

#[cfg(test)]
mod region_tests {
    use super::*;
    use hausdorff_planar::algebra::{PointRegion, disjoint, filter_contained};

    fn pt(x: f64, y: f64, z: f64) -> S2Point {
        S2Point::from_vector(hausdorff_base::space::r3::Vector::new(x, y, z).normalize())
    }

    #[test]
    fn cap_is_a_point_region_with_pruning() {
        let north = pt(0.0, 0.0, 1.0);
        let cap = Cap::from_center_angle(north, Angle::from_degrees(30.0));
        // Membership (inCap_iff_dot): the pole is in, the equator is out.
        assert!(PointRegion::contains(&cap, north));
        assert!(!PointRegion::contains(&cap, pt(1.0, 0.0, 0.0)));
        assert!(!PointRegion::is_empty(&cap));
        assert!(PointRegion::is_empty(&Cap::empty(north)));
        // Pruning (caps_disjoint_of_far): a far cap is provably disjoint.
        let near = Cap::from_center_angle(pt(0.2, 0.0, 1.0), Angle::from_degrees(30.0));
        let far = Cap::from_center_angle(pt(1.0, 0.0, 0.0), Angle::from_degrees(30.0));
        assert!(!disjoint(&cap, &near));
        assert!(disjoint(&cap, &far));
        // The generic filter works on the sphere too.
        let pts = [north, pt(0.1, 0.1, 1.0), pt(1.0, 0.0, 0.0)];
        assert_eq!(filter_contained(&cap, &pts).len(), 2);
    }
}
