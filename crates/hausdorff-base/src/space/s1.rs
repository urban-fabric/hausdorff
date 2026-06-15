//! `s1` — one-dimensional spherical geometry: angular measure. The home for
//! angles in the workspace.
//!
//! Angles live in hausdorff rather than `euclid` on purpose: unlike metres or
//! seconds they carry **no physical unit**, only a choice of representation
//! (radians or degrees) — so they are pure geometry. [`Angle`] is the canonical
//! radians-backed type the geometry computes in (a great-circle distance on the
//! unit sphere *is* an angle); [`AngleDegrees`] is the degrees representation for
//! I/O and human-facing values. Mirrors golang/geo's `s1.Angle` and subsumes
//! `euclid`'s `AngleRadians`/`AngleDegrees`.

#![allow(unsafe_code)]

use std::f64::consts::PI;
use std::ops::{Add, Div, Mul, Neg, Sub};

use bytemuck::{Pod, Zeroable};
use serde::{Deserialize, Serialize};

/// An angle in **radians** — the canonical angular measure for the geometry
/// code. POD layout + transparent serde, so it stores and serializes exactly
/// like the bare `f64` it wraps.
#[repr(transparent)]
#[derive(
    Copy, Clone, PartialEq, PartialOrd, Debug, Default, Pod, Zeroable, Serialize, Deserialize,
)]
#[serde(transparent)]
pub struct Angle(f64);

impl Angle {
    /// Zero radians.
    pub const ZERO: Self = Self(0.0);
    /// π radians (a straight angle).
    pub const PI: Self = Self(std::f64::consts::PI);
    /// 2π radians (a full turn).
    pub const TAU: Self = Self(std::f64::consts::TAU);
    /// A positive-infinite angle (a sentinel for "no/unbounded distance").
    pub const INFINITY: Self = Self(f64::INFINITY);

    /// An angle of `radians`.
    #[must_use]
    pub const fn from_radians(radians: f64) -> Self {
        Self(radians)
    }

    /// An angle of `degrees`.
    #[must_use]
    pub fn from_degrees(degrees: f64) -> Self {
        Self(degrees.to_radians())
    }

    /// The angle in radians.
    #[inline]
    #[must_use]
    pub const fn radians(self) -> f64 {
        self.0
    }

    /// The angle in degrees as a bare `f64`.
    #[inline]
    #[must_use]
    pub fn degrees(self) -> f64 {
        self.0.to_degrees()
    }

    /// The angle as a typed [`AngleDegrees`].
    #[must_use]
    pub fn to_degrees(self) -> AngleDegrees {
        AngleDegrees(self.0.to_degrees())
    }

    /// The radians value narrowed to `f32`.
    #[must_use]
    pub fn to_f32_radians(self) -> f32 {
        self.0 as f32
    }

    /// The absolute value.
    #[inline]
    #[must_use]
    pub fn abs(self) -> Self {
        Self(self.0.abs())
    }

    /// `sin` / `cos` / `tan` of the angle (radians need no extra conversion).
    #[inline]
    #[must_use]
    pub fn sin(self) -> f64 {
        self.0.sin()
    }
    #[inline]
    #[must_use]
    pub fn cos(self) -> f64 {
        self.0.cos()
    }
    #[inline]
    #[must_use]
    pub fn tan(self) -> f64 {
        self.0.tan()
    }

    /// Normalized into `(-π, π]`.
    #[must_use]
    pub fn normalized(self) -> Self {
        let mut r = self.0 % std::f64::consts::TAU;
        if r <= -std::f64::consts::PI {
            r += std::f64::consts::TAU;
        } else if r > std::f64::consts::PI {
            r -= std::f64::consts::TAU;
        }
        Self(r)
    }
}

impl Add for Angle {
    type Output = Angle;
    fn add(self, rhs: Angle) -> Angle {
        Angle(self.0 + rhs.0)
    }
}

impl Sub for Angle {
    type Output = Angle;
    fn sub(self, rhs: Angle) -> Angle {
        Angle(self.0 - rhs.0)
    }
}

impl Neg for Angle {
    type Output = Angle;
    fn neg(self) -> Angle {
        Angle(-self.0)
    }
}

impl Mul<f64> for Angle {
    type Output = Angle;
    fn mul(self, rhs: f64) -> Angle {
        Angle(self.0 * rhs)
    }
}

impl Div<f64> for Angle {
    type Output = Angle;
    fn div(self, rhs: f64) -> Angle {
        Angle(self.0 / rhs)
    }
}

/// An angle in **degrees** — the representation for I/O and human-facing values.
/// Convert to [`Angle`] (radians) before doing trigonometry.
#[repr(transparent)]
#[derive(
    Copy, Clone, PartialEq, PartialOrd, Debug, Default, Pod, Zeroable, Serialize, Deserialize,
)]
#[serde(transparent)]
pub struct AngleDegrees(f64);

impl AngleDegrees {
    /// Zero degrees.
    pub const ZERO: Self = Self(0.0);

    /// An angle of `degrees`.
    #[must_use]
    pub const fn new(degrees: f64) -> Self {
        Self(degrees)
    }

    /// The value in degrees.
    #[inline]
    #[must_use]
    pub const fn as_degrees(self) -> f64 {
        self.0
    }

    /// As a radians-backed [`Angle`].
    #[must_use]
    pub fn to_radians(self) -> Angle {
        Angle(self.0.to_radians())
    }
}

impl From<Angle> for AngleDegrees {
    fn from(a: Angle) -> Self {
        a.to_degrees()
    }
}

impl From<AngleDegrees> for Angle {
    fn from(d: AngleDegrees) -> Self {
        d.to_radians()
    }
}

/// An angular **sector** — a 1-D spherical cap on the unit circle: every angle
/// within an angular `radius` of a `center`. The S1 mirror of the spherical
/// `s2::Cap` (in `hausdorff-sphere`) and the circular analogue of
/// [`r1::Interval`](super::r1::Interval), it is the simplest angular region —
/// `O(1)` membership and pruning — useful for bearing/heading ranges,
/// field-of-view filters, and radial-sector queries.
///
/// A negative radius is the empty sector; a radius of `π` or more is the whole
/// circle. Membership ([`contains_angle`](Self::contains_angle)) is tested
/// **without an `acos`**, by the cosine of the angular gap — `cos(θ − center) ≥
/// cos(radius)` — exactly as `s2::Cap` tests its dot product.
///
/// Proof-backed by `riemann` `Hausdorff/Circular`:
/// - membership is the `acos`-free dot test (`inArc_iff_dot`: for `r ∈ [0, π]`,
///   `circDist θ α ≤ r ↔ cos r ≤ cos(θ − α)`);
/// - [`angular_distance`](Self::angular_distance) is the geodesic `circDist`
///   (`circDist`/`cos_circDist`), so the `clamp` before `acos` is a no-op;
/// - [`intersects`](Self::intersects)/[`contains_sector`](Self::contains_sector)
///   pruning is sound (`arcs_disjoint_of_far` + `circDist_triangle`, the angular
///   triangle inequality);
/// - [`expanded`](Self::expanded) only ever grows the region (`InArc.mono`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sector {
    center: Angle,
    radius: Angle,
}

impl Sector {
    /// The sector centred at `center` spanning `radius` to each side.
    #[must_use]
    pub fn from_center_radius(center: Angle, radius: Angle) -> Self {
        Self { center, radius }
    }

    /// The degenerate sector containing only the angle `a` (radius 0).
    #[must_use]
    pub fn from_angle(a: Angle) -> Self {
        Self {
            center: a,
            radius: Angle::ZERO,
        }
    }

    /// The empty sector centred at `center` (contains nothing).
    #[must_use]
    pub fn empty(center: Angle) -> Self {
        Self {
            center,
            radius: Angle::from_radians(-1.0),
        }
    }

    /// The full sector centred at `center` (the whole circle).
    #[must_use]
    pub fn full(center: Angle) -> Self {
        Self {
            center,
            radius: Angle::PI,
        }
    }

    /// The centre angle.
    #[inline]
    #[must_use]
    pub fn center(self) -> Angle {
        self.center
    }

    /// The angular radius (half-width).
    #[inline]
    #[must_use]
    pub fn radius(self) -> Angle {
        self.radius
    }

    /// True iff the sector contains no angles (negative radius).
    #[inline]
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.radius.radians() < 0.0
    }

    /// True iff the sector is the whole circle (radius `≥ π`).
    #[inline]
    #[must_use]
    pub fn is_full(self) -> bool {
        self.radius.radians() >= PI
    }

    /// The angular width `2·radius` — the arc the sector spans — clamped to a
    /// full turn `2π`; `0` for the empty sector.
    #[must_use]
    pub fn width(self) -> Angle {
        if self.is_empty() {
            Angle::ZERO
        } else {
            Angle::from_radians((2.0 * self.radius.radians()).min(Angle::TAU.radians()))
        }
    }

    /// Geodesic **angular distance** from the centre to `a` — the central angle
    /// `arccos(cos(a − center)) ∈ [0, π]` (the shorter of the two arcs). This is
    /// `riemann` `circDist`; `cos_circDist` proves the argument lands in
    /// `[-1, 1]`, so the `clamp` is a defensive no-op.
    #[must_use]
    pub fn angular_distance(self, a: Angle) -> Angle {
        let d = (a.radians() - self.center.radians()).cos().clamp(-1.0, 1.0);
        Angle::from_radians(d.acos())
    }

    /// True iff `a` lies in the sector — the `acos`-free test `cos(a − center) ≥
    /// cos(radius)` (`riemann` `inArc_iff_dot`).
    #[must_use]
    pub fn contains_angle(self, a: Angle) -> bool {
        if self.is_empty() {
            return false;
        }
        if self.is_full() {
            return true;
        }
        (a.radians() - self.center.radians()).cos() >= self.radius.radians().cos()
    }

    /// True iff this sector and `other` share any angle (`circDist(centres) ≤
    /// r₁ + r₂`). The contrapositive is `riemann` `arcs_disjoint_of_far`: centres
    /// farther apart than the sum of the radii share no angle, so a query prunes.
    #[must_use]
    pub fn intersects(self, other: Sector) -> bool {
        if self.is_empty() || other.is_empty() {
            return false;
        }
        if self.is_full() || other.is_full() {
            return true;
        }
        self.center_distance(other) <= self.radius.radians() + other.radius.radians()
    }

    /// True iff this sector fully contains `other` (`r₁ ≥ circDist(centres) +
    /// r₂`). Soundness — every angle of `other` then lies in `self` — is the
    /// angular triangle inequality (`riemann` `circDist_triangle`).
    #[must_use]
    pub fn contains_sector(self, other: Sector) -> bool {
        if self.is_full() || other.is_empty() {
            return true;
        }
        if self.is_empty() {
            return false;
        }
        self.radius.radians() >= self.center_distance(other) + other.radius.radians()
    }

    /// The **complementary** sector: the antipodal centre `center + π` and the
    /// complementary radius `π − radius`, covering exactly the angles this sector
    /// excludes.
    #[must_use]
    pub fn complement(self) -> Sector {
        Sector {
            center: Angle::from_radians(self.center.radians() + PI),
            radius: Angle::from_radians(PI - self.radius.radians()),
        }
    }

    /// The sector grown by `margin` (radius increased, clamped to a full circle);
    /// a negative `margin` shrinks it. Growing the radius only adds angles
    /// (`riemann` `InArc.mono`).
    #[must_use]
    pub fn expanded(self, margin: Angle) -> Sector {
        Sector {
            center: self.center,
            radius: Angle::from_radians((self.radius.radians() + margin.radians()).min(PI)),
        }
    }

    /// Angular distance between the two centres, in radians.
    fn center_distance(self, other: Sector) -> f64 {
        self.angular_distance(other.center).radians()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(d: f64) -> Angle {
        Angle::from_degrees(d)
    }

    #[test]
    fn membership_wraps_around_zero() {
        // A 30°-radius sector centred at 350°: contains 340..=20 (mod 360).
        let s = Sector::from_center_radius(deg(350.0), deg(30.0));
        assert!(s.contains_angle(deg(350.0))); // centre
        assert!(s.contains_angle(deg(10.0))); // 20° away, across 0°
        assert!(s.contains_angle(deg(340.0))); // 10° away
        assert!(s.contains_angle(deg(20.0))); // on the boundary
        assert!(!s.contains_angle(deg(40.0))); // 50° away — outside
        assert!(!s.contains_angle(deg(180.0))); // antipodal-ish — outside
    }

    #[test]
    fn empty_and_full() {
        let e = Sector::empty(deg(10.0));
        assert!(e.is_empty());
        assert!(!e.contains_angle(deg(10.0))); // even the "centre" is excluded
        let f = Sector::full(deg(10.0));
        assert!(f.is_full());
        assert!(f.contains_angle(deg(200.0)));
        // A point sector contains exactly its centre (to within rounding).
        let p = Sector::from_angle(deg(45.0));
        assert!(p.contains_angle(deg(45.0)));
        assert!(!p.contains_angle(deg(45.001)));
    }

    #[test]
    fn angular_distance_is_the_shorter_arc() {
        let s = Sector::from_center_radius(deg(10.0), deg(5.0));
        // 350° is 20° from 10° (the short way), not 340°.
        assert!((s.angular_distance(deg(350.0)).degrees() - 20.0).abs() < 1e-9);
        assert!(s.angular_distance(deg(10.0)).radians().abs() < 1e-12);
        // Always in [0, π].
        for d in [0.0, 90.0, 179.0, 181.0, 270.0, 359.0] {
            let dist = s.angular_distance(deg(d)).radians();
            assert!((0.0..=PI + 1e-9).contains(&dist));
        }
    }

    #[test]
    fn intersects_and_contains_with_pruning() {
        let a = Sector::from_center_radius(deg(0.0), deg(20.0));
        let b = Sector::from_center_radius(deg(30.0), deg(20.0)); // 30 ≤ 40 → overlap
        let c = Sector::from_center_radius(deg(90.0), deg(20.0)); // 90 > 40 → disjoint
        assert!(a.intersects(b));
        assert!(b.intersects(a));
        assert!(!a.intersects(c)); // arcs_disjoint_of_far prunes
        // Containment: a wide sector swallows a narrow nearby one.
        let wide = Sector::from_center_radius(deg(0.0), deg(60.0));
        let narrow = Sector::from_center_radius(deg(10.0), deg(20.0)); // 10 + 20 ≤ 60
        assert!(wide.contains_sector(narrow));
        assert!(!narrow.contains_sector(wide));
        // Containment ⟹ intersection (a genuine consistency check).
        assert!(wide.intersects(narrow));
    }

    #[test]
    fn contained_sector_points_are_all_inside() {
        // The triangle inequality (circDist_triangle) guarantee: if `wide`
        // contains `narrow`, every angle of `narrow` is in `wide`.
        let wide = Sector::from_center_radius(deg(100.0), deg(50.0));
        let narrow = Sector::from_center_radius(deg(120.0), deg(20.0));
        assert!(wide.contains_sector(narrow));
        for off in [-20.0, -10.0, 0.0, 10.0, 20.0] {
            assert!(wide.contains_angle(deg(120.0 + off)));
        }
    }

    #[test]
    fn expanded_and_complement() {
        let s = Sector::from_center_radius(deg(0.0), deg(20.0));
        let big = s.expanded(deg(30.0));
        assert!((big.radius().degrees() - 50.0).abs() < 1e-9);
        assert!(big.contains_angle(deg(45.0))); // grew to include it
        assert!(!s.contains_angle(deg(45.0)));
        // Expansion only adds angles (InArc.mono): everything in `s` stays in.
        for d in [-15.0, 0.0, 15.0] {
            assert!(big.contains_angle(deg(d)));
        }
        // The complement covers exactly what the sector excludes.
        let comp = s.complement();
        assert!((comp.radius().degrees() - 160.0).abs() < 1e-9);
        assert!(comp.contains_angle(deg(180.0)));
        assert!(!comp.contains_angle(deg(0.0)));
        // Expansion past π saturates to a full circle.
        assert!(s.expanded(Angle::PI).is_full());
    }

    #[test]
    fn width_tracks_radius() {
        let s = Sector::from_center_radius(deg(0.0), deg(20.0));
        assert!((s.width().degrees() - 40.0).abs() < 1e-9);
        assert_eq!(Sector::empty(deg(0.0)).width(), Angle::ZERO);
        // Full sector spans the whole circle.
        assert!((Sector::full(deg(0.0)).width().radians() - Angle::TAU.radians()).abs() < 1e-9);
    }
}
