//! `r1` — one-dimensional Euclidean geometry: the [`Interval`] type, mirroring
//! golang/geo's `r1` package (whose whole content is `Interval`).
//!
//! `r1::Interval` is the real line's analogue of `r2::Rect` and the substrate of
//! the spherical [`s1`](crate::space::s1) interval: a closed range `[lo, hi]` on
//! ℝ, with the canonical *empty* interval `[1, 0]` (any `lo > hi`). It completes
//! the Euclidean dimensional ladder — the line `r1`, the plane `r2`, the space
//! `r3` — and carries no physical unit (metres enter only in [`crate::earth`]).

use hausdorff_geometry_traits::ClosedInterval;
use serde::{Deserialize, Serialize};

/// A closed interval `[lo, hi]` on the real line. `lo > hi` denotes the empty
/// interval; the canonical empty is [`Interval::EMPTY`] = `[1, 0]`.
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Interval {
    pub lo: f64,
    pub hi: f64,
}

/// Tolerance for [`Interval::approx_equal`], as in golang/geo's `r1`.
const EPSILON: f64 = 1e-15;

impl Interval {
    /// The canonical empty interval, `[1, 0]` (so any operation that produces
    /// `lo > hi` is empty, and `union`/`add_point` start from "no extent").
    pub const EMPTY: Self = Self { lo: 1.0, hi: 0.0 };

    /// The interval `[lo, hi]`.
    #[inline]
    #[must_use]
    pub const fn new(lo: f64, hi: f64) -> Self {
        Self { lo, hi }
    }

    /// The empty interval ([`Interval::EMPTY`]).
    #[inline]
    #[must_use]
    pub const fn empty() -> Self {
        Self::EMPTY
    }

    /// The zero-length interval `[p, p]` containing exactly the point `p`.
    #[inline]
    #[must_use]
    pub const fn from_point(p: f64) -> Self {
        Self { lo: p, hi: p }
    }

    /// True iff the interval is empty (`lo > hi`).
    #[inline]
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.lo > self.hi
    }

    /// The low endpoint.
    #[inline]
    #[must_use]
    pub const fn lo(self) -> f64 {
        self.lo
    }

    /// The high endpoint.
    #[inline]
    #[must_use]
    pub const fn hi(self) -> f64 {
        self.hi
    }

    /// The midpoint `(lo + hi) / 2`. Undefined (but finite-arithmetic stable) on
    /// the empty interval.
    #[inline]
    #[must_use]
    pub fn center(self) -> f64 {
        0.5 * (self.lo + self.hi)
    }

    /// The length `hi - lo` (negative for the empty interval).
    #[inline]
    #[must_use]
    pub fn length(self) -> f64 {
        self.hi - self.lo
    }

    /// True iff `p ∈ [lo, hi]` (endpoints included).
    #[inline]
    #[must_use]
    pub fn contains(self, p: f64) -> bool {
        self.lo <= p && p <= self.hi
    }

    /// True iff `p ∈ (lo, hi)` (endpoints excluded).
    #[inline]
    #[must_use]
    pub fn interior_contains(self, p: f64) -> bool {
        self.lo < p && p < self.hi
    }

    /// True iff `other` is wholly contained in `self`. The empty interval is
    /// contained in every interval.
    #[must_use]
    pub fn contains_interval(self, other: Self) -> bool {
        if other.is_empty() {
            return true;
        }
        self.lo <= other.lo && other.hi <= self.hi
    }

    /// True iff `other` is wholly contained in the **interior** of `self`. The
    /// empty interval is interior-contained in every interval.
    #[must_use]
    pub fn interior_contains_interval(self, other: Self) -> bool {
        if other.is_empty() {
            return true;
        }
        self.lo < other.lo && other.hi < self.hi
    }

    /// True iff `self` and `other` share at least one point.
    #[must_use]
    pub fn intersects(self, other: Self) -> bool {
        if self.lo <= other.lo {
            // other.lo ∈ self, and other is non-empty
            other.lo <= self.hi && other.lo <= other.hi
        } else {
            // self.lo ∈ other, and self is non-empty
            self.lo <= other.hi && self.lo <= self.hi
        }
    }

    /// True iff the **interior** of `self` shares a point with `other`.
    #[must_use]
    pub fn interior_intersects(self, other: Self) -> bool {
        other.lo < self.hi && self.lo < other.hi && self.lo < self.hi && other.lo <= other.hi
    }

    /// The intersection of `self` and `other` (empty if they are disjoint).
    #[must_use]
    pub fn intersection(self, other: Self) -> Self {
        Self {
            lo: self.lo.max(other.lo),
            hi: self.hi.min(other.hi),
        }
    }

    /// The smallest interval containing both `self` and `other` (the other if one
    /// is empty).
    #[must_use]
    pub fn union(self, other: Self) -> Self {
        if self.is_empty() {
            return other;
        }
        if other.is_empty() {
            return self;
        }
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }

    /// The smallest interval containing `self` and the point `p`.
    #[must_use]
    pub fn add_point(self, p: f64) -> Self {
        if self.is_empty() {
            return Self::from_point(p);
        }
        if p < self.lo {
            Self { lo: p, hi: self.hi }
        } else if p > self.hi {
            Self { lo: self.lo, hi: p }
        } else {
            self
        }
    }

    /// `self` grown (or, with a negative `margin`, shrunk) by `margin` on each
    /// side. The empty interval stays empty.
    #[must_use]
    pub fn expanded(self, margin: f64) -> Self {
        if self.is_empty() {
            return self;
        }
        Self {
            lo: self.lo - margin,
            hi: self.hi + margin,
        }
    }

    /// `p` clamped into `[lo, hi]` (the closest point of the interval to `p`).
    #[inline]
    #[must_use]
    pub fn clamp_point(self, p: f64) -> f64 {
        self.lo.max(self.hi.min(p))
    }

    /// True iff `self` and `other` are equal to within a fixed `1e-15` tolerance
    /// (golang/geo's `r1` epsilon); an empty interval is approximately equal to
    /// any interval of length `≤ 2·1e-15`.
    #[must_use]
    pub fn approx_equal(self, other: Self) -> bool {
        if self.is_empty() {
            return other.length() <= 2.0 * EPSILON;
        }
        if other.is_empty() {
            return self.length() <= 2.0 * EPSILON;
        }
        (other.lo - self.lo).abs() <= EPSILON && (other.hi - self.hi).abs() <= EPSILON
    }
}

/// `r1::Interval` is the canonical [`ClosedInterval`] — the real-line interval the
/// shared algebra is written over (delegating to the inherent `const fn`
/// methods). The interval laws are machine-checked in `riemann`
/// `Hausdorff/R1Interval`: `mem_inter_iff` (intersection = set intersection),
/// `mem_hull_of_mem_left/right` (union contains both), `mem_of_subset`
/// (containment is sound), `clamp_nearest` (the clamp is the nearest point).
impl ClosedInterval for Interval {
    type Scalar = f64;

    #[inline]
    fn from_bounds(lo: f64, hi: f64) -> Self {
        Self::new(lo, hi)
    }
    #[inline]
    fn lo(self) -> f64 {
        self.lo
    }
    #[inline]
    fn hi(self) -> f64 {
        self.hi
    }
    #[inline]
    fn empty() -> Self {
        Interval::empty()
    }
    #[inline]
    fn is_empty(self) -> bool {
        Interval::is_empty(self)
    }
    #[inline]
    fn contains(self, p: f64) -> bool {
        Interval::contains(self, p)
    }
    #[inline]
    fn contains_interval(self, other: Self) -> bool {
        Interval::contains_interval(self, other)
    }
    #[inline]
    fn intersection(self, other: Self) -> Self {
        Interval::intersection(self, other)
    }
    #[inline]
    fn union(self, other: Self) -> Self {
        Interval::union(self, other)
    }
    #[inline]
    fn expanded(self, margin: f64) -> Self {
        Interval::expanded(self, margin)
    }
}
