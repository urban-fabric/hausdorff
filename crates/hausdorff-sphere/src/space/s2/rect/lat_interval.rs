//! `LatInterval` — the closed latitude interval (in radians) used as the
//! latitude coordinate of an `s2::Rect`. Ported from golang/geo's `r1.Interval`.

use hausdorff_geometry_traits::ClosedInterval;

/// A closed latitude interval in radians; `lo > hi` is the empty interval.
///
/// The `Rect` latitude coordinate, sharing the [`ClosedInterval`] algebra with
/// `r1::Interval`. Proof-backed by `riemann` `Hausdorff/SphereRect` (`LatInterval`):
/// `mem_inter_iff` (intersection = set intersection), `union_mem_left/right`
/// (union contains both, empty-absorbing), `subset_sound` (`contains_interval` is
/// sound, empty contained in everything), `expanded_mem_of_{mem,near}` (a
/// non-negative margin only grows it, capturing the margin band), and
/// `fromPointPair_minimal` (the two-point hull is least).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct LatInterval {
    pub(super) lo: f64,
    pub(super) hi: f64,
}

impl LatInterval {
    pub(super) const fn new(lo: f64, hi: f64) -> Self {
        Self { lo, hi }
    }

    pub(super) const fn empty() -> Self {
        Self { lo: 1.0, hi: 0.0 }
    }

    pub(super) fn from_point_pair(a: f64, b: f64) -> Self {
        Self {
            lo: a.min(b),
            hi: a.max(b),
        }
    }

    pub(super) fn is_empty(self) -> bool {
        self.lo > self.hi
    }

    pub(super) fn contains(self, p: f64) -> bool {
        self.lo <= p && p <= self.hi
    }

    pub(super) fn contains_interval(self, o: Self) -> bool {
        o.is_empty() || (self.lo <= o.lo && o.hi <= self.hi)
    }

    pub(super) fn union(self, o: Self) -> Self {
        if self.is_empty() {
            o
        } else if o.is_empty() {
            self
        } else {
            Self {
                lo: self.lo.min(o.lo),
                hi: self.hi.max(o.hi),
            }
        }
    }

    pub(super) fn intersection(self, o: Self) -> Self {
        Self {
            lo: self.lo.max(o.lo),
            hi: self.hi.min(o.hi),
        }
    }

    pub(super) fn expanded(self, margin: f64) -> Self {
        if self.is_empty() {
            self
        } else {
            Self {
                lo: self.lo - margin,
                hi: self.hi + margin,
            }
        }
    }
}

/// `LatInterval` shares the [`ClosedInterval`] algebra with `r1::Interval` (both
/// golang `r1.Interval` ports) — delegating to the inherent methods. Laws cited
/// on the struct (`riemann` `Hausdorff/SphereRect`).
impl ClosedInterval for LatInterval {
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
        LatInterval::empty()
    }
    #[inline]
    fn is_empty(self) -> bool {
        LatInterval::is_empty(self)
    }
    #[inline]
    fn contains(self, p: f64) -> bool {
        LatInterval::contains(self, p)
    }
    #[inline]
    fn contains_interval(self, other: Self) -> bool {
        LatInterval::contains_interval(self, other)
    }
    #[inline]
    fn intersection(self, other: Self) -> Self {
        LatInterval::intersection(self, other)
    }
    #[inline]
    fn union(self, other: Self) -> Self {
        LatInterval::union(self, other)
    }
    #[inline]
    fn expanded(self, margin: f64) -> Self {
        LatInterval::expanded(self, margin)
    }
    #[inline]
    fn from_point_pair(a: f64, b: f64) -> Self {
        LatInterval::from_point_pair(a, b)
    }
}
