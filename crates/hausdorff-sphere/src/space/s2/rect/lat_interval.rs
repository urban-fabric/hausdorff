//! `LatInterval` — the closed latitude interval (in radians) used as the
//! latitude coordinate of an `s2::Rect`. Ported from golang/geo's `r1.Interval`.

/// A closed latitude interval in radians; `lo > hi` is the empty interval.
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
