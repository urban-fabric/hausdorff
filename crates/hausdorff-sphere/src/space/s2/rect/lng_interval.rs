//! `LngInterval` — the wrap-aware longitude interval (in radians) used as the
//! longitude coordinate of an `s2::Rect`. Ported from golang/geo's `s1.Interval`:
//! an *inverted* interval (`lo > hi`) wraps through ±π.

use std::f64::consts::{PI, TAU};

use super::DBL_EPSILON;

/// A longitude interval in radians on the circle `[-π, π]`, after golang/geo's
/// `s1.Interval`: `lo > hi` is *inverted* (wraps through ±π), the empty
/// interval is `[π, −π]`, and the full circle is `[−π, π]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct LngInterval {
    pub(super) lo: f64,
    pub(super) hi: f64,
}

/// The distance from `a` to `b` travelling counter-clockwise (never negative).
fn positive_distance(a: f64, b: f64) -> f64 {
    let d = b - a;
    if d >= 0.0 { d } else { (b + PI) - (a - PI) }
}

/// IEEE-style remainder of `x` by `y` (result in `[-y/2, y/2]`).
fn remainder(x: f64, y: f64) -> f64 {
    let mut r = x % y;
    if r > 0.5 * y {
        r -= y;
    } else if r < -0.5 * y {
        r += y;
    }
    r
}

impl LngInterval {
    /// An interval from endpoints, mapping a `-π` endpoint to `π` so the
    /// representation is canonical (golang/geo's `IntervalFromEndpoints`).
    pub(super) fn from_endpoints(lo: f64, hi: f64) -> Self {
        let mut i = Self { lo, hi };
        if lo == -PI && hi != PI {
            i.lo = PI;
        }
        if hi == -PI && lo != PI {
            i.hi = PI;
        }
        i
    }

    /// The minimal interval containing both points (the shorter way around).
    pub(super) fn from_point_pair(mut a: f64, mut b: f64) -> Self {
        if a == -PI {
            a = PI;
        }
        if b == -PI {
            b = PI;
        }
        if positive_distance(a, b) <= PI {
            Self { lo: a, hi: b }
        } else {
            Self { lo: b, hi: a }
        }
    }

    pub(super) const fn empty() -> Self {
        Self { lo: PI, hi: -PI }
    }

    pub(super) const fn full() -> Self {
        Self { lo: -PI, hi: PI }
    }

    pub(super) fn is_empty(self) -> bool {
        self.lo == PI && self.hi == -PI
    }

    pub(super) fn is_full(self) -> bool {
        self.lo == -PI && self.hi == PI
    }

    fn is_inverted(self) -> bool {
        self.lo > self.hi
    }

    pub(super) fn length(self) -> f64 {
        let mut l = self.hi - self.lo;
        if l >= 0.0 {
            return l;
        }
        l += TAU;
        if l > 0.0 { l } else { -1.0 }
    }

    fn fast_contains(self, p: f64) -> bool {
        if self.is_inverted() {
            (p >= self.lo || p <= self.hi) && !self.is_empty()
        } else {
            p >= self.lo && p <= self.hi
        }
    }

    pub(super) fn contains(self, p: f64) -> bool {
        let p = if p == -PI { PI } else { p };
        self.fast_contains(p)
    }

    /// Whether this interval contains all of `o` (golang's `ContainsInterval`).
    ///
    /// Proof-backed by `riemann` `Hausdorff/S1Interval.subset_sound`: across all
    /// four inverted/normal combinations (and the full/empty edge cases) a
    /// passing test implies every longitude of `o` is in `self` — the
    /// antimeridian-correct soundness `Rect::contains_rect` rests on.
    pub(super) fn contains_interval(self, o: Self) -> bool {
        if self.is_inverted() {
            if o.is_inverted() {
                return o.lo >= self.lo && o.hi <= self.hi;
            }
            return (o.lo >= self.lo || o.hi <= self.hi) && !self.is_empty();
        }
        if o.is_inverted() {
            return self.is_full() || o.is_empty();
        }
        o.lo >= self.lo && o.hi <= self.hi
    }

    pub(super) fn union(self, o: Self) -> Self {
        if o.is_empty() {
            return self;
        }
        if self.fast_contains(o.lo) {
            if self.fast_contains(o.hi) {
                if self.contains_interval(o) {
                    return self;
                }
                return Self::full();
            }
            return Self {
                lo: self.lo,
                hi: o.hi,
            };
        }
        if self.fast_contains(o.hi) {
            return Self {
                lo: o.lo,
                hi: self.hi,
            };
        }
        if self.is_empty() || o.fast_contains(self.lo) {
            return o;
        }
        let dlo = positive_distance(o.hi, self.lo);
        let dhi = positive_distance(self.hi, o.lo);
        if dlo < dhi {
            Self {
                lo: o.lo,
                hi: self.hi,
            }
        } else {
            Self {
                lo: self.lo,
                hi: o.hi,
            }
        }
    }

    pub(super) fn intersection(self, o: Self) -> Self {
        if o.is_empty() {
            return Self::empty();
        }
        if self.fast_contains(o.lo) {
            if self.fast_contains(o.hi) {
                if o.length() < self.length() {
                    return o;
                }
                return self;
            }
            return Self {
                lo: o.lo,
                hi: self.hi,
            };
        }
        if self.fast_contains(o.hi) {
            return Self {
                lo: self.lo,
                hi: o.hi,
            };
        }
        if o.fast_contains(self.lo) {
            return self;
        }
        Self::empty()
    }

    pub(super) fn expanded(self, margin: f64) -> Self {
        if self.is_empty() {
            return self;
        }
        if self.length() + 2.0 * margin + 2.0 * DBL_EPSILON >= TAU {
            return Self::full();
        }
        let mut lo = remainder(self.lo - margin, TAU);
        let hi = remainder(self.hi + margin, TAU);
        if lo <= -PI {
            lo = PI;
        }
        Self { lo, hi }
    }

    /// The interval midpoint, on the wrapped circle for inverted intervals.
    pub(super) fn center(self) -> f64 {
        let c = 0.5 * (self.lo + self.hi);
        if !self.is_inverted() {
            c
        } else if c <= 0.0 {
            c + PI
        } else {
            c - PI
        }
    }
}
