//! `s2::Rect` — a latitude/longitude rectangle on the sphere, ported from the
//! essentials of golang/geo's `s2.Rect` (with its `r1.Interval` latitude and
//! `s1.Interval` longitude internals). The longitude interval is wrap-aware: an
//! *inverted* interval (`lo > hi`) crosses the antimeridian, so a rectangle
//! like `lng ∈ [170°, −170°]` is first-class. Implements [`Region`], so a
//! [`super::RegionCoverer`] can cover it with cells; the cell tests run on the
//! cell's conservative lat/lng bound (golang/geo's `Cell.RectBound`), which
//! over-covers near cell corners but never misses a true intersection.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

use crate::earth::GeoPoint;
use crate::space::r3::Vector;
use crate::space::s1::Angle;

use super::coords;
use super::region::Region;
use super::{Cap, S2CellId, S2Point};

const DBL_EPSILON: f64 = f64::EPSILON;

/// A closed latitude interval in radians; `lo > hi` is the empty interval.
#[derive(Clone, Copy, Debug, PartialEq)]
struct LatInterval {
    lo: f64,
    hi: f64,
}

impl LatInterval {
    const fn new(lo: f64, hi: f64) -> Self {
        Self { lo, hi }
    }

    const fn empty() -> Self {
        Self { lo: 1.0, hi: 0.0 }
    }

    fn from_point_pair(a: f64, b: f64) -> Self {
        Self {
            lo: a.min(b),
            hi: a.max(b),
        }
    }

    fn is_empty(self) -> bool {
        self.lo > self.hi
    }

    fn contains(self, p: f64) -> bool {
        self.lo <= p && p <= self.hi
    }

    fn contains_interval(self, o: Self) -> bool {
        o.is_empty() || (self.lo <= o.lo && o.hi <= self.hi)
    }

    fn union(self, o: Self) -> Self {
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

    fn intersection(self, o: Self) -> Self {
        Self {
            lo: self.lo.max(o.lo),
            hi: self.hi.min(o.hi),
        }
    }

    fn expanded(self, margin: f64) -> Self {
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

/// A longitude interval in radians on the circle `[-π, π]`, after golang/geo's
/// `s1.Interval`: `lo > hi` is *inverted* (wraps through ±π), the empty
/// interval is `[π, −π]`, and the full circle is `[−π, π]`.
#[derive(Clone, Copy, Debug, PartialEq)]
struct LngInterval {
    lo: f64,
    hi: f64,
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
    fn from_endpoints(lo: f64, hi: f64) -> Self {
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
    fn from_point_pair(mut a: f64, mut b: f64) -> Self {
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

    const fn empty() -> Self {
        Self { lo: PI, hi: -PI }
    }

    const fn full() -> Self {
        Self { lo: -PI, hi: PI }
    }

    fn is_empty(self) -> bool {
        self.lo == PI && self.hi == -PI
    }

    fn is_full(self) -> bool {
        self.lo == -PI && self.hi == PI
    }

    fn is_inverted(self) -> bool {
        self.lo > self.hi
    }

    fn length(self) -> f64 {
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

    fn contains(self, p: f64) -> bool {
        let p = if p == -PI { PI } else { p };
        self.fast_contains(p)
    }

    fn contains_interval(self, o: Self) -> bool {
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

    fn union(self, o: Self) -> Self {
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

    fn intersection(self, o: Self) -> Self {
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

    fn expanded(self, margin: f64) -> Self {
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
    fn center(self) -> f64 {
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

/// A latitude/longitude rectangle: the product of a latitude interval and a
/// wrap-aware longitude interval, both stored in radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rect {
    lat: LatInterval,
    lng: LngInterval,
}

impl Rect {
    /// The empty rectangle (contains nothing).
    #[must_use]
    pub fn empty() -> Self {
        Self {
            lat: LatInterval::empty(),
            lng: LngInterval::empty(),
        }
    }

    /// The full rectangle (the whole sphere).
    #[must_use]
    pub fn full() -> Self {
        Self {
            lat: LatInterval::new(-FRAC_PI_2, FRAC_PI_2),
            lng: LngInterval::full(),
        }
    }

    /// A rectangle from corner latitudes/longitudes in **degrees**. A
    /// `lng_lo > lng_hi` pair wraps across the antimeridian (e.g.
    /// `from_degrees(10.0, 170.0, 20.0, -170.0)` spans 20° of longitude
    /// centered on ±180°).
    #[must_use]
    pub fn from_degrees(lat_lo: f64, lng_lo: f64, lat_hi: f64, lng_hi: f64) -> Self {
        Self {
            lat: LatInterval::new(lat_lo.to_radians(), lat_hi.to_radians()),
            lng: LngInterval::from_endpoints(lng_lo.to_radians(), lng_hi.to_radians()),
        }
    }

    /// True iff the rectangle contains no points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.lat.is_empty() || self.lng.is_empty()
    }

    /// True iff the rectangle is the whole sphere.
    #[must_use]
    pub fn is_full(&self) -> bool {
        self.lat == LatInterval::new(-FRAC_PI_2, FRAC_PI_2) && self.lng.is_full()
    }

    /// True iff the bounds are well formed: latitudes within `[-90°, 90°]`,
    /// longitudes within `[-180°, 180°]`, and the two intervals empty only
    /// together.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.lat.lo.abs() <= FRAC_PI_2
            && self.lat.hi.abs() <= FRAC_PI_2
            && self.lng.lo.abs() <= PI
            && self.lng.hi.abs() <= PI
            && self.lat.is_empty() == self.lng.is_empty()
    }

    /// The low corner `(lat_lo, lng_lo)` as a WGS84 coordinate.
    #[must_use]
    pub fn lo(&self) -> GeoPoint {
        GeoPoint::from_lat_lon(self.lat.lo.to_degrees(), self.lng.lo.to_degrees())
    }

    /// The high corner `(lat_hi, lng_hi)` as a WGS84 coordinate.
    #[must_use]
    pub fn hi(&self) -> GeoPoint {
        GeoPoint::from_lat_lon(self.lat.hi.to_degrees(), self.lng.hi.to_degrees())
    }

    /// True iff the rectangle contains the WGS84 coordinate.
    #[must_use]
    pub fn contains_lat_lng(&self, p: GeoPoint) -> bool {
        self.lat.contains(p.lat.to_radians()) && self.lng.contains(p.lon.to_radians())
    }

    /// True iff this rectangle fully contains `other`.
    #[must_use]
    pub fn contains_rect(&self, other: &Rect) -> bool {
        if other.is_empty() {
            return true;
        }
        self.lat.contains_interval(other.lat) && self.lng.contains_interval(other.lng)
    }

    /// The smallest rectangle containing both rectangles (longitudes take the
    /// shorter way around when the inputs are disjoint).
    #[must_use]
    pub fn union(&self, other: &Rect) -> Rect {
        Rect {
            lat: self.lat.union(other.lat),
            lng: self.lng.union(other.lng),
        }
    }

    /// The intersection of the two rectangles ([`Rect::empty`] if they are
    /// disjoint in either coordinate).
    #[must_use]
    pub fn intersection(&self, other: &Rect) -> Rect {
        let lat = self.lat.intersection(other.lat);
        let lng = self.lng.intersection(other.lng);
        if lat.is_empty() || lng.is_empty() {
            return Rect::empty();
        }
        Rect { lat, lng }
    }

    /// The rectangle grown by the given margins (radians), latitude clamped to
    /// the poles.
    fn expanded(&self, lat_margin: f64, lng_margin: f64) -> Rect {
        let lat = self
            .lat
            .expanded(lat_margin)
            .intersection(LatInterval::new(-FRAC_PI_2, FRAC_PI_2));
        if lat.is_empty() {
            return Rect::empty();
        }
        Rect {
            lat,
            lng: self.lng.expanded(lng_margin),
        }
    }

    /// If the rectangle reaches a pole, widen the longitude span to the full
    /// circle (every longitude passes through a pole).
    fn polar_closure(&self) -> Rect {
        if self.lat.lo == -FRAC_PI_2 || self.lat.hi == FRAC_PI_2 {
            Rect {
                lat: self.lat,
                lng: LngInterval::full(),
            }
        } else {
            *self
        }
    }

    /// The four corner points (lo/lo, lo/hi, hi/hi, hi/lo) as unit vectors.
    fn corner_points(&self) -> [S2Point; 4] {
        let p = |lat: f64, lng: f64| {
            S2Point::from_geo_point(GeoPoint::from_lat_lon(lat.to_degrees(), lng.to_degrees()))
        };
        [
            p(self.lat.lo, self.lng.lo),
            p(self.lat.lo, self.lng.hi),
            p(self.lat.hi, self.lng.hi),
            p(self.lat.hi, self.lng.lo),
        ]
    }

    /// A conservative lat/lng bound of `cell` — a `Rect` superset of the cell
    /// (golang/geo's `Cell.RectBound`). For level ≥ 1 cells the latitude and
    /// longitude extremes are attained at two diagonally opposite vertex pairs
    /// (the cell lies within one quadrant of its face); face cells use fixed
    /// per-face bounds. Both cases are expanded by a small epsilon and closed
    /// over the poles, so the bound never under-covers.
    fn cell_bound(cell: S2CellId) -> Rect {
        if cell.level() == 0 {
            let pole_min_lat = (1.0f64 / 3.0).sqrt().asin() - 0.5 * DBL_EPSILON;
            let bound = match cell.face() {
                0 => Rect {
                    lat: LatInterval::new(-FRAC_PI_4, FRAC_PI_4),
                    lng: LngInterval::from_endpoints(-FRAC_PI_4, FRAC_PI_4),
                },
                1 => Rect {
                    lat: LatInterval::new(-FRAC_PI_4, FRAC_PI_4),
                    lng: LngInterval::from_endpoints(FRAC_PI_4, 3.0 * FRAC_PI_4),
                },
                2 => Rect {
                    lat: LatInterval::new(pole_min_lat, FRAC_PI_2),
                    lng: LngInterval::full(),
                },
                3 => Rect {
                    lat: LatInterval::new(-FRAC_PI_4, FRAC_PI_4),
                    lng: LngInterval::from_endpoints(3.0 * FRAC_PI_4, -3.0 * FRAC_PI_4),
                },
                4 => Rect {
                    lat: LatInterval::new(-FRAC_PI_4, FRAC_PI_4),
                    lng: LngInterval::from_endpoints(-3.0 * FRAC_PI_4, -FRAC_PI_4),
                },
                _ => Rect {
                    lat: LatInterval::new(-FRAC_PI_2, -pole_min_lat),
                    lng: LngInterval::full(),
                },
            };
            return bound.expanded(DBL_EPSILON, 0.0).polar_closure();
        }

        let (face, (u_lo, u_hi), (v_lo, v_hi)) = cell.face_uv_bounds();
        let u = [u_lo, u_hi];
        let v = [v_lo, v_hi];
        // The face's u/v tangent axes point along z only on faces 3/4 (u) and
        // 0/1 (v); that decides which diagonal carries the latitude extremes.
        let u_axis_z = matches!(face, 3 | 4);
        let v_axis_z = matches!(face, 0 | 1);
        let i = usize::from(if u_axis_z {
            u_lo + u_hi > 0.0
        } else {
            u_lo + u_hi < 0.0
        });
        let j = usize::from(if v_axis_z {
            v_lo + v_hi > 0.0
        } else {
            v_lo + v_hi < 0.0
        });
        let lat = LatInterval::from_point_pair(
            corner_lat(face, u[i], v[j]),
            corner_lat(face, u[1 - i], v[1 - j]),
        );
        let lng = LngInterval::from_point_pair(
            corner_lng(face, u[i], v[1 - j]),
            corner_lng(face, u[1 - i], v[j]),
        );
        Rect { lat, lng }
            .expanded(2.0 * DBL_EPSILON, 2.0 * DBL_EPSILON)
            .polar_closure()
    }
}

/// Latitude (radians) of the face corner at `(u, v)`.
fn corner_lat(face: u8, u: f64, v: f64) -> f64 {
    let p = coords::face_uv_to_xyz(face, u, v);
    p.z.atan2(p.x.hypot(p.y))
}

/// Longitude (radians) of the face corner at `(u, v)`.
fn corner_lng(face: u8, u: f64, v: f64) -> f64 {
    let p = coords::face_uv_to_xyz(face, u, v);
    p.y.atan2(p.x)
}

impl Region for Rect {
    fn contains_point(&self, p: S2Point) -> bool {
        self.contains_lat_lng(p.to_geo_point())
    }

    /// Port of golang/geo's `Rect.CapBound`: the smaller of a polar cap (around
    /// the nearer pole, reaching the far latitude edge) and — when the
    /// longitude span is at most a half turn — a cap centered on the rectangle
    /// midpoint that reaches all four corners.
    fn bounding_cap(&self) -> Cap {
        if self.is_empty() {
            return Cap::empty(S2Point::from_vector(Vector::new(1.0, 0.0, 0.0)));
        }
        // Padding absorbs the cos() round trip in Cap::contains_point so a
        // boundary point computed at exactly the cap radius stays inside.
        let pad = 1e-14;
        let (pole_z, pole_angle) = if self.lat.lo + self.lat.hi < 0.0 {
            (-1.0, FRAC_PI_2 + self.lat.hi)
        } else {
            (1.0, FRAC_PI_2 - self.lat.lo)
        };
        let pole_cap = Cap::from_center_angle(
            S2Point::from_vector(Vector::new(0.0, 0.0, pole_z)),
            Angle::from_radians(pole_angle + pad),
        );
        if self.lng.length() <= PI {
            let mid_lat = 0.5 * (self.lat.lo + self.lat.hi);
            let center = S2Point::from_geo_point(GeoPoint::from_lat_lon(
                mid_lat.to_degrees(),
                self.lng.center().to_degrees(),
            ));
            let radius = self
                .corner_points()
                .iter()
                .map(|c| center.vector().angle(c.vector()).radians())
                .fold(0.0_f64, f64::max);
            let mid_cap = Cap::from_center_angle(center, Angle::from_radians(radius + pad));
            if mid_cap.height() < pole_cap.height() {
                return mid_cap;
            }
        }
        pole_cap
    }

    fn intersects_cell(&self, cell: S2CellId) -> bool {
        if self.is_empty() {
            return false;
        }
        !self.intersection(&Rect::cell_bound(cell)).is_empty()
    }

    fn contains_cell(&self, cell: S2CellId) -> bool {
        self.contains_rect(&Rect::cell_bound(cell))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::s2::RegionCoverer;

    fn gp(lat: f64, lon: f64) -> GeoPoint {
        GeoPoint::from_lat_lon(lat, lon)
    }

    #[test]
    fn from_degrees_contains_and_excludes() {
        let r = Rect::from_degrees(40.0, -5.0, 45.0, 5.0);
        assert!(r.is_valid() && !r.is_empty() && !r.is_full());
        assert!(r.contains_lat_lng(gp(42.0, 0.0)));
        assert!(r.contains_lat_lng(gp(40.0, -5.0)), "boundary is closed");
        assert!(!r.contains_lat_lng(gp(39.9, 0.0)));
        assert!(!r.contains_lat_lng(gp(42.0, 5.1)));
        assert!((r.lo().lat - 40.0).abs() < 1e-12 && (r.hi().lon - 5.0).abs() < 1e-12);
    }

    #[test]
    fn empty_and_full() {
        assert!(Rect::empty().is_empty());
        assert!(Rect::empty().is_valid());
        assert!(!Rect::empty().contains_lat_lng(gp(0.0, 0.0)));
        assert!(Rect::full().is_full());
        for &(lat, lon) in &[(0.0, 0.0), (90.0, 0.0), (-90.0, 0.0), (15.0, 180.0)] {
            assert!(Rect::full().contains_lat_lng(gp(lat, lon)));
        }
        assert!(Rect::full().contains_rect(&Rect::from_degrees(1.0, 2.0, 3.0, 4.0)));
    }

    #[test]
    fn antimeridian_wrap_contains_and_set_ops() {
        let r = Rect::from_degrees(10.0, 170.0, 20.0, -170.0);
        assert!(r.is_valid());
        assert!(r.contains_lat_lng(gp(15.0, 175.0)));
        assert!(r.contains_lat_lng(gp(15.0, -175.0)));
        assert!(r.contains_lat_lng(gp(15.0, 180.0)));
        assert!(r.contains_lat_lng(gp(15.0, -180.0)));
        assert!(!r.contains_lat_lng(gp(15.0, 0.0)));
        assert!(!r.contains_lat_lng(gp(15.0, -160.0)));
        assert!(!r.contains_lat_lng(gp(25.0, 175.0)));

        let other = Rect::from_degrees(12.0, -175.0, 18.0, -150.0);
        let inter = r.intersection(&other);
        assert!(inter.contains_lat_lng(gp(15.0, -172.0)));
        assert!(!inter.contains_lat_lng(gp(15.0, -160.0)));
        assert!(!inter.contains_lat_lng(gp(15.0, 175.0)));

        let union = r.union(&other);
        assert!(union.contains_lat_lng(gp(15.0, 175.0)));
        assert!(union.contains_lat_lng(gp(15.0, -160.0)));
        assert!(!union.contains_lat_lng(gp(15.0, 0.0)));
        assert!(union.contains_rect(&r) && union.contains_rect(&other));

        let disjoint = Rect::from_degrees(12.0, -10.0, 18.0, 10.0);
        assert!(r.intersection(&disjoint).is_empty());
    }

    #[test]
    fn union_intersection_of_plain_rects() {
        let a = Rect::from_degrees(0.0, 0.0, 10.0, 10.0);
        let b = Rect::from_degrees(5.0, 5.0, 15.0, 15.0);
        let u = a.union(&b);
        let i = a.intersection(&b);
        assert!(u.contains_rect(&a) && u.contains_rect(&b) && u.contains_rect(&i));
        assert!(a.contains_rect(&i) && b.contains_rect(&i));
        assert!(i.contains_lat_lng(gp(7.0, 7.0)));
        assert!(!i.contains_lat_lng(gp(3.0, 3.0)));
        assert!(u.contains_lat_lng(gp(3.0, 3.0)) && u.contains_lat_lng(gp(13.0, 13.0)));
        assert_eq!(a.union(&Rect::empty()), a);
        assert!(a.intersection(&Rect::empty()).is_empty());
    }

    #[test]
    fn cell_bound_contains_cell_samples() {
        for &(lat, lon, level) in &[
            (0.0, 0.0, 3),
            (37.7749, -122.4194, 10),
            (-33.8688, 151.2093, 6),
            (84.0, 5.0, 4),
            (-84.0, 5.0, 4),
            (0.1, 179.9, 8),
            (45.0, 45.0, 1),
        ] {
            let cell = S2CellId::from_geo_point_at_level(gp(lat, lon), level);
            let bound = Rect::cell_bound(cell);
            let mut samples: Vec<S2Point> = vec![cell.center_point()];
            let verts = cell.vertices();
            samples.extend_from_slice(&verts);
            for k in 0..4 {
                let mid = verts[k]
                    .vector()
                    .add(verts[(k + 1) % 4].vector())
                    .normalize();
                samples.push(S2Point::from_vector(mid));
            }
            for s in samples {
                assert!(
                    bound.contains_lat_lng(s.to_geo_point()),
                    "cell {cell} sample escaped its bound"
                );
            }
        }
        for face in 0..6u8 {
            let cell = S2CellId::from_face(face);
            let bound = Rect::cell_bound(cell);
            assert!(bound.contains_lat_lng(cell.center()));
            for v in cell.vertices() {
                assert!(bound.contains_lat_lng(v.to_geo_point()));
            }
        }
    }

    #[test]
    fn region_impl_covers_rect_and_respects_strict_contains() {
        let r = Rect::from_degrees(40.0, -5.0, 45.0, 5.0);
        let cover = RegionCoverer::at_level(7).covering(&r);
        assert!(!cover.is_empty());
        for &c in &cover {
            assert!(r.intersects_cell(c));
        }
        for &(lat, lon) in &[(40.5, -4.5), (42.5, 0.0), (44.9, 4.9)] {
            let leaf = S2CellId::from_geo_point(gp(lat, lon));
            assert!(
                cover.iter().any(|c| c.contains(leaf)),
                "interior point ({lat},{lon}) not covered"
            );
        }
        let outside = S2CellId::from_geo_point(gp(50.0, 20.0));
        assert!(!cover.iter().any(|c| c.contains(outside)));

        let deep_inside = S2CellId::from_geo_point_at_level(gp(42.5, 0.0), 10);
        assert!(r.contains_cell(deep_inside));
        assert!(r.intersects_cell(deep_inside));
        let far = S2CellId::from_geo_point_at_level(gp(-40.0, 120.0), 10);
        assert!(!r.intersects_cell(far));
        assert!(!r.contains_cell(far));
    }

    #[test]
    fn antimeridian_rect_adaptive_cover() {
        let r = Rect::from_degrees(10.0, 170.0, 20.0, -170.0);
        let cover = RegionCoverer {
            min_level: 0,
            max_level: 12,
            max_cells: 32,
        }
        .adaptive_covering(&r);
        assert!(!cover.is_empty());
        assert!(cover.len() <= 32);
        for &(lat, lon) in &[(15.0, 178.0), (15.0, -178.0), (11.0, 171.0), (19.0, -171.0)] {
            assert!(
                cover.contains_point(S2Point::from_geo_point(gp(lat, lon))),
                "wrapped point ({lat},{lon}) not covered"
            );
        }
    }

    #[test]
    fn bounding_cap_contains_rect_samples() {
        for r in [
            Rect::from_degrees(40.0, -5.0, 45.0, 5.0),
            Rect::from_degrees(10.0, 170.0, 20.0, -170.0),
            Rect::from_degrees(70.0, -180.0, 90.0, 180.0),
            Rect::from_degrees(-90.0, -180.0, -70.0, 180.0),
        ] {
            let cap = r.bounding_cap();
            for p in r.corner_points() {
                assert!(cap.contains_point(p), "corner escaped the bounding cap");
            }
            let mid = GeoPoint::from_lat_lon(
                0.5 * (r.lo().lat + r.hi().lat),
                r.lng.center().to_degrees(),
            );
            assert!(cap.contains_point(S2Point::from_geo_point(mid)));
        }
        assert!(Rect::empty().bounding_cap().is_empty());
    }
}
