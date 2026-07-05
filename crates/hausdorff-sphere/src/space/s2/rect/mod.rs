//! `s2::Rect` — a latitude/longitude rectangle on the sphere, ported from the
//! essentials of golang/geo's `s2.Rect` (with its `r1.Interval` latitude and
//! `s1.Interval` longitude internals). The longitude interval is wrap-aware: an
//! *inverted* interval (`lo > hi`) crosses the antimeridian, so a rectangle
//! like `lng ∈ [170°, −170°]` is first-class. Implements [`Region`], so a
//! [`super::RegionCoverer`] can cover it with cells; the cell tests run on the
//! cell's conservative lat/lng bound (golang/geo's `Cell.RectBound`), which
//! over-covers near cell corners but never misses a true intersection.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use hausdorff_base::earth::GeoPoint;

use super::coords;
use super::region::Region;
use super::{S2CellId, S2Point};

mod lat_interval;
mod lng_interval;
mod region_impl;

#[cfg(test)]
mod tests;

use lat_interval::LatInterval;
use lng_interval::LngInterval;

pub(super) const DBL_EPSILON: f64 = f64::EPSILON;

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

    /// The smallest rectangle containing a single WGS84 point (a degenerate
    /// rectangle whose corners coincide). [`Rect::empty`] grown by
    /// [`add_point`](Self::add_point) reaches this.
    #[must_use]
    pub fn from_point(p: GeoPoint) -> Self {
        Self {
            lat: LatInterval::from_point_pair(p.lat.to_radians(), p.lat.to_radians()),
            lng: LngInterval::from_point_pair(p.lon.to_radians(), p.lon.to_radians()),
        }
    }

    /// This rectangle grown to also contain `p`. Latitude grows by the closed
    /// interval; **longitude grows the shorter way around the circle** (golang/geo
    /// `s1.Interval.Union` of the point), so a point just west of the antimeridian
    /// added to one just east yields a narrow *inverted* span rather than the whole
    /// globe. Adding to [`Rect::empty`] yields [`from_point`](Self::from_point).
    #[must_use]
    pub fn add_point(self, p: GeoPoint) -> Self {
        self.union(&Rect::from_point(p))
    }

    /// The smallest rectangle containing every point of `points` — the **vertex
    /// bound**, folding each point through [`add_point`](Self::add_point). Empty
    /// for an empty slice.
    ///
    /// The longitude span is antimeridian-correct (each point joins the shorter
    /// way, so a ring straddling ±180° gets a tight *inverted* longitude interval,
    /// not the full circle). The latitude span is the vertex min/max — it does
    /// **not** add the poleward bulge a great-circle edge can reach between two
    /// vertices (golang/geo's edge-aware `RectBounder`). So this bounds a geometry
    /// by its vertices, matching the planar `bounding_rect` in latitude while fixing
    /// the longitude wrap; it is a conservative superset only for edges whose
    /// latitude extreme is at a vertex (all edges away from the poles). Callers
    /// needing the edge-exact bound near the poles should widen latitude to the
    /// pole.
    #[must_use]
    pub fn from_geo_points(points: &[GeoPoint]) -> Self {
        points
            .iter()
            .fold(Rect::empty(), |bound, &p| bound.add_point(p))
    }

    /// True iff the rectangle contains no points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.lat.is_empty() || self.lng.is_empty()
    }

    /// True iff the longitude span crosses the antimeridian — its low corner
    /// longitude is **east** of (greater than) its high corner longitude, the
    /// wraparound case a planar `xmin ≤ xmax` box cannot represent. `false` for the
    /// empty and full rectangles.
    #[must_use]
    pub fn is_lng_inverted(&self) -> bool {
        !self.is_empty() && !self.lng.is_full() && self.lng.lo > self.lng.hi
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
    ///
    /// Proof-backed by `riemann` `Hausdorff/SphereRect.contains_rect_sound`: a
    /// passing test implies every point of `other` is in `self` — the latitude
    /// half by closed-interval reasoning (`LatInterval.subset_sound`), the
    /// longitude half by the antimeridian-correct `S1Interval.subset_sound`. So
    /// rectangle containment on the sphere is sound in both coordinates.
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
