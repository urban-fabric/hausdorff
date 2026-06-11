//! `r2::prepared` — a [`PreparedGeometry`]: index a geometry once, then answer
//! many spatial predicates against it cheaply.
//!
//! After GEOS's `PreparedGeometry` — the standard accelerator for spatial joins,
//! where one geometry is tested against a stream of others ("which rows fall in
//! this polygon", "distance from each row to this boundary"). The geometry's
//! boundary segments are bulk-loaded into an [`StrTree`] **once**; each query then
//! visits only nearby segments instead of all of them, so a batch amortizes the
//! index build (the amortization the per-call `euclidean_distance` cost model
//! can't capture).

use super::min_distance::{seg_box, segment_gap_sq, segments};
use super::{Coord, Geometry, Intersects, PointBatch, Rect, StrTree};

/// A geometry indexed for repeated spatial predicates against many others.
pub struct PreparedGeometry<'a> {
    geom: &'a Geometry,
    segments: Vec<(Coord, Coord)>,
    index: StrTree<usize>,
    bbox: Option<(Coord, Coord)>,
}

impl<'a> PreparedGeometry<'a> {
    /// Index `geom`'s boundary segments for fast repeated queries.
    #[must_use]
    pub fn new(geom: &'a Geometry) -> Self {
        let segs = segments(geom);
        let items: Vec<(Rect, usize)> = segs
            .iter()
            .enumerate()
            .map(|(i, &(p, q))| {
                let (lo, hi) = seg_box(p, q);
                (Rect::new(lo, hi), i)
            })
            .collect();
        let bbox = bbox_of(&segs);
        let index = StrTree::new(items);
        Self {
            geom,
            segments: segs,
            index,
            bbox,
        }
    }

    /// The geometry this index was built over.
    #[must_use]
    pub fn geometry(&self) -> &Geometry {
        self.geom
    }

    /// Minimum Euclidean distance to `other` — `0.0` if they intersect. Visits
    /// only the prepared segments near each of `other`'s segments.
    #[must_use]
    pub fn distance(&self, other: &Geometry) -> f64 {
        let other_segs = segments(other);
        if self.segments.is_empty() || other_segs.is_empty() {
            return f64::INFINITY;
        }
        if self.boxes_touch(&other_segs) && self.geom.intersects(other) {
            return 0.0;
        }
        let mut min_sq = f64::INFINITY;
        for &(a1, a2) in &other_segs {
            let (lo, hi) = seg_box(a1, a2);
            let d = min_sq.sqrt();
            let q = Rect::new(
                Coord::new(lo.x - d, lo.y - d),
                Coord::new(hi.x + d, hi.y + d),
            );
            self.index.query_each(q, |_, &j| {
                let (b1, b2) = self.segments[j];
                min_sq = min_sq.min(segment_gap_sq(a1, a2, b1, b2));
            });
        }
        min_sq.sqrt()
    }

    /// Whether `other` is within `d` (`ST_DWithin`), using the prepared index.
    #[must_use]
    pub fn within(&self, other: &Geometry, d: f64) -> bool {
        if d < 0.0 {
            return false;
        }
        let other_segs = segments(other);
        if self.segments.is_empty() || other_segs.is_empty() {
            return false;
        }
        if self.boxes_touch(&other_segs) && self.geom.intersects(other) {
            return true;
        }
        let d_sq = d * d;
        let mut found = false;
        for &(a1, a2) in &other_segs {
            if found {
                break;
            }
            let (lo, hi) = seg_box(a1, a2);
            let q = Rect::new(
                Coord::new(lo.x - d, lo.y - d),
                Coord::new(hi.x + d, hi.y + d),
            );
            self.index.query_each(q, |_, &j| {
                let (b1, b2) = self.segments[j];
                if segment_gap_sq(a1, a2, b1, b2) <= d_sq {
                    found = true;
                }
            });
        }
        found
    }

    /// Whether the prepared geometry intersects `other`. A bounding-box reject in
    /// front of the exact [`Intersects`] kernel.
    #[must_use]
    pub fn intersects(&self, other: &Geometry) -> bool {
        self.boxes_touch(&segments(other)) && self.geom.intersects(other)
    }

    /// Whether `point` lies inside the prepared **areal** geometry, by indexed
    /// crossing-number ray casting — only the prepared segments crossing a
    /// rightward ray at the point's latitude are visited. (For non-areal prepared
    /// geometries the result is not meaningful.)
    #[must_use]
    pub fn contains_point(&self, point: Coord) -> bool {
        let Some((lo, hi)) = self.bbox else {
            return false;
        };
        if point.x < lo.x || point.x > hi.x || point.y < lo.y || point.y > hi.y {
            return false;
        }
        let ray = Rect::new(point, Coord::new(f64::INFINITY, point.y));
        let mut crossings = 0u32;
        self.index.query_each(ray, |_, &j| {
            let (a, b) = self.segments[j];
            if ray_crosses(point, a, b) {
                crossings += 1;
            }
        });
        crossings % 2 == 1
    }

    /// For each point in `points`, whether it lies inside this prepared **areal**
    /// geometry — a vectorized [`contains_point`](Self::contains_point) over a
    /// columnar batch: the point-in-polygon join.
    #[must_use]
    pub fn contains_points(&self, points: &PointBatch) -> Vec<bool> {
        (0..points.len())
            .map(|i| self.contains_point(points.coord(i)))
            .collect()
    }

    /// For each geometry in `others`, whether it is within `d` of this prepared
    /// geometry (`ST_DWithin`) — a vectorized [`within`](Self::within), the
    /// distance-join predicate.
    #[must_use]
    pub fn within_each(&self, others: &[Geometry], d: f64) -> Vec<bool> {
        others.iter().map(|o| self.within(o, d)).collect()
    }

    /// For each geometry in `others`, whether it intersects this prepared
    /// geometry — a vectorized [`intersects`](Self::intersects).
    #[must_use]
    pub fn intersects_each(&self, others: &[Geometry]) -> Vec<bool> {
        others.iter().map(|o| self.intersects(o)).collect()
    }

    /// Minimum distance from this prepared geometry to each geometry in `others`
    /// — a vectorized [`distance`](Self::distance).
    #[must_use]
    pub fn distance_each(&self, others: &[Geometry]) -> Vec<f64> {
        others.iter().map(|o| self.distance(o)).collect()
    }

    /// Bounding-box pre-test: could the prepared geometry touch `other`'s segments?
    fn boxes_touch(&self, other_segs: &[(Coord, Coord)]) -> bool {
        match (self.bbox, bbox_of(other_segs)) {
            (Some(s), Some(o)) => box_gap_sq(s, o) == 0.0,
            _ => false,
        }
    }
}

/// Does a rightward horizontal ray from `p` cross segment `a`–`b`? Half-open in
/// `y` (`[min_y, max_y)`) so a vertex on the ray is counted once.
fn ray_crosses(p: Coord, a: Coord, b: Coord) -> bool {
    let (lo, hi) = if a.y <= b.y { (a, b) } else { (b, a) };
    if p.y < lo.y || p.y >= hi.y {
        return false;
    }
    let t = (p.y - lo.y) / (hi.y - lo.y);
    lo.x + t * (hi.x - lo.x) > p.x
}

fn bbox_of(segs: &[(Coord, Coord)]) -> Option<(Coord, Coord)> {
    if segs.is_empty() {
        return None;
    }
    let mut lo = Coord::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Coord::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(p, q) in segs {
        for c in [p, q] {
            lo = Coord::new(lo.x.min(c.x), lo.y.min(c.y));
            hi = Coord::new(hi.x.max(c.x), hi.y.max(c.y));
        }
    }
    Some((lo, hi))
}

fn box_gap_sq(a: (Coord, Coord), b: (Coord, Coord)) -> f64 {
    let dx = (b.0.x - a.1.x).max(a.0.x - b.1.x).max(0.0);
    let dy = (b.0.y - a.1.y).max(a.0.y - b.1.y).max(0.0);
    dx * dx + dy * dy
}
