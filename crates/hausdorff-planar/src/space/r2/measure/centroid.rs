//! Planar [`Centroid`], mirroring `geo`'s `Centroid` trait. Polygons use the
//! area-weighted centroid (holes subtract), degenerate shapes fall back to the
//! mean of their vertices.

use crate::space::r2::util::ring_vertices;
use crate::space::r2::{
    Coord, Euclidean, Length, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point,
    Polygon, Rect, Triangle,
};

/// The centroid of a geometry, or `None` when it has no coordinates.
pub trait Centroid {
    fn centroid(&self) -> Option<Point>;
}

/// Area-weighted centroid of a single ring, returning `(centroid, signed_area)`.
/// `signed_area` is `0` for a degenerate ring, in which case the centroid is the
/// vertex mean.
///
/// Proof-backed by `riemann` `ExactArith/PolygonCentroid` (the shoelace-weighted
/// `cx`/`cy` accumulators collapse on a triangle to the vertex average over
/// `3·area` — `centroidX_triangle`/`centroidY_triangle`) and `ExactArith/Centroid`
/// (`centroid_in_triangle`: the area-weighted centroid lies inside the triangle).
fn ring_centroid(ring: &[Coord]) -> (Coord, f64) {
    let v = ring_vertices(ring);
    let n = v.len();
    if n == 0 {
        return (Coord::new(0.0, 0.0), 0.0);
    }
    if n < 3 {
        let mean = v.iter().fold(Coord::new(0.0, 0.0), |a, c| a + *c) / n as f64;
        return (mean, 0.0);
    }
    // Accumulated about `v[0]` and shifted back at the end. The sums are
    // translation-covariant in exact arithmetic, but not in `f64`: the raw form
    // builds `cx` from terms of size `|coord|³` to land on an answer of size
    // `|coord|·area`, so a small polygon far from the origin loses the centroid
    // entirely — at an offset of `1e10` a unit square's centroid came back
    // 3.3e9 away from the square. See
    // [`util::signed_area2`](crate::space::r2::util::signed_area2) for the same hazard in the
    // area alone.
    let o = v[0];
    let mut a2 = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let p = v[i];
        let q = v[(i + 1) % n];
        let (px, py) = (p.x - o.x, p.y - o.y);
        let (qx, qy) = (q.x - o.x, q.y - o.y);
        let cross = px * qy - qx * py;
        a2 += cross;
        cx += (px + qx) * cross;
        cy += (py + qy) * cross;
    }
    if a2 == 0.0 {
        let mean = v.iter().fold(Coord::new(0.0, 0.0), |a, c| a + *c) / n as f64;
        return (mean, 0.0);
    }
    (
        Coord::new(o.x + cx / (3.0 * a2), o.y + cy / (3.0 * a2)),
        a2 / 2.0,
    )
}

impl Centroid for Polygon {
    fn centroid(&self) -> Option<Point> {
        if self.exterior().0.is_empty() {
            return None;
        }
        let (c_ext, a_ext) = ring_centroid(&self.exterior().0);
        let mut area = a_ext.abs();
        let mut cx = c_ext.x * area;
        let mut cy = c_ext.y * area;
        for hole in self.interiors() {
            let (c_h, a_h) = ring_centroid(&hole.0);
            let w = a_h.abs();
            area -= w;
            cx -= c_h.x * w;
            cy -= c_h.y * w;
        }
        if area.abs() < f64::EPSILON {
            return Some(Point(c_ext));
        }
        Some(Point::new(cx / area, cy / area))
    }
}

impl Centroid for MultiPolygon {
    fn centroid(&self) -> Option<Point> {
        let mut area = 0.0;
        let mut cx = 0.0;
        let mut cy = 0.0;
        for poly in self.iter() {
            let (c, a) = ring_centroid(&poly.exterior().0);
            let mut w = a.abs();
            for hole in poly.interiors() {
                w -= (ring_centroid(&hole.0).1).abs();
            }
            if let Some(pc) = poly.centroid() {
                area += w;
                cx += pc.x() * w;
                cy += pc.y() * w;
            }
            let _ = c;
        }
        if area.abs() < f64::EPSILON {
            let pts: Vec<Point> = self.iter().filter_map(Centroid::centroid).collect();
            return MultiPoint::new(pts).centroid();
        }
        Some(Point::new(cx / area, cy / area))
    }
}

impl Centroid for LineString {
    fn centroid(&self) -> Option<Point> {
        match self.0.len() {
            0 => None,
            1 => Some(Point(self.0[0])),
            _ => {
                let mut len = 0.0;
                let mut cx = 0.0;
                let mut cy = 0.0;
                for w in self.0.windows(2) {
                    let d = w[1] - w[0];
                    let seg = (d.x * d.x + d.y * d.y).sqrt();
                    let mid = (w[0] + w[1]) / 2.0;
                    len += seg;
                    cx += mid.x * seg;
                    cy += mid.y * seg;
                }
                if len == 0.0 {
                    Some(Point(self.0[0]))
                } else {
                    Some(Point::new(cx / len, cy / len))
                }
            }
        }
    }
}

impl Centroid for Point {
    fn centroid(&self) -> Option<Point> {
        Some(*self)
    }
}

impl Centroid for MultiPoint {
    fn centroid(&self) -> Option<Point> {
        if self.is_empty() {
            return None;
        }
        let n = self.0.len() as f64;
        let sum = self.iter().fold(Coord::new(0.0, 0.0), |a, p| a + p.0);
        Some(Point(sum / n))
    }
}

impl Centroid for Line {
    fn centroid(&self) -> Option<Point> {
        Some(Point((self.start + self.end) / 2.0))
    }
}

impl Centroid for Rect {
    fn centroid(&self) -> Option<Point> {
        Some(self.center())
    }
}

impl Centroid for Triangle {
    fn centroid(&self) -> Option<Point> {
        Some(Point((self.0 + self.1 + self.2) / 3.0))
    }
}

impl Centroid for MultiLineString {
    /// Length-weighted mean of the member centroids; if every member is
    /// degenerate (zero length), the mean of the member centroids.
    fn centroid(&self) -> Option<Point> {
        let mut acc = Coord::new(0.0, 0.0);
        let mut total = 0.0;
        for ls in self.iter() {
            if let Some(c) = ls.centroid() {
                let len = Euclidean.length(ls);
                acc = acc + c.0 * len;
                total += len;
            }
        }
        if total > 0.0 {
            Some(Point(acc / total))
        } else {
            let pts: Vec<Point> = self.iter().filter_map(Centroid::centroid).collect();
            MultiPoint::new(pts).centroid()
        }
    }
}
