//! Planar [`Centroid`], mirroring `geo`'s `Centroid` trait. Polygons use the
//! area-weighted centroid (holes subtract), degenerate shapes fall back to the
//! mean of their vertices.

use super::util::ring_vertices;
use super::{Coord, LineString, MultiPoint, MultiPolygon, Point, Polygon};

/// The centroid of a geometry, or `None` when it has no coordinates.
pub trait Centroid {
    fn centroid(&self) -> Option<Point>;
}

/// Area-weighted centroid of a single ring, returning `(centroid, signed_area)`.
/// `signed_area` is `0` for a degenerate ring, in which case the centroid is the
/// vertex mean.
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
    let mut a2 = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let p = v[i];
        let q = v[(i + 1) % n];
        let cross = p.x * q.y - q.x * p.y;
        a2 += cross;
        cx += (p.x + q.x) * cross;
        cy += (p.y + q.y) * cross;
    }
    if a2 == 0.0 {
        let mean = v.iter().fold(Coord::new(0.0, 0.0), |a, c| a + *c) / n as f64;
        return (mean, 0.0);
    }
    (Coord::new(cx / (3.0 * a2), cy / (3.0 * a2)), a2 / 2.0)
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
            return Some(Point(c_ext)); // degenerate: fall back to exterior centroid
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
            // all-degenerate: mean of exterior centroids
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
                // length-weighted midpoint of the segments
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

#[cfg(test)]
mod tests {
    use super::*;

    fn poly(coords: &[(f64, f64)]) -> Polygon {
        Polygon::new(LineString::from(coords.to_vec()), vec![])
    }

    #[test]
    fn square_centroid_is_center() {
        let p = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
        let c = p.centroid().unwrap();
        assert!((c.x() - 1.0).abs() < 1e-12);
        assert!((c.y() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn centroid_is_winding_independent() {
        let ccw = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
        let cw = poly(&[(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)]);
        assert_eq!(ccw.centroid(), cw.centroid());
    }

    #[test]
    fn multipoint_centroid_is_mean() {
        let mp = MultiPoint::from(vec![(0.0, 0.0), (2.0, 0.0), (1.0, 3.0)]);
        let c = mp.centroid().unwrap();
        assert!((c.x() - 1.0).abs() < 1e-12);
        assert!((c.y() - 1.0).abs() < 1e-12);
    }
}
