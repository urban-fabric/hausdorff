//! Ramer–Douglas–Peucker line/polygon simplification — mirror of `geo`'s
//! `Simplify`. Drops vertices whose perpendicular distance to the working chord
//! (a bounded segment) is at or below `epsilon`, in the geometry's own planar
//! units.

use super::{Coord, LineString, MultiLineString, MultiPolygon, Polygon};

/// Simplify a geometry with Douglas–Peucker at tolerance `epsilon`.
pub trait Simplify {
    fn simplify(&self, epsilon: f64) -> Self;
}

/// Distance from `p` to the bounded segment `a`–`b` (matching geo's
/// point-to-`Line` metric, which clamps the projection to the segment).
fn point_segment_distance(p: Coord, a: Coord, b: Coord) -> f64 {
    let ab = b - a;
    let len2 = ab.x * ab.x + ab.y * ab.y;
    let proj = if len2 == 0.0 {
        a
    } else {
        let t = (((p - a).x * ab.x + (p - a).y * ab.y) / len2).clamp(0.0, 1.0);
        a + ab * t
    };
    let d = p - proj;
    (d.x * d.x + d.y * d.y).sqrt()
}

/// Douglas–Peucker over a coordinate run; endpoints are always retained.
fn rdp(coords: &[Coord], epsilon: f64) -> Vec<Coord> {
    if coords.len() < 3 {
        return coords.to_vec();
    }
    let (first, last) = (coords[0], coords[coords.len() - 1]);

    let mut max_dist = 0.0;
    let mut index = 0;
    for (i, &c) in coords.iter().enumerate().take(coords.len() - 1).skip(1) {
        let d = point_segment_distance(c, first, last);
        if d > max_dist {
            max_dist = d;
            index = i;
        }
    }

    if max_dist > epsilon {
        let mut left = rdp(&coords[..=index], epsilon);
        let right = rdp(&coords[index..], epsilon);
        left.pop();
        left.extend(right);
        left
    } else {
        vec![first, last]
    }
}

impl Simplify for LineString {
    fn simplify(&self, epsilon: f64) -> Self {
        LineString::new(rdp(&self.0, epsilon))
    }
}

impl Simplify for MultiLineString {
    fn simplify(&self, epsilon: f64) -> Self {
        MultiLineString::new(self.iter().map(|ls| ls.simplify(epsilon)).collect())
    }
}

impl Simplify for Polygon {
    fn simplify(&self, epsilon: f64) -> Self {
        Polygon::new(
            self.exterior().simplify(epsilon),
            self.interiors()
                .iter()
                .map(|r| r.simplify(epsilon))
                .collect(),
        )
    }
}

impl Simplify for MultiPolygon {
    fn simplify(&self, epsilon: f64) -> Self {
        MultiPolygon::new(self.iter().map(|p| p.simplify(epsilon)).collect())
    }
}
