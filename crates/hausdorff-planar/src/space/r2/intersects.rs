//! Native intersection predicate — do two geometries share at least one point?
//! The in-house equivalent of `geo`'s `Intersects` on `Geometry`, built from
//! robust segment crossing + point-in-polygon (no DE-9IM matrix needed).

use super::util::{point_in_ring, ring_vertices};
use super::{Coord, CoordsIter, Geometry, Polygon, Rect};

/// True iff `self` and `other` share at least one point.
pub trait Intersects {
    fn intersects(&self, other: &Self) -> bool;
}

impl Intersects for Rect {
    /// Axis-aligned rectangle overlap (touching counts) — the broad-phase test
    /// behind range queries and spatial joins.
    fn intersects(&self, other: &Rect) -> bool {
        self.min().x <= other.max().x
            && other.min().x <= self.max().x
            && self.min().y <= other.max().y
            && other.min().y <= self.max().y
    }
}

impl Intersects for Geometry {
    fn intersects(&self, other: &Geometry) -> bool {
        let (sa, sb) = (segments(self), segments(other));
        for &(a1, a2) in &sa {
            for &(b1, b2) in &sb {
                if segments_intersect(a1, a2, b1, b2) {
                    return true;
                }
            }
        }
        let (pa, pb) = (polygons(self), polygons(other));
        if !pb.is_empty() {
            for c in self.coords_iter() {
                if pb.iter().any(|p| poly_contains(p, c)) {
                    return true;
                }
            }
        }
        if !pa.is_empty() {
            for c in other.coords_iter() {
                if pa.iter().any(|p| poly_contains(p, c)) {
                    return true;
                }
            }
        }
        let ca: Vec<Coord> = self.coords_iter().collect();
        let cb: Vec<Coord> = other.coords_iter().collect();
        if sa.is_empty() || sb.is_empty() {
            for &x in &ca {
                if cb.contains(&x) {
                    return true;
                }
                if sb.iter().any(|&(b1, b2)| point_on_segment(x, b1, b2)) {
                    return true;
                }
            }
            for &x in &cb {
                if sa.iter().any(|&(a1, a2)| point_on_segment(x, a1, a2)) {
                    return true;
                }
            }
        }
        false
    }
}

fn orient(a: Coord, b: Coord, c: Coord) -> i32 {
    let s = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    if s > 0.0 {
        1
    } else if s < 0.0 {
        -1
    } else {
        0
    }
}

fn on_bbox(a: Coord, b: Coord, p: Coord) -> bool {
    p.x <= a.x.max(b.x) && p.x >= a.x.min(b.x) && p.y <= a.y.max(b.y) && p.y >= a.y.min(b.y)
}

/// Robust segment intersection (CLRS), including collinear/touching cases.
pub(crate) fn segments_intersect(a1: Coord, a2: Coord, b1: Coord, b2: Coord) -> bool {
    let d1 = orient(b1, b2, a1);
    let d2 = orient(b1, b2, a2);
    let d3 = orient(a1, a2, b1);
    let d4 = orient(a1, a2, b2);
    if d1 != d2 && d3 != d4 {
        return true;
    }
    (d1 == 0 && on_bbox(b1, b2, a1))
        || (d2 == 0 && on_bbox(b1, b2, a2))
        || (d3 == 0 && on_bbox(a1, a2, b1))
        || (d4 == 0 && on_bbox(a1, a2, b2))
}

fn point_on_segment(p: Coord, a: Coord, b: Coord) -> bool {
    orient(a, b, p) == 0 && on_bbox(a, b, p)
}

fn poly_contains(p: &Polygon, c: Coord) -> bool {
    point_in_ring(c, &p.exterior().0) && !p.interiors().iter().any(|h| point_in_ring(c, &h.0))
}

fn ring_segments(coords: &[Coord], out: &mut Vec<(Coord, Coord)>) {
    let v = ring_vertices(coords);
    let n = v.len();
    if n < 2 {
        return;
    }
    for i in 0..n {
        out.push((v[i], v[(i + 1) % n]));
    }
}

fn line_segments(coords: &[Coord], out: &mut Vec<(Coord, Coord)>) {
    for w in coords.windows(2) {
        out.push((w[0], w[1]));
    }
}

fn poly_segments(p: &Polygon, out: &mut Vec<(Coord, Coord)>) {
    ring_segments(&p.exterior().0, out);
    for h in p.interiors() {
        ring_segments(&h.0, out);
    }
}

fn segments(g: &Geometry) -> Vec<(Coord, Coord)> {
    let mut out = Vec::new();
    match g {
        Geometry::Point(_) | Geometry::MultiPoint(_) => {}
        Geometry::Line(l) => out.push((l.start, l.end)),
        Geometry::LineString(ls) => line_segments(&ls.0, &mut out),
        Geometry::MultiLineString(mls) => {
            for ls in mls.iter() {
                line_segments(&ls.0, &mut out);
            }
        }
        Geometry::Polygon(p) => poly_segments(p, &mut out),
        Geometry::MultiPolygon(mp) => {
            for p in mp.iter() {
                poly_segments(p, &mut out);
            }
        }
        Geometry::Rect(r) => poly_segments(&r.to_polygon(), &mut out),
        Geometry::Triangle(t) => poly_segments(&t.to_polygon(), &mut out),
        Geometry::GeometryCollection(gs) => {
            for sub in gs {
                out.extend(segments(sub));
            }
        }
    }
    out
}

fn polygons(g: &Geometry) -> Vec<Polygon> {
    match g {
        Geometry::Polygon(p) => vec![p.clone()],
        Geometry::MultiPolygon(mp) => mp.0.clone(),
        Geometry::Rect(r) => vec![r.to_polygon()],
        Geometry::Triangle(t) => vec![t.to_polygon()],
        Geometry::GeometryCollection(gs) => gs.iter().flat_map(polygons).collect(),
        _ => Vec::new(),
    }
}
