//! Point-location predicates and interior sampling: classifying a point against
//! a polygon or polygon set ([`locate_in_polygon`] / [`locate_areal`]), the
//! segment/ring on-boundary tests ([`point_on_segment`], [`point_on_ring`]), and
//! the strictly-interior sample generators ([`representative_interior_point`],
//! [`sample_interior`]) the areal overlap check consumes.

use super::super::util::{orient, point_in_ring, ring_vertices};
use super::super::{Coord, Polygon};
use super::{Decomposed, Loc};

/// Classify a point against a set of polygons (areal location only).
pub(super) fn locate_areal(p: Coord, polys: &[Polygon]) -> Loc {
    let mut inside = false;
    let mut on = false;
    for poly in polys {
        match locate_in_polygon(p, poly) {
            Loc::Interior => inside = true,
            Loc::Boundary => on = true,
            Loc::Exterior => {}
        }
    }
    if inside {
        Loc::Interior
    } else if on {
        Loc::Boundary
    } else {
        Loc::Exterior
    }
}

/// Calls `f` with strictly-interior sample points of `poly` (bottom-vertex
/// bisector, centroid, and each edge midpoint nudged inward), each validated as
/// interior to `owner`'s areal region through `owner`'s prebuilt index — so no
/// throwaway per-polygon index is rebuilt. Stops as soon as `f` returns `true`,
/// so the caller's inside/outside short-circuit bounds how many samples are
/// generated (and located) instead of materializing all `~2n` up front.
pub(super) fn sample_interior(
    poly: &Polygon,
    owner: &Decomposed,
    mut f: impl FnMut(Coord) -> bool,
) {
    let v = ring_vertices(&poly.exterior().0);
    let n = v.len();
    if n < 3 {
        return;
    }
    let inside = |p: Coord| owner.locate_areal_only(p) == Loc::Interior;
    if let Some(p) = representative_interior_point(poly)
        && f(p)
    {
        return;
    }
    let (cx, cy) = (
        v.iter().map(|c| c.x).sum::<f64>() / n as f64,
        v.iter().map(|c| c.y).sum::<f64>() / n as f64,
    );
    let centroid = Coord::new(cx, cy);
    if inside(centroid) && f(centroid) {
        return;
    }
    for i in 0..n {
        let (a, b) = (v[i], v[(i + 1) % n]);
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let len = dx.hypot(dy);
        if len < 1e-12 {
            continue;
        }
        let mid = Coord::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);
        let (nx, ny) = (-dy / len, dx / len);
        let step = len * 1e-3;
        for s in [step, -step] {
            let p = Coord::new(mid.x + nx * s, mid.y + ny * s);
            if inside(p) {
                if f(p) {
                    return;
                }
                break;
            }
        }
    }
}

/// A point guaranteed strictly inside the polygon (bottom-vertex bisector).
fn representative_interior_point(poly: &Polygon) -> Option<Coord> {
    let v = ring_vertices(&poly.exterior().0);
    let n = v.len();
    if n < 3 {
        return None;
    }
    let mut k = 0;
    for i in 1..n {
        if v[i].y < v[k].y || (v[i].y == v[k].y && v[i].x < v[k].x) {
            k = i;
        }
    }
    let cur = v[k];
    let prev = v[(k + n - 1) % n];
    let next = v[(k + 1) % n];
    let u = |a: Coord| {
        let l = (a.x * a.x + a.y * a.y).sqrt();
        if l < 1e-12 {
            (0.0, 0.0)
        } else {
            (a.x / l, a.y / l)
        }
    };
    let (u1x, u1y) = u(Coord::new(prev.x - cur.x, prev.y - cur.y));
    let (u2x, u2y) = u(Coord::new(next.x - cur.x, next.y - cur.y));
    let (mut bx, mut by) = (u1x + u2x, u1y + u2y);
    let bl = (bx * bx + by * by).sqrt();
    if bl < 1e-12 {
        bx = -u1y;
        by = u1x;
    } else {
        bx /= bl;
        by /= bl;
    }
    let elen = (prev.x - cur.x)
        .hypot(prev.y - cur.y)
        .min((next.x - cur.x).hypot(next.y - cur.y));
    let step = elen * 1e-3;
    let cand = Coord::new(cur.x + bx * step, cur.y + by * step);
    if locate_in_polygon(cand, poly) == Loc::Interior {
        Some(cand)
    } else {
        let alt = Coord::new(cur.x - bx * step, cur.y - by * step);
        (locate_in_polygon(alt, poly) == Loc::Interior).then_some(alt)
    }
}

/// Classify a point against a single polygon (ring boundary vs interior).
pub(super) fn locate_in_polygon(p: Coord, poly: &Polygon) -> Loc {
    if point_on_ring(p, &poly.exterior().0) {
        return Loc::Boundary;
    }
    for h in poly.interiors() {
        if point_on_ring(p, &h.0) {
            return Loc::Boundary;
        }
    }
    if !point_in_ring(p, &poly.exterior().0) {
        return Loc::Exterior;
    }
    for h in poly.interiors() {
        if point_in_ring(p, &h.0) {
            return Loc::Exterior;
        }
    }
    Loc::Interior
}

fn on_bbox(a: Coord, b: Coord, p: Coord) -> bool {
    let eps = 1e-9 * (1.0 + a.x.abs() + b.x.abs() + a.y.abs() + b.y.abs());
    p.x <= a.x.max(b.x) + eps
        && p.x >= a.x.min(b.x) - eps
        && p.y <= a.y.max(b.y) + eps
        && p.y >= a.y.min(b.y) - eps
}

/// Is `p` on segment `a–b`? Collinearity uses a scale-relative tolerance so
/// computed intersection points and midpoints land on the edges they came from.
pub(super) fn point_on_segment(p: Coord, a: Coord, b: Coord) -> bool {
    let o = orient(a, b, p);
    let scale = (b.x - a.x).hypot(b.y - a.y) * (p.x - a.x).hypot(p.y - a.y);
    o.abs() <= 1e-9 * scale.max(1.0) && on_bbox(a, b, p)
}

pub(super) fn point_on_segment_interior(p: Coord, a: Coord, b: Coord) -> bool {
    p != a && p != b && point_on_segment(p, a, b)
}

fn point_on_ring(p: Coord, ring: &[Coord]) -> bool {
    let v = ring_vertices(ring);
    let n = v.len();
    for i in 0..n {
        if point_on_segment(p, v[i], v[(i + 1) % n]) {
            return true;
        }
    }
    false
}
