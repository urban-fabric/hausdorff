//! Point-location predicates and interior sampling: classifying a point against
//! a polygon or polygon set ([`locate_in_polygon`] / [`locate_areal`]), the
//! segment/ring on-boundary tests ([`point_on_segment`], [`point_on_ring`]), and
//! the strictly-interior sample generators ([`representative_interior_point`],
//! [`sample_interior`]) the areal overlap check consumes.

use super::{Decomposed, Loc};
use crate::space::r2::util::{orient, point_in_ring, ring_vertices};
use crate::space::r2::{Coord, Polygon};

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

/// The largest coordinate magnitude among the arguments — the scale at which
/// forming differences between them rounds.
fn magnitude(a: Coord, b: Coord, p: Coord) -> f64 {
    a.x.abs()
        .max(a.y.abs())
        .max(b.x.abs())
        .max(b.y.abs())
        .max(p.x.abs())
        .max(p.y.abs())
}

/// Is `p` within the bounding box of segment `a–b`, allowing for rounding?
///
/// **The slack is relative to the segment's own span, not to how far the segment
/// is from the origin.** It was `1e-9 · (|aₓ| + |bₓ| + |a_y| + |b_y|)`, which is a
/// *length* scaled by position: at a coordinate magnitude of `1e10` that is a box
/// grown by 40 units, so a vertex a whole unit beyond the segment's end still
/// counted as inside it. Two triangles meeting at a single point on the line
/// `y = x` therefore had each one's far vertex read as lying on the other's edge,
/// `relate` reported a 1-dimensional boundary intersection, and a valid
/// `MultiPolygon` was rejected with "polygons touch on a line" — which in turn
/// made `make_valid` fail with `StillInvalid` on geometry it had repaired
/// correctly.
///
/// The second term is the part a position-scaled slack was reaching for and
/// getting wrong by seven orders of magnitude: a *computed* point (an
/// intersection, a midpoint) can land just outside the box it belongs to, and that
/// error is a few ULPs at the coordinates' magnitude — `8 · ε · |coord|`, the same
/// bound the arrangement's vertex pool snaps with, not `1e-9 · |coord|`.
fn on_bbox(a: Coord, b: Coord, p: Coord) -> bool {
    let span = (b.x - a.x).abs().max((b.y - a.y).abs());
    let eps = (1e-9 * span).max(8.0 * f64::EPSILON * magnitude(a, b, p));
    p.x <= a.x.max(b.x) + eps
        && p.x >= a.x.min(b.x) - eps
        && p.y <= a.y.max(b.y) + eps
        && p.y >= a.y.min(b.y) - eps
}

/// Is `p` on segment `a–b`? Collinearity uses a scale-relative tolerance so
/// computed intersection points and midpoints land on the edges they came from.
///
/// **Both terms of the bound are areas, because `orient` is one.** The previous
/// bound was `1e-9 · max(span · reach, 1.0)`, and that floor of `1.0` is an
/// *absolute* area: for a shape whose features are smaller than about `3e-5` units
/// the floor dominates, every `orient` falls under it, and every triple of points
/// reads as collinear. Features that small are ordinary — `1e-5` of a degree is
/// about a metre — and the same two triangles were rejected at a feature size of
/// `1e-9`.
///
/// So the floor is now the rounding error actually incurred rather than a
/// constant: forming `b − a` and `p − a` rounds at `ε · |coord|`, and those
/// differences are then multiplied by something of size `span` or `reach`, which
/// makes the resulting slack an area that shrinks with the geometry instead of
/// swallowing it.
pub(super) fn point_on_segment(p: Coord, a: Coord, b: Coord) -> bool {
    let o = orient(a, b, p);
    let span = (b.x - a.x).hypot(b.y - a.y);
    let reach = (p.x - a.x).hypot(p.y - a.y);
    let eps = 1e-9 * span * reach + 8.0 * f64::EPSILON * magnitude(a, b, p) * span.max(reach);
    o.abs() <= eps && on_bbox(a, b, p)
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
