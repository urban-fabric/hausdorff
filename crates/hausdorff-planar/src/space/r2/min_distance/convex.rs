//! The `O(n+m)` convex fast path for [`min_distance`](super): when both operands
//! are convex areal shapes, their gap is the distance from the origin to their
//! Minkowski difference, found by the 2-simplex GJK distance algorithm — no
//! `O(n·m)` segment scan. Falls back (returns `None`) if GJK doesn't converge.

use hausdorff_traits::is_convex_ring;

use super::super::{Coord, Geometry, Intersects};
use super::boxes_overlap;

/// Resolve the minimum distance between `a` and `b` via the convex fast path, or
/// `None` if it doesn't apply — either operand non-convex, or GJK failed to
/// converge — in which case the caller falls back to the exact segment scan.
///
/// Both operands convex ⇒ the gap is found in `O(n+m)` with no segment
/// decomposition: a disjoint bounding-box pair is certainly separated, while an
/// overlapping pair needs the exact [`Intersects`] test to rule out a true
/// intersection (gap `0`) before GJK runs.
///
/// Proof-backed by `riemann` `Hausdorff/MinDistance` — the *premise* GJK computes
/// against (not its convergence): the gap equals the origin's closest approach to
/// the **Minkowski difference** `A ⊖ B` (`gap_eq_minkowski`), the origin lies in
/// `A ⊖ B` iff the shapes intersect (`zero_mem_sub_iff`, justifying the
/// `intersects ⇒ Some(0.0)` short-circuit via `gap_eq_zero_of_mem`), and the gap
/// under-bounds every realised pair distance (`gap_le_dist`) — a disjoint-box
/// pair is a genuine separation certificate.
pub(super) fn convex_gap(a: &Geometry, b: &Geometry) -> Option<f64> {
    let (ca, cb) = (convex_shape(a)?, convex_shape(b)?);
    if boxes_overlap(loop_bbox(ca.verts()), loop_bbox(cb.verts())) && a.intersects(b) {
        return Some(0.0);
    }
    convex_min_distance(ca.verts(), cb.verts())
}

/// `(min, max)` corner of a flat vertex loop (caller ensures non-empty).
fn loop_bbox(verts: &[Coord]) -> (Coord, Coord) {
    let mut lo = Coord::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Coord::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &c in verts {
        lo = Coord::new(lo.x.min(c.x), lo.y.min(c.y));
        hi = Coord::new(hi.x.max(c.x), hi.y.max(c.y));
    }
    (lo, hi)
}

/// A convex operand's boundary vertices for the fast path. Polygon rings are
/// borrowed (no copy); the small fixed shapes are materialised on the stack.
enum Convex<'a> {
    Borrowed(&'a [Coord]),
    Owned([Coord; 4], usize),
}

impl Convex<'_> {
    /// The open vertex loop (no repeated closing vertex).
    fn verts(&self) -> &[Coord] {
        match self {
            Convex::Borrowed(s) => s,
            Convex::Owned(buf, n) => &buf[..*n],
        }
    }
}

/// If `g` is a single, non-degenerate **convex** areal shape, return its
/// boundary vertices. `None` for anything holed, non-areal, degenerate, or
/// non-convex — those take the brute/index path. Allocation-free: polygon rings
/// are borrowed in place.
fn convex_shape(g: &Geometry) -> Option<Convex<'_>> {
    let shape = match g {
        Geometry::Polygon(p) if p.interiors().is_empty() => {
            let ring = p.exterior().0.as_slice();
            let end = usize::from(ring.len() >= 2 && ring.first() == ring.last());
            Convex::Borrowed(&ring[..ring.len() - end])
        }
        Geometry::Triangle(t) => Convex::Owned([t.0, t.1, t.2, t.0], 3),
        Geometry::Rect(r) => {
            let (lo, hi) = (r.min(), r.max());
            Convex::Owned([lo, Coord::new(hi.x, lo.y), hi, Coord::new(lo.x, hi.y)], 4)
        }
        _ => return None,
    };
    // A non-degenerate convex polygon's vertex hull equals the shape — exactly
    // what GJK on the Minkowski difference needs. The consistent-turn-sign test
    // is the generic `Orientation` convexity classifier (shared with `s2::Loop`
    // and `r2::ConvexPolygon`).
    is_convex_ring(shape.verts()).then_some(shape)
}

/// Minimum distance between two **disjoint** convex vertex loops via the
/// 2-simplex GJK distance algorithm on their Minkowski difference `A ⊖ B`:
/// the gap equals the distance from the origin to that convex set. `O(n+m)`.
///
/// Returns `None` if it fails to converge within a generous iteration budget, so
/// the caller can fall back to the exact brute/index scan — correctness never
/// depends on GJK converging.
fn convex_min_distance(va: &[Coord], vb: &[Coord]) -> Option<f64> {
    let mink = |d: Coord| support(va, d) - support(vb, -d);

    let mut dir = vb[0] - va[0];
    if dir.x == 0.0 && dir.y == 0.0 {
        dir = Coord::new(1.0, 0.0);
    }
    let mut a = mink(dir);
    let mut b = mink(-dir);
    let budget = 2 * (va.len() + vb.len()) + 8;
    for _ in 0..budget {
        let p = closest_on_segment_to_origin(a, b);
        let p_len2 = p.dot(p);
        if p_len2 <= 1e-24 {
            return Some(0.0);
        }
        let p_len = p_len2.sqrt();
        let unit = p * (-1.0 / p_len);
        let c = mink(unit);
        let progress = c.dot(unit) - p.dot(unit);
        if progress < 1e-12 * (1.0 + p_len) {
            return Some(p_len);
        }
        if a.dot(a) <= b.dot(b) {
            b = c;
        } else {
            a = c;
        }
    }
    None
}

/// Vertex of `verts` farthest in direction `d` (the support function).
fn support(verts: &[Coord], d: Coord) -> Coord {
    let mut best = verts[0];
    let mut best_dot = best.dot(d);
    for &v in &verts[1..] {
        let vd = v.dot(d);
        if vd > best_dot {
            best_dot = vd;
            best = v;
        }
    }
    best
}

/// Closest point to the origin on the bounded segment `a`–`b`.
fn closest_on_segment_to_origin(a: Coord, b: Coord) -> Coord {
    let ab = b - a;
    let len2 = ab.dot(ab);
    if len2 <= 0.0 {
        return a;
    }
    let t = ((-a).dot(ab) / len2).clamp(0.0, 1.0);
    a + ab * t
}
