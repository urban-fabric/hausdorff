//! `r2::min_distance` — the **minimum** Euclidean distance between two whole
//! geometries: the `ST_Distance` kernel.
//!
//! This is distinct from the other distance notions in the crate:
//! - [`Distance`](super::Distance) ([`line_measures`](super::line_measures)) is
//!   point-to-point only.
//! - [`set_distance`](super::set_distance) gives the *max-min* Hausdorff and the
//!   Fréchet curve distance.
//!
//! Here we want the smallest separation between two geometries — `0` if they
//! intersect, otherwise the closest approach of their edges/vertices. It is the
//! kernel galois wraps as `ST_Distance` / `ST_DWithin` / `ST_ClosestPoint`.
//!
//! Every geometry is reduced to a set of bounded segments (a 0-dimensional point
//! is a degenerate segment `p–p`), so the point/segment/ring cases all collapse
//! into one segment-to-segment routine. Brute force is `O(n·m)` in the edge
//! counts — correct and competitive with `geo` for moderate inputs. A spatial
//! index (AABB / STR-tree, à la CGAL/GEOS) would accelerate large geometries;
//! the segment decomposition here is the substrate such an index would sit on.

use super::{Coord, Geometry, Intersects, LineString, Polygon, Rect, StrTree, Triangle};

/// Minimum Euclidean distance between two geometries. `0.0` if they intersect;
/// `f64::INFINITY` if either geometry is empty.
///
/// Works in *squared* distance internally (one final `sqrt`), skips the
/// intersection test when the bounding boxes are disjoint, and prunes segment
/// pairs whose bounding boxes are already farther than the running best.
#[must_use]
pub fn euclidean_distance(a: &Geometry, b: &Geometry) -> f64 {
    let (sa, sb) = (segments(a), segments(b));
    if sa.is_empty() || sb.is_empty() {
        return f64::INFINITY;
    }
    let (ba, bb) = (bbox_of(&sa), bbox_of(&sb));
    // The full intersection test is only needed when the boxes overlap; disjoint
    // boxes can neither cross nor contain, so the gap is the boundary distance.
    if boxes_overlap(ba, bb) && a.intersects(b) {
        return 0.0;
    }
    // A tiny cost model picks the cheaper strategy from cheap input features (see
    // [`prefer_indexed`]): brute-with-prune for separated/small operands, an
    // STR-tree for large operands whose boxes overlap (where the prune is weak).
    let min_sq = if prefer_indexed(sa.len(), sb.len(), box_gap_sq(ba, bb) == 0.0) {
        // Index the larger segment set, iterate the smaller (fewer queries).
        if sa.len() <= sb.len() {
            indexed_min_sq(&sa, &sb)
        } else {
            indexed_min_sq(&sb, &sa)
        }
    } else {
        brute_min_sq(&sa, bb, &sb)
    };
    min_sq.sqrt()
}

/// Cost model: should `euclidean_distance` use the spatial index instead of the
/// brute scan?
///
/// The brute scan's bbox prune is near-linear when the operands are **separated**
/// (disjoint bounding boxes), so the index can only win when the boxes overlap —
/// then the prune is weak and the scan tends to `O(n·m)`. There, the index's
/// `~6·m` build + `~6·n·log m` query beats `n·m` once both operands are large.
/// The constants are rough relative weights; the comparison, not their exact
/// values, drives the choice.
fn prefer_indexed(n: usize, m: usize, boxes_overlap: bool) -> bool {
    if !boxes_overlap {
        return false;
    }
    let weights = crate::config::config().cost_model.min_distance;
    let (nf, mf) = (n as f64, m as f64);
    let brute_cost = nf * mf;
    let index_cost =
        weights.index_build_weight * mf + weights.index_query_weight * nf * (m.max(2) as f64).log2();
    index_cost < brute_cost
}

/// Brute-force squared min gap with bbox pruning.
fn brute_min_sq(sa: &[(Coord, Coord)], bb: (Coord, Coord), sb: &[(Coord, Coord)]) -> f64 {
    // Precompute B's segment boxes once (reused across every A-segment).
    let sb_boxes: Vec<(Coord, Coord)> = sb.iter().map(|&(p, q)| seg_box(p, q)).collect();
    let mut min_sq = f64::INFINITY;
    for &(a1, a2) in sa {
        let abox = seg_box(a1, a2);
        // Coarse prune: skip this whole A-segment if it's farther than the best
        // from B's entire bounding box.
        if box_gap_sq(abox, bb) >= min_sq {
            continue;
        }
        for (j, &(b1, b2)) in sb.iter().enumerate() {
            if box_gap_sq(abox, sb_boxes[j]) >= min_sq {
                continue; // this pair can't beat the running best
            }
            min_sq = min_sq.min(segment_gap_sq(a1, a2, b1, b2));
            if min_sq == 0.0 {
                return 0.0;
            }
        }
    }
    min_sq
}

/// Whether two geometries are within `d` of each other — the `ST_DWithin`
/// kernel. Early-exits as soon as a pair is within `d`, so it is typically much
/// faster than computing the full [`euclidean_distance`].
#[must_use]
pub fn euclidean_within(a: &Geometry, b: &Geometry, d: f64) -> bool {
    if d < 0.0 {
        return false;
    }
    let (sa, sb) = (segments(a), segments(b));
    if sa.is_empty() || sb.is_empty() {
        return false;
    }
    let (ba, bb) = (bbox_of(&sa), bbox_of(&sb));
    if boxes_overlap(ba, bb) && a.intersects(b) {
        return true;
    }
    let d_sq = d * d;
    let sb_boxes: Vec<(Coord, Coord)> = sb.iter().map(|&(p, q)| seg_box(p, q)).collect();
    for &(a1, a2) in &sa {
        let abox = seg_box(a1, a2);
        if box_gap_sq(abox, bb) > d_sq {
            continue;
        }
        for (j, &(b1, b2)) in sb.iter().enumerate() {
            if box_gap_sq(abox, sb_boxes[j]) > d_sq {
                continue;
            }
            if segment_gap_sq(a1, a2, b1, b2) <= d_sq {
                return true;
            }
        }
    }
    false
}

/// Index `b_segs` in an STR-tree, then for each A-segment query only the B
/// segments within the running best distance — `O(n·log m)` for large inputs.
fn indexed_min_sq(a_segs: &[(Coord, Coord)], b_segs: &[(Coord, Coord)]) -> f64 {
    let items: Vec<(Rect, usize)> = b_segs
        .iter()
        .enumerate()
        .map(|(i, &(p, q))| {
            let (lo, hi) = seg_box(p, q);
            (Rect::new(lo, hi), i)
        })
        .collect();
    let tree = StrTree::new(items);
    let mut min_sq = f64::INFINITY;
    for &(a1, a2) in a_segs {
        let (alo, ahi) = seg_box(a1, a2);
        // Query A's box expanded by the running best distance; B segments outside
        // it can't beat it. (Infinite on the first A-segment → visits all of B.)
        let d = min_sq.sqrt();
        let q = Rect::new(
            Coord::new(alo.x - d, alo.y - d),
            Coord::new(ahi.x + d, ahi.y + d),
        );
        let mut local = min_sq;
        tree.query_each(q, |_, &j| {
            let (b1, b2) = b_segs[j];
            local = local.min(segment_gap_sq(a1, a2, b1, b2));
        });
        min_sq = local;
        if min_sq == 0.0 {
            return 0.0;
        }
    }
    min_sq
}

/// Distance from point `p` to the bounded segment `a`–`b` (the projection is
/// clamped to the segment).
#[must_use]
pub fn point_segment_distance(p: Coord, a: Coord, b: Coord) -> f64 {
    point_segment_distance_sq(p, a, b).sqrt()
}

/// Minimum distance between two bounded segments — `0.0` if they intersect. This
/// is the standalone primitive: it performs the crossing test itself, unlike the
/// internal squared-gap routine used on the already-non-intersecting hot path.
#[must_use]
pub fn segment_segment_distance(a1: Coord, a2: Coord, b1: Coord, b2: Coord) -> f64 {
    if segments_intersect(a1, a2, b1, b2) {
        return 0.0;
    }
    segment_gap_sq(a1, a2, b1, b2).sqrt()
}

/// Squared distance from `p` to the bounded segment `a`–`b` (no `sqrt`).
fn point_segment_distance_sq(p: Coord, a: Coord, b: Coord) -> f64 {
    let ab = b - a;
    let len2 = ab.x * ab.x + ab.y * ab.y;
    let proj = if len2 == 0.0 {
        a
    } else {
        let t = (((p - a).x * ab.x + (p - a).y * ab.y) / len2).clamp(0.0, 1.0);
        a + ab * t
    };
    let d = p - proj;
    d.x * d.x + d.y * d.y
}

/// Smallest **squared** endpoint-to-segment gap between two non-crossing segments.
pub(crate) fn segment_gap_sq(a1: Coord, a2: Coord, b1: Coord, b2: Coord) -> f64 {
    point_segment_distance_sq(a1, b1, b2)
        .min(point_segment_distance_sq(a2, b1, b2))
        .min(point_segment_distance_sq(b1, a1, a2))
        .min(point_segment_distance_sq(b2, a1, a2))
}

// ---- axis-aligned bounding-box helpers (for pruning) -----------------------

/// `(min, max)` corner of all endpoints in a segment list (caller ensures
/// non-empty).
fn bbox_of(segs: &[(Coord, Coord)]) -> (Coord, Coord) {
    let mut lo = Coord::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Coord::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(p, q) in segs {
        for c in [p, q] {
            lo = Coord::new(lo.x.min(c.x), lo.y.min(c.y));
            hi = Coord::new(hi.x.max(c.x), hi.y.max(c.y));
        }
    }
    (lo, hi)
}

/// `(min, max)` corner of a single segment.
pub(crate) fn seg_box(a: Coord, b: Coord) -> (Coord, Coord) {
    (
        Coord::new(a.x.min(b.x), a.y.min(b.y)),
        Coord::new(a.x.max(b.x), a.y.max(b.y)),
    )
}

/// Do two AABBs overlap (touching counts)?
fn boxes_overlap(a: (Coord, Coord), b: (Coord, Coord)) -> bool {
    a.0.x <= b.1.x && b.0.x <= a.1.x && a.0.y <= b.1.y && b.0.y <= a.1.y
}

/// Squared gap between two AABBs (`0` if they overlap) — a lower bound on the
/// distance between anything inside them.
fn box_gap_sq(a: (Coord, Coord), b: (Coord, Coord)) -> f64 {
    let dx = (b.0.x - a.1.x).max(a.0.x - b.1.x).max(0.0);
    let dy = (b.0.y - a.1.y).max(a.0.y - b.1.y).max(0.0);
    dx * dx + dy * dy
}

/// Do the bounded segments `p1p2` and `p3p4` intersect (cross or touch)?
fn segments_intersect(p1: Coord, p2: Coord, p3: Coord, p4: Coord) -> bool {
    let d1 = orient(p3, p4, p1);
    let d2 = orient(p3, p4, p2);
    let d3 = orient(p1, p2, p3);
    let d4 = orient(p1, p2, p4);
    if (d1 > 0.0) != (d2 > 0.0) && (d3 > 0.0) != (d4 > 0.0) {
        return true; // proper crossing
    }
    (d1 == 0.0 && in_box(p3, p4, p1))
        || (d2 == 0.0 && in_box(p3, p4, p2))
        || (d3 == 0.0 && in_box(p1, p2, p3))
        || (d4 == 0.0 && in_box(p1, p2, p4))
}

/// Twice the signed area of triangle `abc` (orientation: >0 left, <0 right).
fn orient(a: Coord, b: Coord, c: Coord) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// Is `p` inside the bounding box of segment `a`–`b` (used for collinear cases)?
fn in_box(a: Coord, b: Coord, p: Coord) -> bool {
    p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
}

// ---- geometry → bounded segments -------------------------------------------

pub(crate) fn segments(g: &Geometry) -> Vec<(Coord, Coord)> {
    let mut out = Vec::new();
    push_segments(g, &mut out);
    out
}

fn push_segments(g: &Geometry, out: &mut Vec<(Coord, Coord)>) {
    match g {
        Geometry::Point(p) => {
            let c = Coord::new(p.x(), p.y());
            out.push((c, c));
        }
        Geometry::Line(l) => out.push((l.start, l.end)),
        Geometry::LineString(ls) => push_line_string(ls, out),
        Geometry::Polygon(poly) => push_polygon(poly, out),
        Geometry::MultiPoint(mp) => {
            for p in mp.iter() {
                let c = Coord::new(p.x(), p.y());
                out.push((c, c));
            }
        }
        Geometry::MultiLineString(mls) => {
            for ls in mls.iter() {
                push_line_string(ls, out);
            }
        }
        Geometry::MultiPolygon(mp) => {
            for poly in mp.iter() {
                push_polygon(poly, out);
            }
        }
        Geometry::Rect(r) => push_rect(r, out),
        Geometry::Triangle(t) => push_triangle(t, out),
        Geometry::GeometryCollection(c) => {
            for g in c.iter() {
                push_segments(g, out);
            }
        }
    }
}

fn push_line_string(ls: &LineString, out: &mut Vec<(Coord, Coord)>) {
    match ls.0.as_slice() {
        [] => {}
        [only] => out.push((*only, *only)),
        cs => {
            for w in cs.windows(2) {
                out.push((w[0], w[1]));
            }
        }
    }
}

fn push_polygon(poly: &Polygon, out: &mut Vec<(Coord, Coord)>) {
    push_line_string(poly.exterior(), out);
    for hole in poly.interiors() {
        push_line_string(hole, out);
    }
}

fn push_rect(r: &Rect, out: &mut Vec<(Coord, Coord)>) {
    let (lo, hi) = (r.min(), r.max());
    let corners = [
        Coord::new(lo.x, lo.y),
        Coord::new(hi.x, lo.y),
        Coord::new(hi.x, hi.y),
        Coord::new(lo.x, hi.y),
    ];
    for i in 0..4 {
        out.push((corners[i], corners[(i + 1) % 4]));
    }
}

fn push_triangle(t: &Triangle, out: &mut Vec<(Coord, Coord)>) {
    out.push((t.0, t.1));
    out.push((t.1, t.2));
    out.push((t.2, t.0));
}
