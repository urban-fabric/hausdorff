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
//! counts — correct and competitive with `geo` for moderate inputs. Two
//! accelerators sit on top of it:
//!
//! - **Convex fast path:** when both operands are convex areal shapes, their gap
//!   is found in `O(n+m)` by GJK on the Minkowski difference (no `O(n·m)` scan).
//! - **Spatial index:** for large operands whose bounding boxes overlap (where
//!   the bbox prune is weak), an [`StrTree`] over one operand's segments turns the
//!   scan into `O(n·log m)`. A tiny cost model picks between brute and index.
//!
//! Repeated queries (a spatial join, a batch of distances) can reuse a
//! [`MinDistanceWorkspace`] to keep the per-call segment/bbox buffers allocated
//! across calls. (For the index-once-query-many shape, see
//! [`PreparedGeometry`](super::PreparedGeometry).)

use super::util::orient;
use super::{Coord, Geometry, Intersects, LineString, Polygon, Rect, StrTree, Triangle};

/// Minimum Euclidean distance between two geometries. `0.0` if they intersect;
/// `f64::INFINITY` if either geometry is empty.
///
/// Works in *squared* distance internally (one final `sqrt`), takes an `O(n+m)`
/// convex fast path when both operands are convex, skips the intersection test
/// when the bounding boxes are disjoint, and prunes segment pairs whose bounding
/// boxes are already farther than the running best.
#[must_use]
pub fn euclidean_distance(a: &Geometry, b: &Geometry) -> f64 {
    MinDistanceWorkspace::new().distance(a, b)
}

/// Whether two geometries are within `d` of each other — the `ST_DWithin`
/// kernel. Early-exits as soon as a pair is within `d`, so it is typically much
/// faster than computing the full [`euclidean_distance`].
#[must_use]
pub fn euclidean_within(a: &Geometry, b: &Geometry, d: f64) -> bool {
    MinDistanceWorkspace::new().within(a, b, d)
}

/// Reusable scratch buffers for the minimum-distance kernel. A batch of queries
/// (a spatial join, a fan of distances) can reuse one workspace so the per-call
/// segment and bounding-box vectors are allocated once instead of per call —
/// the throughput win the single-shot [`euclidean_distance`] can't capture.
///
/// ```
/// use hausdorff::space::r2::MinDistanceWorkspace;
/// use hausdorff::{Geometry, Point};
///
/// let mut ws = MinDistanceWorkspace::new();
/// let target = Geometry::Point(Point::new(0.0, 0.0));
/// for k in 1..=3 {
///     let p = Geometry::Point(Point::new(k as f64, 0.0));
///     assert_eq!(ws.distance(&target, &p), k as f64);
/// }
/// ```
#[derive(Debug, Default, Clone)]
pub struct MinDistanceWorkspace {
    sa: Vec<(Coord, Coord)>,
    sb: Vec<(Coord, Coord)>,
    sb_boxes: Vec<(Coord, Coord)>,
}

impl MinDistanceWorkspace {
    /// A new workspace with empty buffers.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Minimum Euclidean distance between `a` and `b`, reusing this workspace's
    /// buffers. Equivalent to [`euclidean_distance`] but without re-allocating.
    #[must_use]
    pub fn distance(&mut self, a: &Geometry, b: &Geometry) -> f64 {
        // Convex fast path first, so two convex operands resolve their gap in
        // O(n+m) and skip the segment decomposition (and its buffers) entirely.
        if let Some(d) = convex_gap(a, b) {
            return d;
        }
        self.sa.clear();
        self.sb.clear();
        push_segments(a, &mut self.sa);
        push_segments(b, &mut self.sb);
        if self.sa.is_empty() || self.sb.is_empty() {
            return f64::INFINITY;
        }
        let (ba, bb) = (bbox_of(&self.sa), bbox_of(&self.sb));
        // The full intersection test is only needed when the boxes overlap;
        // disjoint boxes can neither cross nor contain, so the gap is the
        // boundary distance.
        if boxes_overlap(ba, bb) && a.intersects(b) {
            return 0.0;
        }
        // A tiny cost model picks the cheaper strategy from cheap input features
        // (see [`prefer_indexed`]): brute-with-prune for separated/small operands,
        // an STR-tree for large operands whose boxes overlap (prune is weak).
        let min_sq = if prefer_indexed(self.sa.len(), self.sb.len(), box_gap_sq(ba, bb) == 0.0) {
            // Index the larger segment set, iterate the smaller (fewer queries).
            if self.sa.len() <= self.sb.len() {
                indexed_min_sq(&self.sa, &self.sb)
            } else {
                indexed_min_sq(&self.sb, &self.sa)
            }
        } else {
            brute_min_sq(&self.sa, bb, &self.sb, &mut self.sb_boxes)
        };
        min_sq.sqrt()
    }

    /// Whether `a` and `b` are within `d` (`ST_DWithin`), reusing this
    /// workspace's buffers. Equivalent to [`euclidean_within`].
    #[must_use]
    pub fn within(&mut self, a: &Geometry, b: &Geometry, d: f64) -> bool {
        if d < 0.0 {
            return false;
        }
        // Convex operands: resolve the exact gap (O(n+m)) and compare it.
        if let Some(gap) = convex_gap(a, b) {
            return gap <= d;
        }
        self.sa.clear();
        self.sb.clear();
        push_segments(a, &mut self.sa);
        push_segments(b, &mut self.sb);
        if self.sa.is_empty() || self.sb.is_empty() {
            return false;
        }
        let (ba, bb) = (bbox_of(&self.sa), bbox_of(&self.sb));
        if boxes_overlap(ba, bb) && a.intersects(b) {
            return true;
        }
        within_brute(&self.sa, bb, &self.sb, &mut self.sb_boxes, d * d)
    }
}

/// Cost model: should the kernel use the spatial index instead of the brute scan?
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

/// Brute-force squared min gap. `sb_boxes` is reused scratch for B's segment
/// boxes; it is rebuilt each call (cleared, then filled).
fn brute_min_sq(
    sa: &[(Coord, Coord)],
    bb: (Coord, Coord),
    sb: &[(Coord, Coord)],
    sb_boxes: &mut Vec<(Coord, Coord)>,
) -> f64 {
    // Small inputs: the per-pair bbox bookkeeping (building B's boxes, two
    // `box_gap_sq` tests per pair) costs more than the few segment gaps it would
    // save, so run a tight unpruned scan. The threshold is machine-tunable.
    if sa.len() * sb.len() <= crate::config::config().cost_model.min_distance.prune_pair_threshold {
        return brute_min_sq_unpruned(sa, sb);
    }
    sb_boxes.clear();
    sb_boxes.extend(sb.iter().map(|&(p, q)| seg_box(p, q)));
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

/// Tight unpruned `O(n·m)` squared min gap, for small operand pairs.
fn brute_min_sq_unpruned(sa: &[(Coord, Coord)], sb: &[(Coord, Coord)]) -> f64 {
    let mut min_sq = f64::INFINITY;
    for &(a1, a2) in sa {
        for &(b1, b2) in sb {
            min_sq = min_sq.min(segment_gap_sq(a1, a2, b1, b2));
            if min_sq == 0.0 {
                return 0.0;
            }
        }
    }
    min_sq
}

/// Brute `ST_DWithin`: early-exit the moment any pair is within `d` (squared
/// `d_sq`). Mirrors [`brute_min_sq`]'s skip-prune/prune split.
fn within_brute(
    sa: &[(Coord, Coord)],
    bb: (Coord, Coord),
    sb: &[(Coord, Coord)],
    sb_boxes: &mut Vec<(Coord, Coord)>,
    d_sq: f64,
) -> bool {
    if sa.len() * sb.len() <= crate::config::config().cost_model.min_distance.prune_pair_threshold {
        for &(a1, a2) in sa {
            for &(b1, b2) in sb {
                if segment_gap_sq(a1, a2, b1, b2) <= d_sq {
                    return true;
                }
            }
        }
        return false;
    }
    sb_boxes.clear();
    sb_boxes.extend(sb.iter().map(|&(p, q)| seg_box(p, q)));
    for &(a1, a2) in sa {
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

// ---- convex fast path (O(n+m) gap via GJK) ---------------------------------

/// Resolve the minimum distance between `a` and `b` via the convex fast path, or
/// `None` if it doesn't apply — either operand non-convex, or GJK failed to
/// converge — in which case the caller falls back to the exact segment scan.
///
/// Both operands convex ⇒ the gap is found in `O(n+m)` with no segment
/// decomposition: a disjoint bounding-box pair is certainly separated, while an
/// overlapping pair needs the exact [`Intersects`] test to rule out a true
/// intersection (gap `0`) before GJK runs.
fn convex_gap(a: &Geometry, b: &Geometry) -> Option<f64> {
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
            // Drop the repeated closing vertex so the loop is open.
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
    is_convex(shape.verts()).then_some(shape)
}

/// Whether `verts` (an open polygon loop) is a non-degenerate convex polygon:
/// every consecutive turn has the same sign — allowing collinear (zero) turns —
/// with at least one strict turn (so a collinear, zero-area loop is rejected).
/// For a simple polygon this guarantees the vertex set's hull equals the shape,
/// which is what the GJK fast path needs.
fn is_convex(verts: &[Coord]) -> bool {
    let n = verts.len();
    if n < 3 {
        return false;
    }
    let mut sign = 0i8;
    for i in 0..n {
        let turn = orient(verts[i], verts[(i + 1) % n], verts[(i + 2) % n]);
        if turn > 0.0 {
            if sign < 0 {
                return false;
            }
            sign = 1;
        } else if turn < 0.0 {
            if sign > 0 {
                return false;
            }
            sign = -1;
        }
    }
    sign != 0
}

/// Minimum distance between two **disjoint** convex vertex loops via the
/// 2-simplex GJK distance algorithm on their Minkowski difference `A ⊖ B`:
/// the gap equals the distance from the origin to that convex set. `O(n+m)`.
///
/// Returns `None` if it fails to converge within a generous iteration budget, so
/// the caller can fall back to the exact brute/index scan — correctness never
/// depends on GJK converging.
fn convex_min_distance(va: &[Coord], vb: &[Coord]) -> Option<f64> {
    // Support vertex of the Minkowski difference in direction `d`.
    let mink = |d: Coord| support(va, d) - support(vb, -d);

    let mut dir = vb[0] - va[0];
    if dir.x == 0.0 && dir.y == 0.0 {
        dir = Coord::new(1.0, 0.0);
    }
    // Maintain a 2-point simplex (a segment of the Minkowski difference) and walk
    // it toward the origin, replacing the farther endpoint with a better support.
    let mut a = mink(dir);
    let mut b = mink(-dir);
    let budget = 2 * (va.len() + vb.len()) + 8;
    for _ in 0..budget {
        let p = closest_on_segment_to_origin(a, b);
        let p_len2 = p.dot(p);
        if p_len2 <= 1e-24 {
            return Some(0.0); // origin lies on the simplex → shapes touch
        }
        let p_len = p_len2.sqrt();
        let unit = p * (-1.0 / p_len); // unit direction from P toward the origin
        let c = mink(unit);
        // Progress = how far the new support reaches past P toward the origin.
        // `p.dot(unit)` is `-p_len`, so this is `c.dot(unit) + p_len ≥ 0`.
        let progress = c.dot(unit) - p.dot(unit);
        if progress < 1e-12 * (1.0 + p_len) {
            return Some(p_len);
        }
        // Keep the endpoint nearer the origin; replace the other with `c`.
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
        return a; // degenerate segment
    }
    // t = (O - a)·ab / |ab|², clamped to the segment.
    let t = ((-a).dot(ab) / len2).clamp(0.0, 1.0);
    a + ab * t
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
