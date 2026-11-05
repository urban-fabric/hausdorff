//! Topology-preserving simplification — Douglas–Peucker that refuses a removal
//! which would change what the geometry *is*.
//!
//! Plain [`Simplify`](crate::space::r2::Simplify) is a per-ring operation: it picks the
//! vertices that approximate a chord well and discards the rest, with no notion
//! of the geometry those vertices belong to. On a single open line that is exactly
//! right. On polygonal data it is quietly destructive in two ways:
//!
//! - **a ring can simplify into itself.** Removing enough of a concave bay lets
//!   the replacement chord cross another part of the same ring, and the result is
//!   a self-intersecting polygon — invalid, and rejected by anything downstream
//!   that validates.
//! - **a ring can collapse or invert.** A thin sliver or a small island can lose
//!   enough vertices that its signed area flips sign or goes to zero, so the
//!   polygon turns inside out or vanishes.
//!
//! Both matter most for exactly the data people most want to simplify: census
//! tracts, admin boundaries, parcel layers. Those are *coverages* — adjacent
//! polygons sharing edges — and simplifying each ring independently is what opens
//! gaps and overlaps along the shared boundaries.
//!
//! ## What "topology-preserving" does and does not promise here
//!
//! It promises, per ring: **the simplified ring is simple** (no self-intersection)
//! and **its orientation and non-degeneracy are preserved**. That is the guarantee
//! JTS's `TopologyPreservingSimplifier` makes, and it is what makes the output
//! valid rather than merely smaller.
//!
//! It does **not** promise that two rings which shared an edge still share it.
//! Preserving a *coverage* is a different and stronger operation — it needs all
//! the rings at once, so that a shared boundary is simplified once and the result
//! substituted into both neighbours. That is `CoverageSimplify` in recent JTS and
//! PostGIS, and it is not this. Stating the difference matters because the names
//! are close enough to assume: if you are simplifying an admin layer and need the
//! borders to still match, this is not sufficient on its own.
//!
//! ## How the guarantee is enforced
//!
//! The structure is ordinary Douglas–Peucker with a veto. Where plain RDP accepts
//! a chord as soon as every dropped vertex is within `epsilon`, this additionally
//! asks whether the chord crosses any part of the ring it is not replacing. If it
//! does, the removal is refused and the recursion splits anyway — so the output is
//! never worse than the input topologically, at the cost of keeping more vertices
//! than `epsilon` alone would.
//!
//! The check is against the *original* ring rather than the partially-simplified
//! one. That is deliberate and it is the conservative choice: the original is a
//! superset of every chord the final result could contain, so a chord that clears
//! it cannot conflict with whatever the rest of the recursion decides to keep.
//! Checking against a moving target would make the result depend on traversal
//! order.

use crate::space::r2::{Coord, Geometry, LineString, MultiLineString, MultiPolygon, Polygon};

/// Simplify while keeping each ring simple and correctly oriented.
pub trait SimplifyTopologyPreserving {
    /// Douglas–Peucker at tolerance `epsilon`, refusing any vertex removal that
    /// would make a ring self-intersect, collapse, or invert.
    ///
    /// Always returns geometry at least as topologically sound as the input, which
    /// plain `simplify` does not.
    fn simplify_topology_preserving(&self, epsilon: f64) -> Self;
}

/// Whether the chord `a`–`b` crosses any segment of `ring` other than the run it
/// would replace.
///
/// `skip` is the half-open index range of the run being collapsed; its own
/// segments are excluded because a chord is of course allowed to cross the bay it
/// is cutting off.
fn chord_is_clear(ring: &[Coord], a: Coord, b: Coord, skip: std::ops::Range<usize>) -> bool {
    let n = ring.len();
    for i in 0..n.saturating_sub(1) {
        // The run being replaced, plus the segments immediately adjoining the
        // chord's own endpoints — those share a vertex with it, which is a touch
        // and not a crossing.
        if skip.contains(&i) {
            continue;
        }
        let (c, d) = (ring[i], ring[i + 1]);
        if c == a || c == b || d == a || d == b {
            continue;
        }
        if segments_cross(a, b, c, d) {
            return false;
        }
    }
    true
}

/// A proper crossing test in `f64`, falling back on the sign of the orientation
/// determinants.
///
/// Deliberately *not* the exact `i128` predicate. Those require coordinates
/// already on an integer grid, and simplification is applied to arbitrary `f64`
/// data — scaling it onto a grid to decide a veto would change which vertices
/// survive depending on the grid chosen. A false "they cross" here is harmless:
/// it keeps a vertex that could have been dropped, which is the safe direction for
/// an operation whose contract is "never make topology worse".
fn segments_cross(a: Coord, b: Coord, c: Coord, d: Coord) -> bool {
    let orient = |p: Coord, q: Coord, r: Coord| {
        let v = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
        if v > 0.0 {
            1
        } else if v < 0.0 {
            -1
        } else {
            0
        }
    };
    let (d1, d2) = (orient(a, b, c), orient(a, b, d));
    let (d3, d4) = (orient(c, d, a), orient(c, d, b));
    // Strict crossing only: each segment's endpoints must lie strictly on
    // opposite sides of the other. Collinear overlap and endpoint touching are
    // not crossings, and treating them as such would veto almost every removal
    // on a ring with collinear runs.
    d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0 && d1 != d2 && d3 != d4
}

/// Twice the signed area of a ring — the sign is the winding.
///
/// Shifted to a local origin before summing, for the conditioning reason given on
/// [`util::signed_area2`](crate::space::r2::util::signed_area2): the untranslated form loses
/// the whole answer to cancellation once coordinates pass about `1e8`. This one is
/// compared before and after a simplification, so a wrong magnitude silently
/// changes which removals are accepted.
fn signed_area2(ring: &[Coord]) -> f64 {
    let Some(o) = ring.first().copied() else {
        return 0.0;
    };
    let term = |a: Coord, b: Coord| (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
    let mut s = 0.0;
    for w in ring.windows(2) {
        s += term(w[0], w[1]);
    }
    if ring.len() >= 2 && ring[0] != ring[ring.len() - 1] {
        s += term(ring[ring.len() - 1], ring[0]);
    }
    s
}

/// Douglas–Peucker with the crossing veto.
///
/// `offset` is where this run starts in the original ring, so the veto knows which
/// segments it is allowed to cross.
fn rdp_guarded(
    original: &[Coord],
    offset: usize,
    run: &[Coord],
    epsilon: f64,
    closed: bool,
) -> Vec<Coord> {
    if run.len() < 3 {
        return run.to_vec();
    }
    let (first, last) = (run[0], run[run.len() - 1]);

    let mut max_dist = 0.0;
    let mut index = 0;
    for (i, &c) in run.iter().enumerate().take(run.len() - 1).skip(1) {
        let d = point_segment_distance(c, first, last);
        if d > max_dist {
            max_dist = d;
            index = i;
        }
    }

    // The two reasons to refuse a collapse that plain RDP would accept.
    let within_tolerance = max_dist <= epsilon;
    let topologically_safe =
        !closed || chord_is_clear(original, first, last, offset..(offset + run.len() - 1));

    if within_tolerance && topologically_safe {
        return vec![first, last];
    }
    // If the chord is within tolerance but unsafe, there is no split point chosen
    // by distance — every vertex is inside epsilon. Split at the midpoint so the
    // recursion still makes progress and the veto cannot loop.
    if index == 0 {
        index = run.len() / 2;
    }

    let mut left = rdp_guarded(original, offset, &run[..=index], epsilon, closed);
    let right = rdp_guarded(original, offset + index, &run[index..], epsilon, closed);
    left.pop();
    left.extend(right);
    left
}

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

/// Simplify one closed ring, keeping it simple, non-degenerate, and correctly
/// wound.
///
/// A ring is never reduced below four coordinates (three distinct vertices plus
/// the repeat): below that it bounds no area, and a "simplified" ring with no area
/// is a deleted feature rather than a smaller one.
fn simplify_ring(ring: &[Coord], epsilon: f64) -> Vec<Coord> {
    if ring.len() <= 4 {
        return ring.to_vec();
    }
    let before = signed_area2(ring);
    let mut out = rdp_guarded(ring, 0, ring, epsilon, true);
    if out.first() != out.last()
        && let Some(&f) = out.first()
    {
        out.push(f);
    }
    if out.len() < 4 {
        return ring.to_vec();
    }
    let after = signed_area2(&out);
    // An inverted or collapsed ring is worse than an unsimplified one, so the
    // original is returned rather than a wrong answer. This is the backstop behind
    // the per-chord veto: the veto prevents self-intersection, and this catches the
    // remaining way a ring can go wrong.
    if after == 0.0 || (after > 0.0) != (before > 0.0) {
        return ring.to_vec();
    }
    out
}

/// Simplify an open line. There is no ring to self-intersect against, so this is
/// plain Douglas–Peucker — stated rather than left implicit, because a caller
/// reaching for the topology-preserving variant on a linestring should know it
/// gets no extra guarantee.
fn simplify_open(coords: &[Coord], epsilon: f64) -> Vec<Coord> {
    rdp_guarded(coords, 0, coords, epsilon, false)
}

impl SimplifyTopologyPreserving for LineString {
    fn simplify_topology_preserving(&self, epsilon: f64) -> Self {
        let closed = self.0.len() >= 4 && self.0[0] == self.0[self.0.len() - 1];
        if closed {
            LineString::new(simplify_ring(&self.0, epsilon))
        } else {
            LineString::new(simplify_open(&self.0, epsilon))
        }
    }
}

impl SimplifyTopologyPreserving for MultiLineString {
    fn simplify_topology_preserving(&self, epsilon: f64) -> Self {
        MultiLineString::new(
            self.iter()
                .map(|ls| ls.simplify_topology_preserving(epsilon))
                .collect(),
        )
    }
}

impl SimplifyTopologyPreserving for Polygon {
    fn simplify_topology_preserving(&self, epsilon: f64) -> Self {
        Polygon::new(
            LineString::new(simplify_ring(&self.exterior().0, epsilon)),
            self.interiors()
                .iter()
                .map(|h| LineString::new(simplify_ring(&h.0, epsilon)))
                .collect(),
        )
    }
}

impl SimplifyTopologyPreserving for MultiPolygon {
    fn simplify_topology_preserving(&self, epsilon: f64) -> Self {
        MultiPolygon::new(
            self.iter()
                .map(|p| p.simplify_topology_preserving(epsilon))
                .collect(),
        )
    }
}

impl SimplifyTopologyPreserving for Geometry {
    fn simplify_topology_preserving(&self, epsilon: f64) -> Self {
        match self {
            Geometry::LineString(ls) => {
                Geometry::LineString(ls.simplify_topology_preserving(epsilon))
            }
            Geometry::MultiLineString(mls) => {
                Geometry::MultiLineString(mls.simplify_topology_preserving(epsilon))
            }
            Geometry::Polygon(p) => Geometry::Polygon(p.simplify_topology_preserving(epsilon)),
            Geometry::MultiPolygon(mp) => {
                Geometry::MultiPolygon(mp.simplify_topology_preserving(epsilon))
            }
            Geometry::GeometryCollection(gs) => Geometry::GeometryCollection(
                gs.iter()
                    .map(|g| g.simplify_topology_preserving(epsilon))
                    .collect(),
            ),
            // Points, lines, rects, and triangles have no vertices to spare: a
            // point cannot be simplified, and the other three are already their
            // own minimal description.
            other => other.clone(),
        }
    }
}
