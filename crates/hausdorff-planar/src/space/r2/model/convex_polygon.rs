//! `r2::ConvexPolygon` — a planar polygon **validated convex once**, after which
//! every convex-only fast algorithm is *total*: no per-call convexity check, no
//! `Option`. The parse-don't-validate companion to [`Polygon`]'s runtime
//! [`is_convex`](Polygon::is_convex): a general polygon that *happens* to be
//! convex still gets the fast paths transparently (e.g. the GJK
//! [`euclidean_distance`](crate::space::r2::euclidean_distance) path), while a `ConvexPolygon` makes the
//! invariant a compile-time guarantee and unlocks algorithms that are only valid
//! (or only fast) for convex regions.
//!
//! This is the planar mirror of `s2::ConvexLoop`: the *same* generic
//! [`Orientation`](crate::algebra::Orientation)-based machinery —
//! [`convex_winding`] to validate, [`contains_convex`] for the O(log n)
//! containment — instantiated at [`Coord`] instead of `S2Point`. The convexity
//! and containment laws are machine-checked once in `riemann`
//! (`Geometry/Orientation.OrientGeom.leftOfAll_convex`, planar
//! `ExactArith/ConvexPolygon.leftOfAllStrict_convex`).

use hausdorff_traits::{contains_convex, convex_winding};

use crate::space::r2::util::{orient, signed_area2};
use crate::space::r2::{Coord, LineString, Polygon};

/// A planar polygon known to be convex and wound counter-clockwise, stored as an
/// open ring (no repeated closing vertex).
///
/// Construct with [`ConvexPolygon::from_vertices`] / [`from_polygon`](Self::from_polygon),
/// which validate convexity and normalize the winding. Once built, the convexity
/// invariant holds for the lifetime of the value.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvexPolygon {
    /// Invariant: a counter-clockwise convex ring with `≥ 3` vertices and no
    /// closing duplicate.
    ring: Vec<Coord>,
}

impl ConvexPolygon {
    /// Validate a vertex ring as a convex polygon, normalizing to
    /// counter-clockwise. A trailing vertex equal to the first (a closed ring) is
    /// dropped. Returns `None` if the ring is non-convex (a reflex vertex) or
    /// degenerate (`< 3` distinct vertices, or all-collinear).
    #[must_use]
    pub fn from_vertices(ring: Vec<Coord>) -> Option<Self> {
        let mut ring = ring;
        if ring.len() >= 2 && ring.first() == ring.last() {
            ring.pop();
        }
        match convex_winding(&ring) {
            Some(s) if s > 0.0 => Some(Self { ring }),
            // Clockwise convex: reverse to the canonical CCW ring (same region).
            Some(_) => {
                ring.reverse();
                debug_assert_eq!(convex_winding(&ring), Some(1.0));
                Some(Self { ring })
            }
            None => None,
        }
    }

    /// Validate a [`Polygon`] as convex: it must have **no holes** (a holed
    /// polygon is never convex) and a convex exterior ring. See
    /// [`from_vertices`](Self::from_vertices).
    #[must_use]
    pub fn from_polygon(p: &Polygon) -> Option<Self> {
        if !p.interiors().is_empty() {
            return None;
        }
        Self::from_vertices(p.exterior().0.clone())
    }

    /// The vertices, counter-clockwise, as an open ring (no closing duplicate).
    #[must_use]
    pub fn vertices(&self) -> &[Coord] {
        &self.ring
    }

    /// The number of vertices (= the number of edges).
    #[must_use]
    pub fn num_vertices(&self) -> usize {
        self.ring.len()
    }

    /// Area — the (positive) shoelace area of the CCW ring.
    #[must_use]
    pub fn area(&self) -> f64 {
        signed_area2(&self.ring).abs() / 2.0
    }

    /// True iff `p` lies in the interior — **O(log n)** via triangle-fan binary
    /// search, the convex-only speedup over the general polygon's O(n) even-odd
    /// ray cast. Computes the same predicate (inside ⇔ left of every edge); see
    /// [`contains_convex`]. Boundary results (exactly on an edge/vertex) are not
    /// guaranteed.
    #[must_use]
    pub fn contains(&self, p: Coord) -> bool {
        contains_convex(&self.ring, p)
    }

    /// True iff every vertex of `other` lies in `self` — exact convex containment.
    /// A convex region is closed under chords (`riemann`
    /// `OrientGeom.leftOfAll_convex`), so it contains another convex polygon iff
    /// it contains that polygon's vertices. `O(m log n)`.
    #[must_use]
    pub fn contains_polygon(&self, other: &ConvexPolygon) -> bool {
        other.ring.iter().all(|&q| self.contains(q))
    }

    /// Back to a general (closed, hole-free) [`Polygon`].
    #[must_use]
    pub fn to_polygon(&self) -> Polygon {
        Polygon::new(LineString::new(self.ring.clone()), Vec::new())
    }

    /// True iff `self` and `other` are **disjoint** (share no interior point) — the
    /// **separating-axis theorem** for convex polygons, an `O((n+m)²)` overlap test
    /// with no edge-crossing scan. Two convex polygons are disjoint iff some edge
    /// of one has the *entire* other polygon strictly on its outer side: that
    /// edge's line is a separating line.
    ///
    /// Convex-only and proof-backed: a convex polygon is the intersection of its
    /// edge half-planes (`riemann` `Geometry/Orientation.OrientGeom.leftOfAll_convex`),
    /// so if every vertex of `other` lies strictly outside one of `self`'s edge
    /// half-planes (`orient < 0` — the right/outer side of a CCW edge), then `self`
    /// (⊆ that half-plane) and `other` cannot meet. A boundary touch counts as
    /// intersecting (the test is strict).
    #[must_use]
    pub fn disjoint(&self, other: &ConvexPolygon) -> bool {
        separating_edge(&self.ring, &other.ring) || separating_edge(&other.ring, &self.ring)
    }

    /// True iff `self` and `other` overlap — the negation of [`disjoint`](Self::disjoint).
    #[must_use]
    pub fn intersects(&self, other: &ConvexPolygon) -> bool {
        !self.disjoint(other)
    }

    /// The **intersection** of two convex polygons — itself convex — by the
    /// Sutherland–Hodgman algorithm: clip `self`'s ring successively against each
    /// edge half-plane of `other`. `None` if the overlap is empty or degenerate
    /// (a point/segment, i.e. `< 3` vertices or zero area). `O(n·m)`.
    ///
    /// A convex-only algorithm: clipping a convex region by a half-plane yields a
    /// convex region (each edge half-plane is convex, and the intersection of
    /// convex sets is convex — `riemann`
    /// `Geometry/Orientation.OrientGeom.leftOfAll_convex`), so the result of
    /// clipping by *every* edge of the convex `other` is `self ∩ other`, still
    /// convex. The crossing point of a subject edge with a clip line is exact
    /// because the orientation is affine along the edge (`riemann`
    /// `ExactArith/OrientationLaws.triDet_lerp`).
    #[must_use]
    pub fn intersection(&self, other: &ConvexPolygon) -> Option<ConvexPolygon> {
        let mut output = self.ring.clone();
        let m = other.ring.len();
        for j in 0..m {
            if output.is_empty() {
                return None;
            }
            output = clip_by_halfplane(&output, other.ring[j], other.ring[(j + 1) % m]);
        }
        dedup_cycle(&mut output);
        ConvexPolygon::from_vertices(output)
    }
}

/// The crossing point of segment `p→q` with the line through `a→b`, assuming `p`
/// and `q` lie on strictly opposite sides. Exact: the orientation is affine along
/// the segment (`triDet_lerp`), so the root `t = d_p / (d_p − d_q)` of the signed
/// distance is the true crossing.
fn line_crossing(p: Coord, q: Coord, a: Coord, b: Coord) -> Coord {
    let dp = orient(a, b, p);
    let dq = orient(a, b, q);
    let t = dp / (dp - dq);
    p + (q - p) * t
}

/// Sutherland–Hodgman clip of the closed `subject` ring to the left half-plane of
/// the directed line `a→b` (the inner side of a CCW edge): keep vertices with
/// `orient(a, b, ·) ≥ 0`, inserting crossing points where edges enter/leave.
fn clip_by_halfplane(subject: &[Coord], a: Coord, b: Coord) -> Vec<Coord> {
    let inside = |p: Coord| orient(a, b, p) >= 0.0;
    let n = subject.len();
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..n {
        let cur = subject[i];
        let nxt = subject[(i + 1) % n];
        let (cur_in, nxt_in) = (inside(cur), inside(nxt));
        if cur_in {
            out.push(cur);
            if !nxt_in {
                out.push(line_crossing(cur, nxt, a, b));
            }
        } else if nxt_in {
            out.push(line_crossing(cur, nxt, a, b));
        }
    }
    out
}

/// SAT half: does some directed edge of the CCW ring `a` put *all* of `b`
/// strictly on its outer (right) side? Such an edge line separates the two convex
/// polygons. `O(|a|·|b|)`.
fn separating_edge(a: &[Coord], b: &[Coord]) -> bool {
    let n = a.len();
    (0..n).any(|i| {
        let (p, q) = (a[i], a[(i + 1) % n]);
        b.iter().all(|&x| orient(p, q, x) < 0.0)
    })
}

/// Drop consecutive duplicate vertices (including the wrap seam) — clipping can
/// emit a crossing point coincident with a kept vertex.
fn dedup_cycle(ring: &mut Vec<Coord>) {
    ring.dedup();
    while ring.len() >= 2 && ring.first() == ring.last() {
        ring.pop();
    }
}

#[cfg(test)]
mod tests;
