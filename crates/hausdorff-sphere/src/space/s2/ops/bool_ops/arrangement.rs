//! The arrangement: one shared vertex pool, and every arc cut against every arc
//! of the other operand.
//!
//! This is the step that makes the rest of the algorithm structural rather than
//! numeric. After noding, no arc crosses another arc's interior — every meeting
//! is at a shared pool vertex — so classification becomes a point test and
//! assembly becomes a graph walk. Neither of those needs another geometric
//! predicate, which is why the delicate arithmetic is confined here.
//!
//! **Why a pool rather than raw coordinates.** Tracing matches an arc's end to
//! the next arc's start, and an intersection point computed from edge *(a, b)*
//! against *(c, d)* is not bit-identical to the same point computed from *(c, d)*
//! against *(a, b)* — the cross products are evaluated in a different order. Raw
//! `f64` comparison would see two points a few ulp apart and the walk would dead
//! end. Interning every point into a pool with a snap radius gives exact
//! identity, so a shared boundary between the two operands really is shared.
//!
//! The pool lookup is linear, so noding is `O((n·m) + k²)` in the arc counts and
//! the number of intersections. That is deliberate for a first implementation and
//! is the obvious thing to index later — the cell machinery in this crate is
//! built for exactly that prune. It is called out rather than hidden because a
//! caller feeding two 100 000-vertex coverings should know what they are buying.

use crate::space::s2::{Crossing, S2Point, crossing_sign, intersection};

/// Points within this angular distance are the same point.
///
/// Chosen as a compromise between the two ways this fails. Too small and an
/// intersection computed from each side of a shared edge interns as two vertices,
/// breaking the walk. Too large and genuinely distinct vertices of a fine mesh
/// collapse, welding geometry that should stay separate.
///
/// `1e-12` radians is about 6 µm on Earth — far below any real-world vertex
/// spacing, and comfortably above the error of the non-extended-precision
/// `intersection` construction except for edges meeting at a grazing angle, which
/// is the documented weak spot of that routine and is handled by snapping to an
/// existing vertex when one is near.
pub(super) const SNAP_RADIUS: f64 = 1e-12;

/// Which operand an arc came from. Classification asks a different question of
/// each, so the provenance has to survive noding.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Side {
    A,
    B,
}

/// A directed arc of the arrangement, with its endpoints as pool indices.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) struct Arc {
    pub start: usize,
    pub end: usize,
    pub side: Side,
}

/// Interns points so that geometric coincidence becomes index equality.
pub(super) struct VertexPool {
    points: Vec<S2Point>,
}

impl VertexPool {
    pub(super) fn new() -> Self {
        Self { points: Vec::new() }
    }

    /// The index of `p`, reusing an existing vertex within [`SNAP_RADIUS`].
    ///
    /// Returning the *existing* point rather than averaging is what makes this a
    /// snap: a vertex of the input is never moved by an intersection landing near
    /// it, so the output's shared edges are the input's edges exactly.
    pub(super) fn intern(&mut self, p: S2Point) -> usize {
        for (i, q) in self.points.iter().enumerate() {
            if q.distance(p) <= SNAP_RADIUS {
                return i;
            }
        }
        self.points.push(p);
        self.points.len() - 1
    }

    pub(super) fn point(&self, i: usize) -> S2Point {
        self.points[i]
    }
}

/// One operand's boundary as a list of directed edges, already oriented so the
/// region lies to the left.
pub(super) type Boundary = Vec<(S2Point, S2Point)>;

/// Cut every edge of `a` against every edge of `b` and vice versa, returning the
/// noded arcs over a shared pool.
///
/// Both operands are cut, not just one. Cutting only `a` would leave `b`'s edges
/// crossing `a`'s interior, and those arcs are needed for every operation whose
/// result boundary includes part of `b`.
pub(super) fn node(a: &Boundary, b: &Boundary) -> (VertexPool, Vec<Arc>) {
    let mut pool = VertexPool::new();
    let mut arcs = Vec::new();
    split_against(&mut pool, &mut arcs, a, b, Side::A);
    split_against(&mut pool, &mut arcs, b, a, Side::B);
    (pool, arcs)
}

/// Split each edge of `edges` at its crossings with `others`, appending the
/// resulting arcs.
fn split_against(
    pool: &mut VertexPool,
    arcs: &mut Vec<Arc>,
    edges: &Boundary,
    others: &Boundary,
    side: Side,
) {
    for &(from, to) in edges {
        // Degenerate edge: a zero-length arc bounds nothing and would intern to a
        // self-loop that the walk could never leave.
        if from.distance(to) <= SNAP_RADIUS {
            continue;
        }

        // Collect the split points along this edge, as (distance from `from`,
        // point). The endpoints are always splits.
        let mut cuts: Vec<(f64, S2Point)> = vec![(0.0, from), (from.distance(to), to)];

        for &(c, d) in others {
            if c.distance(d) <= SNAP_RADIUS {
                continue;
            }
            // A T-junction — a vertex of the other edge sitting on this one — is
            // found by *distance*, independently of what `crossing_sign` makes of
            // the pair. The vertex itself is the split, already an endpoint of the
            // other edge, so interning it is enough and no point is constructed.
            //
            // This check used to live under `crossing_sign`'s `Maybe` arm, which
            // meant it only ran when the three points were *exactly* coplanar, and
            // that is a different and much stricter test than the pool's. The two
            // notions of "same point" then disagreed. A meridian edge at longitude
            // 0 is exactly the `y = 0` plane, because `sin 0° = 0` exactly, so a
            // vertex on it is an exact degeneracy and got its split; the same edge
            // at longitude 90° is not, because `cos 90° = 6.1e-17` rather than
            // zero, so a vertex on *it* came back as a definite `DoNotCross` and
            // got no split — a deviation of about 1e-17 radians, five orders of
            // magnitude inside `SNAP_RADIUS`. That left the T-junction this
            // module's own comment warns about: one side of the junction had a
            // vertex and the other did not, the traced walk ran into a vertex it
            // could not leave, and `trace` dropped *every* ring. A 90°-wide lune
            // unioned with a polar cap whose vertices lay on the lune's meridians
            // returned nothing at all — both operands silently lost.
            //
            // Deciding it on `distance_to_arc` puts the split on the same footing
            // as the pool's `intern`: within `SNAP_RADIUS` is the same point, at
            // every longitude. The distance is clamped at the endpoints, so a
            // candidate past either end is measured to that end and only registers
            // when it effectively *is* that end — in which case the cut is a
            // duplicate and the interning below collapses it.
            for candidate in [c, d] {
                if crate::space::s2::arc::distance_to_arc(candidate, from, to).radians()
                    <= SNAP_RADIUS
                {
                    cuts.push((from.distance(candidate), candidate));
                }
            }

            if matches!(crossing_sign(from, to, c, d), Crossing::Cross) {
                let x = intersection(from, to, c, d);
                cuts.push((from.distance(x), x));
            }
        }

        // Sorting by angular distance from `from` is a valid total order because
        // an arc shorter than π is monotone in that distance. S2 edges are
        // required to be shorter than π, so this holds for every input.
        cuts.sort_by(|x, y| x.0.total_cmp(&y.0));

        let mut previous: Option<usize> = None;
        for (_, p) in cuts {
            let index = pool.intern(p);
            match previous {
                // Consecutive cuts that interned to the same vertex are a
                // duplicate split, not a zero-length arc.
                Some(prev) if prev == index => {}
                Some(prev) => {
                    arcs.push(Arc {
                        start: prev,
                        end: index,
                        side,
                    });
                    previous = Some(index);
                }
                None => previous = Some(index),
            }
        }
    }
}
