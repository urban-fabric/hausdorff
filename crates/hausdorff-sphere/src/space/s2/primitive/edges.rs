//! `s2` great-circle edge crossings — do two spherical edges cross, and where.
//! Robust (built on the orientation [`robust_sign`]) and
//! projection-free. Ported from golang/geo's `edge_crossings.go`. These are the
//! foundation for spherical clipping / polygon boolean operations.

use super::predicates::robust_sign;
use crate::space::s2::S2Point;

/// The result of testing whether two edges cross.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Crossing {
    /// The edges cross at a point interior to both.
    Cross,
    /// The edges share or touch a vertex (a boundary case, not an interior cross).
    Maybe,
    /// The edges do not cross.
    DoNotCross,
}

/// Whether directed edge `a → b` crosses directed edge `c → d`.
///
/// Returns [`Crossing::Cross`] only for a crossing interior to *both* edges,
/// [`Crossing::Maybe`] when an endpoint lies on the other edge's great circle
/// (shared/touching vertex), and [`Crossing::DoNotCross`] otherwise. Symmetric
/// in the natural ways: `crossing_sign(b,a,c,d) == crossing_sign(a,b,c,d)` and
/// `crossing_sign(c,d,a,b) == crossing_sign(a,b,c,d)`.
#[must_use]
pub fn crossing_sign(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Crossing {
    let abc = robust_sign(a, b, c);
    let abd = robust_sign(a, b, d);
    if abc == 0 || abd == 0 {
        return Crossing::Maybe;
    }
    if abc == abd {
        return Crossing::DoNotCross;
    }
    let cda = robust_sign(c, d, a);
    let cdb = robust_sign(c, d, b);
    if cda == 0 || cdb == 0 {
        return Crossing::Maybe;
    }
    if cda == cdb {
        return Crossing::DoNotCross;
    }
    // The orientation signs now say each arc separates the other's endpoints —
    // necessary but NOT sufficient. The two great circles meet at an antipodal
    // pair `±x`; arc `a→b` may cross the `c→d` circle at `x` while arc `c→d`
    // crosses the `a→b` circle at `−x` (e.g. two small polygons on opposite sides
    // of the globe). That is not a shared point. Confirm one of `±x` lies on
    // *both* minor arcs.
    if arcs_share_a_crossing(a, b, c, d) {
        Crossing::Cross
    } else {
        Crossing::DoNotCross
    }
}

/// Whether the great circles of `a→b` and `c→d` meet at a point lying on **both**
/// minor arcs. The intersection directions are `±(a×b)×(c×d)`; a normalized
/// direction `p` is within the minor arc `a→b` iff it is no farther from either
/// endpoint than they are from each other (`p·a ≥ a·b` and `p·b ≥ a·b`).
fn arcs_share_a_crossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> bool {
    let n1 = a.vector().cross(b.vector());
    let n2 = c.vector().cross(d.vector());
    let x = n1.cross(n2);
    if x.dot(x) == 0.0 {
        return false; // coincident great circles: no isolated crossing
    }
    let x = x.normalize();
    let on_both =
        |p: crate::space::s2::S2Point| within_minor_arc(a, b, p) && within_minor_arc(c, d, p);
    on_both(S2Point::from_vector(x)) || on_both(S2Point::from_vector(x.mul(-1.0)))
}

/// Whether unit point `p` (already on the great circle through `a`, `b`) lies on
/// the minor arc `a → b`, endpoints included.
fn within_minor_arc(a: S2Point, b: S2Point, p: S2Point) -> bool {
    let (av, bv, pv) = (a.vector(), b.vector(), p.vector());
    let ab = av.dot(bv);
    pv.dot(av) >= ab && pv.dot(bv) >= ab
}

/// The intersection point of two crossing great-circle edges `a0 → a1` and
/// `b0 → b1`. Only meaningful when [`crossing_sign`] is [`Crossing::Cross`].
///
/// Computed as the normalized cross product of the two edge-plane normals,
/// flipped to the hemisphere of the edges. This is the direct (non
/// extended-precision) construction — accurate except when the two edges meet
/// at an extremely small angle.
#[must_use]
pub fn intersection(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point) -> S2Point {
    let n1 = a0.vector().cross(a1.vector());
    let n2 = b0.vector().cross(b1.vector());
    let mut x = n1.cross(n2).normalize();
    let sum = a0
        .vector()
        .add(a1.vector())
        .add(b0.vector())
        .add(b1.vector());
    if x.dot(sum) < 0.0 {
        x = x.mul(-1.0);
    }
    S2Point::from_vector(x)
}
