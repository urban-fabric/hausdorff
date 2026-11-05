//! Boolean set operations on spherical polygons — union, intersection,
//! difference, and symmetric difference, computed natively on the sphere.
//!
//! This is the operation the planar engine cannot do honestly near the seams.
//! `hausdorff-planar`'s boolean ops are correct in the plane, so applying them to
//! lon/lat treats the antimeridian as a wall and the poles as a line: a polygon
//! spanning either is cut or inverted, and the failure is silent. Projecting first
//! moves the problem rather than removing it, because every projection is wrong
//! about something. Working in unit vectors means there is no seam to be wrong
//! about — a great-circle arc is a plane through the origin and nothing
//! distinguishes ±180°.
//!
//! ## Why not a sweep
//!
//! The planar engine is a Martinez–Rueda sweep, and a sweep needs a **total order**
//! to advance along. The sphere has none: sweeping by longitude reintroduces
//! exactly the antimeridian and pole degeneracies this module exists to avoid, and
//! sweeping by any other parameterization has the same problem somewhere else. A
//! sphere has no global monotone direction.
//!
//! So the architecture here is the *arrangement* approach, which needs no order at
//! all — and it is the one `hausdorff-planar`'s `make_valid` already uses for the
//! same reason:
//!
//! 1. **node** — cut every arc against every arc of the other operand, interning
//!    endpoints into one shared vertex pool, so no arc crosses another's interior;
//! 2. **classify** — keep the arcs the result's boundary is made of, by testing
//!    each arc's midpoint against the other operand;
//! 3. **trace** — walk the survivors into closed rings;
//! 4. **nest** — sort rings into shells and holes by orientation and containment.
//!
//! Each step is exact or structural except the first, which confines all the
//! delicate arithmetic to one place. Steps 2–4 use no geometric predicate beyond a
//! point-in-polygon test.
//!
//! ## What it is built on
//!
//! Nothing here invents a primitive. The robust `crossing_sign` and `robust_sign`
//! predicates, the `intersection` construction, and the proof-backed
//! `Loop::contains_point` were already in this crate; this module is the
//! orchestration that was missing. That matters for trust: the parts that decide
//! *signs* are the parts that were already validated.
//!
//! ## Scope, stated honestly
//!
//! - **Input loops must be simple** — no self-intersections. That is `Loop`'s
//!   existing precondition, not a new one. Self-intersecting input produces a best
//!   effort rather than a panic: a ring that cannot be closed is dropped.
//! - **Arcs must be shorter than π**, which S2 requires anyway. The arrangement
//!   relies on it to sort cuts along an edge by angular distance.
//! - **Noding is `O(n·m)`.** No spatial index is applied yet, so two large
//!   coverings are quadratic. The cell machinery in this crate is the obvious
//!   prune and is a performance change, not a correctness one.
//! - **Results exactly on a boundary** follow `Loop::contains_point`, which does
//!   not guarantee a verdict there. Coincident boundaries and T-junctions are both
//!   handled structurally instead — a shared stretch of boundary by the duplicate
//!   arc table, and a vertex lying on the other operand's edge by splitting that
//!   edge at it, decided on distance so that the snap radius means the same thing
//!   in noding as it does in the vertex pool.
//! - **An operand larger than a hemisphere is only partly supported.** Union and
//!   difference are correct where the result is simply connected — their
//!   orientation is read off the `Polygon` type and a single-ring result is taken
//!   as a shell, neither of which measures anything.
//!
//!   *Intersection* is not, because classification asks
//!   `Loop::contains_point` which is itself wrong above a hemisphere: it reported a
//!   160°-wide region for a 200°-wide lune. See
//!   `Loop::contains_point_winding` for the measured boundary. A 200° lune
//!   intersected with a 90° lune inside it therefore comes back empty rather than
//!   as the inner lune.
//!
//!   Neither is a result that is over a hemisphere *and* multi-ring, which is the
//!   one case `nest` has to measure a ring to classify it; see the limit note
//!   there. Subtracting an interior quad from a 250° lune returns the complement.
//!   Operands that large are not geographic features, but the limits are real and
//!   are better stated than discovered.
//!
//! ## Oracles
//!
//! The suite checks set algebra in steradians — inclusion–exclusion, idempotence,
//! De Morgan, monotone bounds — exactly as the planar suite does, because those
//! relations tie the four operations to one another and to `area` without needing a
//! reference implementation. It also checks **agreement with the planar engine
//! away from the seams**, where the two must give the same answer, and
//! **disagreement across the antimeridian**, where the planar one is wrong and
//! this one is the reason the module exists.

mod arrangement;
mod assemble;

#[cfg(test)]
mod tests;

use crate::space::s2::{Loop, Polygon, S2Point};
use arrangement::Boundary;

/// Which set operation to compute.
///
/// Symmetric difference is deliberately absent: it is `union(a \ b, b \ a)` by
/// definition and is offered as [`symmetric_difference`] rather than as a fourth
/// orientation rule inside the classifier. See the `assemble` module header.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Operation {
    /// Everything in either operand.
    Union,
    /// Only what is in both.
    Intersection,
    /// What is in `a` and not in `b`.
    Difference,
}

/// Set operations on spherical polygons.
///
/// Implemented for `[Polygon]` rather than `Polygon` because every operation can
/// *produce* more than one component — intersecting two crescents can give two
/// pieces — so a type that holds one polygon could not represent its own results.
/// A single polygon participates as `std::slice::from_ref(&poly)`.
pub trait SphericalBooleanOps {
    /// Everything covered by either operand.
    fn union(&self, other: &Self) -> Vec<Polygon>;
    /// Only what both operands cover.
    fn intersection(&self, other: &Self) -> Vec<Polygon>;
    /// What this operand covers and the other does not.
    fn difference(&self, other: &Self) -> Vec<Polygon>;
    /// What exactly one operand covers.
    fn symmetric_difference(&self, other: &Self) -> Vec<Polygon>;
}

impl SphericalBooleanOps for [Polygon] {
    fn union(&self, other: &Self) -> Vec<Polygon> {
        boolean_op(self, other, Operation::Union)
    }
    fn intersection(&self, other: &Self) -> Vec<Polygon> {
        boolean_op(self, other, Operation::Intersection)
    }
    fn difference(&self, other: &Self) -> Vec<Polygon> {
        boolean_op(self, other, Operation::Difference)
    }
    fn symmetric_difference(&self, other: &Self) -> Vec<Polygon> {
        symmetric_difference(self, other)
    }
}

/// `(a \ b) ∪ (b \ a)` — the definition, rather than a fourth sign convention.
#[must_use]
pub fn symmetric_difference(a: &[Polygon], b: &[Polygon]) -> Vec<Polygon> {
    let left = boolean_op(a, b, Operation::Difference);
    let right = boolean_op(b, a, Operation::Difference);
    boolean_op(&left, &right, Operation::Union)
}

/// Compute `op` over two polygon sets.
#[must_use]
pub fn boolean_op(a: &[Polygon], b: &[Polygon], op: Operation) -> Vec<Polygon> {
    // Empty operands have exact answers, and routing them through the arrangement
    // would produce an empty vertex pool and a vacuous walk.
    if a.is_empty() || b.is_empty() {
        return match op {
            Operation::Union => {
                let mut out = a.to_vec();
                out.extend_from_slice(b);
                out
            }
            Operation::Intersection => Vec::new(),
            Operation::Difference => a.to_vec(),
        };
    }

    // Disjoint bounds short-circuit every operation and are the common case in a
    // spatial join, where most candidate pairs do not meet.
    if bounds_disjoint(a, b) {
        return match op {
            Operation::Union => {
                let mut out = a.to_vec();
                out.extend_from_slice(b);
                out
            }
            Operation::Intersection => Vec::new(),
            Operation::Difference => a.to_vec(),
        };
    }

    let a_boundary = boundary_of(a);
    let b_boundary = boundary_of(b);
    let (pool, arcs) = arrangement::node(&a_boundary, &b_boundary);

    let a_owned = a.to_vec();
    let b_owned = b.to_vec();
    let a_contains = move |p: S2Point| a_owned.iter().any(|poly| poly.contains_point(p));
    let b_contains = move |p: S2Point| b_owned.iter().any(|poly| poly.contains_point(p));

    let kept = assemble::classify(&pool, &arcs, op, &a_contains, &b_contains);
    let rings = assemble::trace(&pool, &kept);
    assemble::nest(rings)
}

/// Whether the two sets' bounding rectangles cannot meet.
///
/// `Rect` offers `intersection` rather than a boolean test, and an empty result
/// is the "cannot meet" answer. Going through the lat/lng bound rather than a
/// spherical cap because the rect bound is already cached per polygon and this is
/// only a prune — a false "they might meet" costs a wasted arrangement, never a
/// wrong answer.
fn bounds_disjoint(a: &[Polygon], b: &[Polygon]) -> bool {
    !a.iter().any(|pa| {
        let ra = pa.rect_bound();
        b.iter()
            .any(|pb| !ra.intersection(&pb.rect_bound()).is_empty())
    })
}

/// Flatten a polygon set to directed edges oriented so the covered region lies to
/// the left of every one.
///
/// Normalizing here rather than trusting the input is what lets the classifier be
/// a pure point test. A shell wound clockwise or a hole wound counter-clockwise
/// would otherwise invert the region silently, and the result would be a plausible
/// polygon of the wrong shape — the hardest kind of error to notice.
fn boundary_of(polys: &[Polygon]) -> Boundary {
    let mut edges = Vec::new();
    for poly in polys {
        push_ring(&mut edges, poly.shell(), true);
        for hole in poly.holes() {
            push_ring(&mut edges, hole, false);
        }
    }
    edges
}

/// Append a ring's edges, oriented so the polygon's interior is on the left.
///
/// **The direction is decided by the ring's role, not by measuring its winding.**
/// A shell is kept as given and a hole is reversed, unconditionally. That is not a
/// shortcut — measuring is actively wrong, and this is where it was wrong.
///
/// `Polygon` already defines its interior as the region to the *left* of the
/// shell's directed edges: that is what `Loop::area` returns and what
/// `Loop::contains_point` tests. So a shell is already region-on-left by
/// definition, and a hole — whose own left-region is the hole, which the polygon
/// excludes — needs exactly one reversal.
///
/// The previous version instead inferred winding from the sign of
/// `polygon_signed_area_vos`, and that sign is **ambiguous for any region larger
/// than a hemisphere**. The van Oosterom sum returns the spherical excess, so a
/// counter-clockwise loop enclosing `A > 2π` reports `A − 4π`, which is negative —
/// indistinguishable from a small clockwise loop. A 200°-wide lune (area 6.98 sr)
/// was therefore read as clockwise, reversed, and silently inverted into its own
/// complement. Its union with anything then came back *smaller* than itself.
///
/// There is no winding sign to recover here: a simple closed curve on a sphere
/// divides it into two regions and neither is intrinsically "the inside". Which
/// one counts is a decision the type makes, and reading it off the type is the
/// only correct way to get it.
fn push_ring(edges: &mut Boundary, ring: &Loop, is_shell: bool) {
    let mut vertices = ring.vertices().to_vec();
    if vertices.len() < 3 {
        return;
    }
    if !is_shell {
        vertices.reverse();
    }
    for i in 0..vertices.len() {
        edges.push((vertices[i], vertices[(i + 1) % vertices.len()]));
    }
}

/// Total area of a polygon set in steradians, holes subtracted.
///
/// Here rather than left to the caller because every set-algebra assertion needs
/// it and summing `Polygon::area` over a `Vec` is the kind of one-liner that gets
/// written four slightly different ways.
#[must_use]
pub fn total_area(polys: &[Polygon]) -> f64 {
    polys.iter().map(Polygon::area).sum()
}
