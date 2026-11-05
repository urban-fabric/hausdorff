//! Boolean ops on the inputs that break sweeps: self-touching rings, holes meeting
//! their shell, duplicate and collinear vertices, slivers, and operands that are
//! each other's subsets.
//!
//! The Martinez–Rueda sweep is the hardest algorithm in this crate to test, because
//! for a general pair of polygons there is no closed form to compare against — and
//! the suite's existing answer to that (closed-form rectangle overlap, a convex
//! clipper, golden constants) only covers inputs those oracles can express. None of
//! them can say what the union of two rings that touch at a single vertex should be.
//!
//! So the oracle here is **set algebra**, which needs no reference implementation at
//! all. The four operations are not independent, and the relations between them have
//! to hold for every input whatsoever:
//!
//! - `|A ∪ B| = |A| + |B| − |A ∩ B|` — inclusion–exclusion, which ties all of
//!   union, intersection, and area together in one equation
//! - `A ∪ A = A ∩ A = A` — idempotence
//! - `A \ B` and `B \ A` are disjoint, and `(A \ B) ∪ (A ∩ B) = A`
//! - `A △ B = (A \ B) ∪ (B \ A)`, and `|A △ B| = |A ∪ B| − |A ∩ B|`
//! - `|A ∩ B| ≤ min(|A|, |B|) ≤ max(|A|, |B|) ≤ |A ∪ B|` — monotone bounds
//! - `A ⊆ B` implies `A ∪ B = B` and `A ∩ B = A`
//!
//! An implementation can only satisfy all of these simultaneously by being right.
//! Getting one operation wrong on a pinched ring breaks at least one relation, and
//! because each relation mentions several operations, a compensating error in two
//! places does not hide: it has to be consistent across every equation at once.
//!
//! The monotone bound is the one that actually caught a defect in this repository
//! before (#44: a union smaller than one of its operands, from pinch geometry), so
//! it is asserted on every pair here rather than only where it seems likely.

use super::*;
use crate::space::r2::{BooleanOps, MultiPolygon};

/// Relative tolerance for an area identity. The areas here are O(1)–O(100), and a
/// sweep that reconstructs a boundary accumulates rounding per vertex, so an exact
/// comparison would fail for reasons that are not defects. Anything structurally
/// wrong is off by a whole piece of geometry, which is orders of magnitude larger.
const TOL: f64 = 1e-9;

fn close(a: f64, b: f64, what: &str) {
    let bound = TOL * a.abs().max(b.abs()).max(1.0);
    assert!(
        (a - b).abs() <= bound,
        "{what}: {a} vs {b} (difference {} exceeds {bound})",
        (a - b).abs()
    );
}

/// Assert every set-algebra relation on one pair. Called for each case below, so a
/// new fixture gets the whole battery rather than whichever check its author
/// thought of.
fn assert_set_algebra(name: &str, a: &MultiPolygon, b: &MultiPolygon) {
    let (area_a, area_b) = (ha(a), ha(b));
    let union = a.union(b);
    let inter = a.intersection(b);
    let diff_ab = a.difference(b);
    let diff_ba = b.difference(a);
    let sym = a.xor(b);

    let (u, i, dab, dba, s) = (ha(&union), ha(&inter), ha(&diff_ab), ha(&diff_ba), ha(&sym));

    for (label, v) in [
        ("A", area_a),
        ("B", area_b),
        ("union", u),
        ("intersection", i),
        ("A\\B", dab),
        ("B\\A", dba),
        ("symdiff", s),
    ] {
        assert!(
            v.is_finite() && v >= -TOL,
            "{name}: {label} area must be finite and non-negative, got {v}"
        );
    }

    // Inclusion–exclusion: the relation that ties everything together.
    close(
        u + i,
        area_a + area_b,
        &format!("{name}: |A u B| + |A n B| = |A| + |B|"),
    );
    // Difference decomposition.
    close(dab + i, area_a, &format!("{name}: |A\\B| + |A n B| = |A|"));
    close(dba + i, area_b, &format!("{name}: |B\\A| + |A n B| = |B|"));
    // Symmetric difference.
    close(s, u - i, &format!("{name}: |A sym B| = |A u B| - |A n B|"));
    close(
        s,
        dab + dba,
        &format!("{name}: |A sym B| = |A\\B| + |B\\A|"),
    );

    // Monotone bounds. This is the shape of the #44 defect: a union that came back
    // smaller than an operand.
    assert!(
        u >= area_a - TOL * area_a.max(1.0) && u >= area_b - TOL * area_b.max(1.0),
        "{name}: union {u} is smaller than an operand ({area_a}, {area_b})"
    );
    assert!(
        i <= area_a + TOL * area_a.max(1.0) && i <= area_b + TOL * area_b.max(1.0),
        "{name}: intersection {i} exceeds an operand ({area_a}, {area_b})"
    );
    assert!(
        i <= u + TOL * u.max(1.0),
        "{name}: intersection {i} exceeds union {u}"
    );

    // Commutativity: the sweep processes the operands in order, so this is not free.
    close(ha(&b.union(a)), u, &format!("{name}: union commutes"));
    close(
        ha(&b.intersection(a)),
        i,
        &format!("{name}: intersection commutes"),
    );
    close(
        ha(&b.xor(a)),
        s,
        &format!("{name}: symmetric difference commutes"),
    );

    // Idempotence on each operand.
    close(ha(&a.union(a)), area_a, &format!("{name}: A u A = A"));
    close(
        ha(&a.intersection(a)),
        area_a,
        &format!("{name}: A n A = A"),
    );
    assert!(
        ha(&a.difference(a)) <= TOL,
        "{name}: A \\ A must be empty, got {}",
        ha(&a.difference(a))
    );
}

/// A ring that touches itself at one vertex — a figure eight pinched at the origin.
/// Two lobes meeting at a single point is the configuration where "is this one
/// region or two" has no locally-decidable answer.
fn pinched_bowtie() -> MultiPolygon {
    hmp(&[&[&[
        (0.0, 0.0),
        (-2.0, 1.0),
        (-2.0, -1.0),
        (0.0, 0.0),
        (2.0, -1.0),
        (2.0, 1.0),
    ]]])
}

/// A square whose hole touches the shell along an edge, so the "interior" is an
/// annulus that is open on one side.
fn hole_touching_shell() -> MultiPolygon {
    hmp(&[&[
        &[(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0)],
        // Hole flush against the left edge.
        &[(0.0, 2.0), (3.0, 2.0), (3.0, 4.0), (0.0, 4.0)],
    ]])
}

/// Consecutive duplicate and collinear vertices, which a sweep must normalize away
/// without changing the region.
fn redundant_vertices() -> MultiPolygon {
    hmp(&[&[&[
        (0.0, 0.0),
        (0.0, 0.0),
        (2.0, 0.0),
        (4.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (2.0, 4.0),
        (0.0, 4.0),
        (0.0, 2.0),
    ]]])
}

/// A sliver: a triangle a few ulp tall. Degenerate in practice but not in theory.
fn sliver() -> MultiPolygon {
    hmp(&[&[&[(0.0, 0.0), (4.0, 0.0), (4.0, 1e-13)]]])
}

fn unit_square(x: f64, y: f64, size: f64) -> MultiPolygon {
    hmp(&[&[&[(x, y), (x + size, y), (x + size, y + size), (x, y + size)]]])
}

#[test]
fn set_algebra_holds_on_degenerate_pairs() {
    let square = unit_square(0.0, 0.0, 4.0);
    let offset = unit_square(2.0, 2.0, 4.0);
    let contained = unit_square(1.0, 1.0, 1.0);
    let disjoint = unit_square(100.0, 100.0, 1.0);
    let touching_edge = unit_square(4.0, 0.0, 4.0);
    let touching_corner = unit_square(4.0, 4.0, 4.0);

    let cases: Vec<(&str, MultiPolygon, MultiPolygon)> = vec![
        ("overlapping squares", square.clone(), offset),
        ("one contains the other", square.clone(), contained),
        ("disjoint", square.clone(), disjoint),
        ("edge-to-edge contact", square.clone(), touching_edge),
        ("corner-to-corner contact", square.clone(), touching_corner),
        ("identical", square.clone(), square.clone()),
        ("pinched bowtie vs square", pinched_bowtie(), square.clone()),
        (
            "pinched bowtie vs itself",
            pinched_bowtie(),
            pinched_bowtie(),
        ),
        ("hole touching shell", hole_touching_shell(), square.clone()),
        (
            "hole touching shell vs itself",
            hole_touching_shell(),
            hole_touching_shell(),
        ),
        ("redundant vertices", redundant_vertices(), square.clone()),
        (
            "redundant vs clean same region",
            redundant_vertices(),
            unit_square(0.0, 0.0, 4.0),
        ),
        ("sliver vs square", sliver(), square.clone()),
        ("sliver vs sliver", sliver(), sliver()),
        (
            "two slivers crossing",
            sliver(),
            hmp(&[&[&[(0.0, -1.0), (1e-13, -1.0), (1e-13, 3.0)]]]),
        ),
    ];

    for (name, a, b) in cases {
        assert_set_algebra(name, &a, &b);
    }
}

#[test]
fn redundant_vertices_describe_the_same_region_as_the_clean_ring() {
    // A sweep that mishandles duplicate or collinear vertices usually still produces
    // *a* polygon — just the wrong one. Comparing areas catches that, and the
    // symmetric difference being empty is the stronger claim: not merely equal in
    // measure, but the same set.
    let messy = redundant_vertices();
    let clean = unit_square(0.0, 0.0, 4.0);
    close(ha(&messy), ha(&clean), "redundant ring has the clean area");
    let sym = messy.xor(&clean);
    assert!(
        ha(&sym) <= TOL,
        "the two descriptions should denote the same region, symmetric difference {}",
        ha(&sym)
    );
}

#[test]
fn containment_collapses_union_and_intersection() {
    // `A subset B` is the case with an exactly known answer, so it is worth pinning
    // directly rather than only through the inequalities.
    let outer = unit_square(0.0, 0.0, 10.0);
    for inner in [
        unit_square(1.0, 1.0, 2.0),
        unit_square(0.0, 0.0, 1.0),  // sharing a corner
        unit_square(0.0, 4.0, 1.0),  // flush against an edge
        unit_square(0.0, 0.0, 10.0), // equal, the boundary of containment
    ] {
        let area_inner = ha(&inner);
        let area_outer = ha(&outer);
        close(
            ha(&outer.union(&inner)),
            area_outer,
            "A u B = B when A is inside B",
        );
        close(
            ha(&outer.intersection(&inner)),
            area_inner,
            "A n B = A when A is inside B",
        );
        close(
            ha(&inner.difference(&outer)),
            0.0,
            "A \\ B is empty when A is inside B",
        );
        close(
            ha(&outer.difference(&inner)),
            area_outer - area_inner,
            "B \\ A is the annulus",
        );
    }
}

#[test]
fn disjoint_operands_do_not_interact() {
    let a = unit_square(0.0, 0.0, 1.0);
    for b in [
        unit_square(5.0, 5.0, 1.0),
        unit_square(1.0 + 1e-12, 0.0, 1.0), // a hair apart
        unit_square(-10.0, -10.0, 3.0),
    ] {
        close(
            ha(&a.union(&b)),
            ha(&a) + ha(&b),
            "disjoint union is the sum of areas",
        );
        assert!(
            ha(&a.intersection(&b)) <= TOL,
            "disjoint operands have empty intersection, got {}",
            ha(&a.intersection(&b))
        );
        close(ha(&a.difference(&b)), ha(&a), "nothing to subtract");
    }
}

#[test]
fn empty_operands_behave_as_the_identity_and_annihilator() {
    let empty = MultiPolygon::new(vec![]);
    let a = unit_square(0.0, 0.0, 3.0);
    close(ha(&a.union(&empty)), ha(&a), "A u empty = A");
    close(ha(&empty.union(&a)), ha(&a), "empty u A = A");
    assert!(ha(&a.intersection(&empty)) <= TOL, "A n empty = empty");
    assert!(ha(&empty.intersection(&a)) <= TOL, "empty n A = empty");
    close(ha(&a.difference(&empty)), ha(&a), "A \\ empty = A");
    assert!(ha(&empty.difference(&a)) <= TOL, "empty \\ A = empty");
    close(ha(&a.xor(&empty)), ha(&a), "A sym empty = A");
    // And empty against empty, which must not divide by anything.
    for op in [
        ha(&empty.union(&empty)),
        ha(&empty.intersection(&empty)),
        ha(&empty.difference(&empty)),
        ha(&empty.xor(&empty)),
    ] {
        assert!(
            op.abs() <= TOL,
            "every op on two empties is empty, got {op}"
        );
    }
}

#[test]
fn a_chain_of_unions_is_associative_in_measure() {
    // Associativity is a different stress than commutativity: it changes which
    // boundaries have already been reconstructed when the next operand arrives, so
    // it exercises the sweep on its own output.
    let a = unit_square(0.0, 0.0, 4.0);
    let b = unit_square(2.0, 2.0, 4.0);
    let c = unit_square(3.0, 0.0, 4.0);
    let left = a.union(&b).union(&c);
    let right = a.union(&b.union(&c));
    close(ha(&left), ha(&right), "(A u B) u C = A u (B u C)");
    let sym = left.xor(&right);
    assert!(
        ha(&sym) <= TOL,
        "the two groupings should be the same region, symmetric difference {}",
        ha(&sym)
    );
    // De Morgan on a bounded universe: (A u B) \ C = (A \ C) u (B \ C).
    let lhs = a.union(&b).difference(&c);
    let rhs = a.difference(&c).union(&b.difference(&c));
    close(ha(&lhs), ha(&rhs), "(A u B) \\ C = (A \\ C) u (B \\ C)");
}
