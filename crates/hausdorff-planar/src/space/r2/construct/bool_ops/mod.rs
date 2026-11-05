//! `r2::bool_ops` — polygon boolean operations (intersection, union,
//! difference, symmetric difference) via the Martinez–Rueda–Feito sweep-line
//! algorithm. Pure, no external deps.
//!
//! Ported from the reference implementation by F. Martínez et al. (2009/2013);
//! the [`BooleanOps`] trait mirrors `geo`'s so call sites migrate by import.
//! Input winding is normalised (exterior CCW, holes CW) so arbitrarily-wound
//! polygons work. Validated against `geo`'s `BooleanOps` over fixed
//! configurations (disjoint, nested, crossing, concave, touching/collinear,
//! identical, holes, multipolygons) plus ~900 randomized rectangle/triangle
//! pairs; result areas and the set-algebra invariants match.
//!
//! Sweep events live in a flat arena (`Vec<Event>`) and reference each other by
//! `u32` index rather than `Rc<RefCell<…>>`: this removes a heap allocation and
//! reference-count per event, and the runtime borrow-checks the comparators paid
//! on every comparison. The event queue is a small index min-heap whose ordering
//! reads the live arena (the comparison keys can change as segments are divided,
//! exactly as the reference algorithm assumes).
//!
//! **T-vertices.** A *T-vertex* — a vertex of one polygon lying exactly on the
//! interior of the other's edge — is handled by the overlay: the touching,
//! integer-grid-tiling, overlapping, and many-vertex staircase configurations are
//! regression-tested against `geo` across all four ops (the `t_vertex_*` cases in
//! `tests/bool_ops.rs`). Beyond those, exact-coincidence robustness is the usual
//! caveat of a floating-point sweep; for purely topological queries,
//! [`crate::space::r2::relate()`] stays clipper-free.

use crate::space::r2::{Coord, MultiPolygon, Polygon};

mod assemble;
mod connect;

mod event;
mod intersection;
mod sweep;
use sweep::run;

mod unary_union;
pub use unary_union::{UnaryUnionOperand, unary_union};

/// Which boolean operation to compute.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Operation {
    Intersection,
    Union,
    Difference,
    Xor,
}

/// Polygon set operations on [`Polygon`]/[`MultiPolygon`].
pub trait BooleanOps {
    /// Region in both operands.
    fn intersection(&self, other: &Self) -> MultiPolygon;
    /// Region in either operand.
    fn union(&self, other: &Self) -> MultiPolygon;
    /// Region in `self` but not `other`.
    fn difference(&self, other: &Self) -> MultiPolygon;
    /// Region in exactly one operand.
    fn xor(&self, other: &Self) -> MultiPolygon;
}

impl BooleanOps for MultiPolygon {
    fn intersection(&self, other: &Self) -> MultiPolygon {
        run(self, other, Operation::Intersection)
    }
    fn union(&self, other: &Self) -> MultiPolygon {
        run(self, other, Operation::Union)
    }
    fn difference(&self, other: &Self) -> MultiPolygon {
        run(self, other, Operation::Difference)
    }
    fn xor(&self, other: &Self) -> MultiPolygon {
        run(self, other, Operation::Xor)
    }
}

impl BooleanOps for Polygon {
    fn intersection(&self, other: &Self) -> MultiPolygon {
        run(&one(self), &one(other), Operation::Intersection)
    }
    fn union(&self, other: &Self) -> MultiPolygon {
        run(&one(self), &one(other), Operation::Union)
    }
    fn difference(&self, other: &Self) -> MultiPolygon {
        run(&one(self), &one(other), Operation::Difference)
    }
    fn xor(&self, other: &Self) -> MultiPolygon {
        run(&one(self), &one(other), Operation::Xor)
    }
}

fn one(p: &Polygon) -> MultiPolygon {
    MultiPolygon::new(vec![p.clone()])
}

fn ring_vertices(ring: &[Coord]) -> &[Coord] {
    if ring.len() >= 2 && ring[0] == ring[ring.len() - 1] {
        &ring[..ring.len() - 1]
    } else {
        ring
    }
}
