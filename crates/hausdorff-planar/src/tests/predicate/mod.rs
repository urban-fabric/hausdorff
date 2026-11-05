//! Geometry in, verdict out.
//!
//! The group with the most self-consistency checks available to it, and they are
//! the strongest oracle here: [`binary`] does not need to know the right answer
//! for a random pair, only that `eval` agrees with the matrix it dispatches to and
//! that `transpose` holds. Relations that must agree pin each other down without a
//! reference implementation.
//!
//! [`exact_predicates`] is the exception that anchors the rest: it checks the
//! integer kernels against `f64` on inputs where `f64` is known to be wrong, so it
//! is the one suite here with an unarguable answer.
//!
//! [`validation`] and [`prepared`] are included because both produce verdicts —
//! "is this geometry legal" and "does this cached form still answer the same way
//! as the uncached one".

mod binary;
mod exact_predicates;
mod prepared;
mod relate;
mod validation;
mod winding;
