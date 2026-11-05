//! The unit-sphere substrate: the point type, the robust predicates, and the
//! quantities computed directly from points.
//!
//! Everything else in `s2` rests on this group, and it rests on nothing here —
//! which is why it is worth isolating. [`point`] holds [`S2Point`,
//! the unit vector that makes the whole module projection-free; [`predicates`]
//! holds the orientation sign that every downstream decision reduces to; and
//! [`measures`] and [`edges`] are the two layers built straight on top of those
//! (areas and angles, crossing and intersection).
//!
//! The ordering inside the group is strict and worth preserving: `predicates`
//! depends on nothing, `measures` and `edges` depend only on `predicates`, and
//! `great_circle` only on the point type. A new module here that needs a *shape*
//! is in the wrong group.

pub(crate) mod edges;
pub(crate) mod great_circle;
pub(crate) mod measures;
pub(crate) mod point;
pub(crate) mod predicates;
