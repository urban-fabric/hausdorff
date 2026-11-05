//! Questions about geometries answered with a verdict rather than a geometry.
//!
//! [`contains`] and [`intersects`] are the boolean pair; [`relate`] is the
//! complete DE-9IM matrix the two of them are special cases of; [`winding`]
//! classifies a ring's orientation; [`prepared`] is not a predicate but the
//! caching layer that makes repeated ones cheap, which is why it sits with them.
//!
//! [`exact`] is the group's reason for existing as a group. Every other predicate
//! here is `f64` and therefore approximate near a degeneracy; `exact` answers the
//! same orientation and in-circle questions in exact arithmetic and *reports
//! overflow* rather than absorbing it. Keeping the two kinds adjacent is
//! deliberate — the choice between them is a real decision about cost and
//! certainty, and it is easier to make when both are in view.
//!
//! The tolerances in [`relate`] are the subtle part and are documented where they
//! are used: a predicate tolerance has to be relative to the geometry and
//! dimensionally consistent with the quantity it bounds, and `relate` has been
//! wrong in both directions at once — see `relate::predicates::point_on_segment`.

pub(crate) mod contains;
pub(crate) mod exact;
pub(crate) mod intersects;
pub(crate) mod prepared;
pub(crate) mod relate;
pub(crate) mod winding;
