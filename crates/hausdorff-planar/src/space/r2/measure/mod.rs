//! Scalars derived from a geometry: areas, lengths, distances, centres, extents.
//!
//! The common shape is a geometry in and a number out, with no geometry returned —
//! which is what separates this group from [`construct`](super::construct). A
//! bounding *rectangle* is a geometry, but [`bounding_rect`] is here anyway
//! because the rectangle is a summary of extent rather than a constructed shape;
//! the same goes for the circle and oriented rectangle in [`bounding_shapes`].
//! The test is whether the output is *about* the input or *made from* it.
//!
//! Two sub-groups are worth distinguishing, because they fail differently.
//! [`area`], [`centroid`] and [`geodesic_area`] are **sums over a ring**, so their
//! enemy is cancellation — see `util::signed_area2` for the conditioning that
//! governs all three. [`min_distance`], [`set_distance`] and [`line_measures`] are
//! **minima over pairs**, so their enemy is cost, and they carry the pruning
//! (`seg_box`, workspaces, early exits) to match.
//!
//! [`line_measures`] holds the metric traits rather than a single algorithm:
//! `Euclidean` is planar, while `Haversine`, `Geodesic` and `Rhumb` delegate to
//! `hausdorff-base`'s validated kernels. They live together because they are
//! interchangeable at the call site, which is the whole point of the trait.

pub(crate) mod area;
pub(crate) mod bounding_rect;
pub(crate) mod bounding_shapes;
pub(crate) mod centroid;
pub(crate) mod extremes;
pub(crate) mod geodesic_area;
pub(crate) mod line_measures;
pub(crate) mod min_distance;
pub(crate) mod set_distance;
