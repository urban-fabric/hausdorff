//! A new geometry built from existing ones.
//!
//! The largest group, and the one where the hard problems are, because producing
//! geometry means deciding topology: which rings are shells, which are holes, what
//! happens where two boundaries touch. Four of these — [`bool_ops`],
//! [`make_valid`], [`polygonize`] and [`buffer`] — are *arrangement* algorithms
//! underneath: node every segment against every other on a shared vertex pool,
//! classify the pieces, walk them into rings, nest the rings. They share
//! [`make_valid`]'s `node` and `trace` rather than each carrying a copy, because
//! an independent noding implementation is an independent set of tolerance
//! decisions, and two sets drift.
//!
//! [`triangulate`] and the hulls ([`convex_hull`], [`concave_hull`]) are the
//! other kind: incremental or sweep constructions with no boolean classification,
//! where the difficulty is degeneracy — cocircular points, collinear runs — rather
//! than topology.
//!
//! [`closest_point`], [`interior_point`] and [`line_locate_interpolate`] return
//! points rather than areas, which makes them look like measures. They are here
//! because the point is *constructed* to satisfy a property (on the boundary,
//! strictly inside, a fraction along), not read off the input.

pub(crate) mod bool_ops;
pub(crate) mod buffer;
pub(crate) mod closest_point;
pub(crate) mod concave_hull;
pub(crate) mod convex_hull;
pub(crate) mod interior_point;
pub(crate) mod line_locate_interpolate;
pub(crate) mod make_valid;
pub(crate) mod polygonize;
pub(crate) mod triangulate;
