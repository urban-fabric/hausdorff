//! The same geometry, moved or resampled — vertices in, vertices out, topology
//! unchanged.
//!
//! That invariant is what separates this group from [`construct`](super::construct)
//! and it is a promise, not a description. [`affine`] and [`isometry`] move every
//! vertex by a map that cannot change which side of an edge a point is on, so
//! validity is preserved by construction. [`densify`] only adds vertices along
//! existing edges, and [`remove_repeated_points`] only drops ones that coincide.
//!
//! [`simplify`] is the exception that proves the rule: Douglas–Peucker *can*
//! introduce a self-intersection, which is why [`simplify_topology`] exists
//! separately and vetoes any removal that would cross a retained edge. The two are
//! not a fast path and a slow path for the same operation — they answer different
//! questions, and a caller who needs a valid result afterwards needs the second.

pub(crate) mod affine;
pub(crate) mod densify;
pub(crate) mod isometry;
pub(crate) mod map_coords;
pub(crate) mod remove_repeated_points;
pub(crate) mod simplify;
pub(crate) mod simplify_topology;
