//! The geometry types themselves, and the traits that walk their coordinates.
//!
//! This is the group everything else depends on and which depends on nothing
//! here: [`Coord`](crate::space::r2::Coord) and the shapes built from it, the
//! [`Geometry`](crate::space::r2::Geometry) enum that unifies them, and the two
//! traits — [`CoordsIter`](crate::space::r2::CoordsIter) and
//! [`HasDimensions`](crate::space::r2::HasDimensions) — that let an algorithm
//! accept any shape without matching on it.
//!
//! [`geometry_dispatch`] is here rather than with the algorithms it forwards to,
//! which is worth saying because it looks misfiled. It contains no geometry: it is
//! the `impl Trait for Geometry` arms that route a heterogeneous value to the
//! concrete implementation, so it belongs to the enum rather than to any one
//! algorithm. Putting it beside the enum keeps "what can a `Geometry` do" in one
//! place; distributing it would scatter twelve one-line `match`es across five
//! groups.
//!
//! [`convex_polygon`] carries its own type because convexity is a *proof* that
//! downstream code consumes — [`ConvexPolygon`](crate::space::r2::ConvexPolygon)
//! unlocks the O(log n) containment and distance paths — so it is a model concern
//! and not a construction one, even though a hull is what usually produces it.

pub(crate) mod convex_polygon;
pub(crate) mod coord;
pub(crate) mod coords_iter;
pub(crate) mod dimensions;
pub(crate) mod geometry;
pub(crate) mod geometry_dispatch;
pub(crate) mod line_string;
pub(crate) mod polygon;
pub(crate) mod primitives;
