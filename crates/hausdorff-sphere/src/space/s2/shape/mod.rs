//! The spherical regions: what a geometry *is* on the sphere.
//!
//! The largest group, because a sphere library is mostly shapes. They form a
//! rough hierarchy of generality — [`cap`] and [`rect`] are the two bounding
//! primitives, [`loops`] is the general spherical ring, [`convex_loop`] is the
//! same thing carrying a convexity proof that unlocks O(log n) containment, and
//! [`polygon`] / [`polyline`] / [`triangle`] / [`simplex`] are built from those.
//!
//! [`region`] is the group's interface to [`cell`](super::cell): it defines the
//! `Region` trait that `cap`, `loops` and `polygon` implement, and the
//! `RegionCoverer` that turns any of them into cells. It sits here rather than in
//! `cell` because the trait is a statement about shapes; `cover` consumes it from
//! the other side.

pub(crate) mod cap;
pub(crate) mod convex_loop;
pub(crate) mod loops;
pub(crate) mod polygon;
pub(crate) mod polyline;
pub(crate) mod rect;
pub(crate) mod region;
pub(crate) mod simplex;
pub(crate) mod triangle;
