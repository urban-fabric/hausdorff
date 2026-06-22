//! `r3` — three-dimensional Euclidean geometry: the [`Vector`] type, mirroring
//! golang/geo's `r3` package.
//!
//! `r3::Vector` is the substrate of the spherical world: an `S2Point` (in `hausdorff-sphere`)
//! is a *unit* `r3::Vector`, and spherical measures/predicates are computed from
//! vector dot/cross products rather than lat/lon, which is what makes them free
//! of pole and antimeridian singularities.

use euclid::numeric::{Boolean as BooleanMask, Compare, Real};

use crate::space::s1::Angle;

mod vector;
pub use vector::{Axis, Vector};

mod aabb;
pub use aabb::Aabb;

mod plane;
pub use plane::Plane;

#[cfg(test)]
mod tests;
