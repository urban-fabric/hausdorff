//! In-house S2 cell-indexing implementation — a from-scratch port of the S2
//! Geometry cell hierarchy (Google's S2), giving the workspace a second
//! hierarchical index system alongside [`crate::h3`].
//!
//! Covers the indexing core: the 64-bit [`S2CellId`] codec + validity, the
//! cube-to-sphere quadratic projection (lat/lon ↔ cell), `parent` / `children`
//! / child positions, `range_min`/`range_max` containment, edge/vertex/all
//! neighbors, cell centers, Hilbert-order traversal, and hex tokens. Cells
//! share [`crate::geo::GeoPoint`] at the lat/lon boundary, exactly like the H3
//! side, so a `point → cell` mapping uses one coordinate type across both
//! systems.
//!
//! Ported from the Go reference implementation (`github.com/golang/geo/s2`),
//! vendored under `project-management/libraries/s2`. Like that implementation
//! we use only the quadratic ST↔UV projection. Full S2 geometry (regions,
//! polygons, coverings) is intentionally out of scope.

mod cellid;
mod coords;
mod point;
mod vector;

pub use cellid::S2CellId;
pub use point::S2Point;

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests;
