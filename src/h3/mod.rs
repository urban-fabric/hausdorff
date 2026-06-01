//! In-house H3 implementation — a from-scratch port of the H3 hexagonal
//! geospatial index system, replacing the external `h3o` crate.
//!
//! Covers the operations the workspace uses: the u64 index codec + validity,
//! the icosahedral FaceIJK gnomonic projection (`for_point` / `center`),
//! `parent` / `children`, child positions, compaction/uncompaction, `grid_disk`
//! (k-ring), `grid_distance`, local IJ, grid paths, directed edges, canonical
//! vertices, exact cell area, centroid polygon cover, cell-set boundary
//! dissolve, `is_pentagon`, and cell/edge boundaries.
//!
//! Every operation is differential-tested against `h3o` (kept as a
//! dev-dependency oracle only) plus a random-index validity sweep — see
//! `tests.rs`. Constant tables (base cells, faces, FaceIJK base-cell lookup) are
//! the standard H3 data, lifted verbatim; the algorithms are ported from the H3
//! C reference / h3o.

mod base_cells;
mod children;
mod coord;
mod face;
mod faceijk;
mod grid;
mod index;
mod polygon;

pub use index::{
    H3Cell, H3CompactionError, H3DirectedEdge, H3LocalIJ, H3Vertex, InvalidH3Cell,
    InvalidH3DirectedEdge, InvalidH3Vertex, MAX_RESOLUTION,
};
pub use polygon::{H3BoundaryRing, H3Polygon, H3PolygonError};

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests;
