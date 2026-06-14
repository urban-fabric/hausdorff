//! Geometry validity per the OGC Simple Feature Access rules, mirroring `geo`'s
//! `Validation` trait: `is_valid()` for the boolean verdict, and
//! `check_validation()` / `validation_errors()` for typed `InvalidXxx` reasons.
//!
//! Checks per type: finite coordinates everywhere; a [`LineString`] needs at
//! least 2 distinct points; polygon rings need at least 4 points (closed, so 3
//! distinct), must be simple (no self-intersection — spikes and crossings are
//! caught), interior rings must lie inside the exterior ring and may touch it
//! or each other only at points (never along a line, never overlapping or
//! nested); [`MultiPolygon`] members must additionally have disjoint interiors.
//! Ring containment/overlap is decided with the native DE-9IM
//! [`relate`](super::relate) machinery, ring simplicity with the crate's
//! segment-intersection predicate. Verdicts are differential-tested against
//! `geo`'s `Validation` in `tests/validation.rs`.

use super::intersects::segments_intersect;
use super::relate::relate;
use super::{
    Coord, Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    Rect, Triangle,
};

mod collections;
mod common;
mod polygonal;
mod primitives;

pub use collections::*;
pub use common::*;
pub use polygonal::*;
pub use primitives::*;

#[cfg(test)]
mod tests;
