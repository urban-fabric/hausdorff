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
//! [`relate`](super::relate()) machinery, ring simplicity with the crate's
//! segment-intersection predicate. Verdicts are differential-tested against
//! `geo`'s `Validation` in `tests/validation.rs`.
//!
//! ## Why this is a namespace and not part of the facade
//!
//! [`space::r2`](super) re-exports every other module's items so each has exactly
//! one path; this module is the single exception, and it earns it twice over. Its
//! thirteen items are one vocabulary — [`RingRole`], [`CoordIndex`] and
//! [`GeometryIndex`] locate a fault, the eleven `Invalid*` enums describe one — and
//! they read better qualified: `validation::InvalidPolygon` says what a bare
//! `InvalidPolygon` does not.
//!
//! The deciding evidence is [`GeometryIndex`], which here means "which member of
//! the collection is at fault". There is an unrelated public `GeometryIndex` in the
//! facade that is a spatial index over geometries. Two different public types, one
//! name — which module paths were quietly hiding, and which a flat namespace cannot
//! hold at all. That is what a separate vocabulary looks like.

use super::{
    Coord, Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    Rect, Triangle,
};
use crate::space::r2::relate;
use crate::space::r2::segments_intersect;

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
