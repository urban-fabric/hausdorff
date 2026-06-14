use std::fmt;

use super::super::remove_repeated_points::RemoveRepeatedPoints;
use super::*;

/// Check whether a geometry is valid and report the reason(s) of invalidity.
pub trait Validation {
    type Error: std::error::Error;

    /// True iff the geometry is valid.
    fn is_valid(&self) -> bool {
        self.check_validation().is_ok()
    }

    /// Every reason of invalidity. Later errors may be obscured by earlier ones
    /// (e.g. overlap checks assume valid member polygons), so re-validate after
    /// repairing.
    fn validation_errors(&self) -> Vec<Self::Error> {
        let mut validation_errors = Vec::new();
        self.visit_validation(Box::new(|problem| {
            validation_errors.push(problem);
            Ok::<(), Self::Error>(())
        }))
        .expect("the error collector never fails");
        validation_errors
    }

    /// The first reason of invalidity, as a `Result`.
    fn check_validation(&self) -> Result<(), Self::Error> {
        self.visit_validation(Box::new(Err))
    }

    /// Run the validity checks, calling `handle_validation_error` for each
    /// failure; the visit stops early when the handler returns `Err`.
    fn visit_validation<T>(
        &self,
        handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T>;
}

/// The role of a ring in a [`Polygon`].
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum RingRole {
    Exterior,
    Interior(usize),
}

impl fmt::Display for RingRole {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            RingRole::Exterior => write!(f, "exterior ring"),
            RingRole::Interior(idx) => write!(f, "interior ring at index {idx}"),
        }
    }
}

/// The position of the problem member in a multi-geometry, starting at 0.
#[derive(Debug, PartialEq, Clone)]
pub struct GeometryIndex(pub usize);

/// The index of the offending coordinate in its geometry.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CoordIndex(pub usize);

pub(super) fn coord_is_finite(c: &Coord) -> bool {
    c.is_finite()
}

/// Fewer distinct (consecutively deduplicated) points than a non-degenerate
/// line (2) or closed ring (4) requires.
pub(super) fn has_too_few_points(ls: &LineString, is_ring: bool) -> bool {
    let min_points = if is_ring { 4 } else { 2 };
    ls.remove_repeated_points().0.len() < min_points
}

/// Does the line string cross or touch itself anywhere except at consecutive
/// segments' shared vertex (and a closed ring's seam)? Mirrors `geo`'s check.
pub(super) fn has_self_intersection(ls: &LineString) -> bool {
    for (i, a) in ls.lines().enumerate() {
        for (j, b) in ls.lines().enumerate() {
            if i != j
                && segments_intersect(a.start, a.end, b.start, b.end)
                && a.start != b.end
                && a.end != b.start
            {
                return true;
            }
        }
    }
    false
}

/// View one ring as a standalone hole-free polygon, for [`relate`]-based ring
/// containment/overlap checks.
pub(super) fn ring_as_polygon(ring: &LineString) -> Geometry {
    Geometry::Polygon(Polygon::new(ring.clone(), vec![]))
}
