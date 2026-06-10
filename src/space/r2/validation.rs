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

use std::fmt;

use super::dimensions::HasDimensions;
use super::intersects::segments_intersect;
use super::relate::relate;
use super::remove_repeated_points::RemoveRepeatedPoints;
use super::util::orient;
use super::{
    Coord, Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    Rect, Triangle,
};

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

fn coord_is_finite(c: &Coord) -> bool {
    c.is_finite()
}

/// Fewer distinct (consecutively deduplicated) points than a non-degenerate
/// line (2) or closed ring (4) requires.
fn has_too_few_points(ls: &LineString, is_ring: bool) -> bool {
    let min_points = if is_ring { 4 } else { 2 };
    ls.remove_repeated_points().0.len() < min_points
}

/// Does the line string cross or touch itself anywhere except at consecutive
/// segments' shared vertex (and a closed ring's seam)? Mirrors `geo`'s check.
fn has_self_intersection(ls: &LineString) -> bool {
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

/// Why a [`Point`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidPoint {
    /// A valid [`Point`] must have finite coordinates.
    NonFiniteCoord,
}

impl fmt::Display for InvalidPoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidPoint::NonFiniteCoord => write!(f, "point has non-finite coordinates"),
        }
    }
}

impl std::error::Error for InvalidPoint {}

impl Validation for Point {
    type Error = InvalidPoint;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        if !coord_is_finite(&self.0) {
            handle_validation_error(InvalidPoint::NonFiniteCoord)?;
        }
        Ok(())
    }
}

/// Why a [`Line`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidLine {
    /// A valid [`Line`] must have two distinct points.
    IdenticalCoords,
    /// A valid [`Line`] must have finite coordinates.
    NonFiniteCoord(CoordIndex),
}

impl fmt::Display for InvalidLine {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidLine::IdenticalCoords => write!(f, "line has identical coordinates"),
            InvalidLine::NonFiniteCoord(idx) => write!(
                f,
                "{} coordinate is non-finite",
                if idx.0 == 0 { "start" } else { "end" }
            ),
        }
    }
}

impl std::error::Error for InvalidLine {}

impl Validation for Line {
    type Error = InvalidLine;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        if !coord_is_finite(&self.start) {
            handle_validation_error(InvalidLine::NonFiniteCoord(CoordIndex(0)))?;
        }
        if !coord_is_finite(&self.end) {
            handle_validation_error(InvalidLine::NonFiniteCoord(CoordIndex(1)))?;
        }
        if self.start == self.end {
            handle_validation_error(InvalidLine::IdenticalCoords)?;
        }
        Ok(())
    }
}

/// Why a [`LineString`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidLineString {
    /// A valid [`LineString`] must have at least 2 distinct points.
    TooFewPoints,
    /// A valid [`LineString`] must have finite coordinates.
    NonFiniteCoord(CoordIndex),
}

impl fmt::Display for InvalidLineString {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidLineString::TooFewPoints => {
                write!(f, "line string must have at least 2 distinct points")
            }
            InvalidLineString::NonFiniteCoord(idx) => {
                write!(f, "coordinate at index {} is non-finite", idx.0)
            }
        }
    }
}

impl std::error::Error for InvalidLineString {}

impl Validation for LineString {
    type Error = InvalidLineString;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        if self.is_empty() {
            return Ok(());
        }
        if has_too_few_points(self, false) {
            handle_validation_error(InvalidLineString::TooFewPoints)?;
        }
        for (coord_idx, coord) in self.0.iter().enumerate() {
            if !coord_is_finite(coord) {
                handle_validation_error(InvalidLineString::NonFiniteCoord(CoordIndex(coord_idx)))?;
            }
        }
        Ok(())
    }
}

/// Why a [`Polygon`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidPolygon {
    /// A ring needs at least 4 points (3 distinct plus the closing point).
    TooFewPointsInRing(RingRole),
    /// A ring crosses or touches itself.
    SelfIntersection(RingRole),
    /// A coordinate is non-finite.
    NonFiniteCoord(RingRole, CoordIndex),
    /// An interior ring must be completely within the exterior ring.
    InteriorRingNotContainedInExteriorRing(RingRole),
    /// Two rings intersect over a 1-dimensional (line) section.
    IntersectingRingsOnALine(RingRole, RingRole),
    /// Two rings intersect over a 2-dimensional (area) section.
    IntersectingRingsOnAnArea(RingRole, RingRole),
}

impl fmt::Display for InvalidPolygon {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidPolygon::TooFewPointsInRing(ring) => {
                write!(f, "{ring} must have at least 3 distinct points")
            }
            InvalidPolygon::SelfIntersection(ring) => write!(f, "{ring} has a self-intersection"),
            InvalidPolygon::NonFiniteCoord(ring, idx) => {
                write!(f, "{ring} has a non-finite coordinate at index {}", idx.0)
            }
            InvalidPolygon::InteriorRingNotContainedInExteriorRing(ring) => {
                write!(f, "{ring} is not contained within the polygon's exterior")
            }
            InvalidPolygon::IntersectingRingsOnALine(ring_1, ring_2) => {
                write!(f, "{ring_1} and {ring_2} intersect on a line")
            }
            InvalidPolygon::IntersectingRingsOnAnArea(ring_1, ring_2) => {
                write!(f, "{ring_1} and {ring_2} intersect on an area")
            }
        }
    }
}

impl std::error::Error for InvalidPolygon {}

impl Validation for Polygon {
    type Error = InvalidPolygon;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        if self.is_empty() {
            return Ok(());
        }

        for (ring_idx, ring) in std::iter::once(self.exterior())
            .chain(self.interiors().iter())
            .enumerate()
        {
            if ring.is_empty() {
                continue;
            }
            let ring_role = if ring_idx == 0 {
                RingRole::Exterior
            } else {
                RingRole::Interior(ring_idx - 1)
            };

            if has_too_few_points(ring, true) {
                handle_validation_error(InvalidPolygon::TooFewPointsInRing(ring_role))?;
            }
            if has_self_intersection(ring) {
                handle_validation_error(InvalidPolygon::SelfIntersection(ring_role))?;
            }
            for (coord_idx, coord) in ring.0.iter().enumerate() {
                if !coord_is_finite(coord) {
                    handle_validation_error(InvalidPolygon::NonFiniteCoord(
                        ring_role,
                        CoordIndex(coord_idx),
                    ))?;
                }
            }
        }

        let has_interiors = self.interiors().iter().any(|ring| !ring.is_empty());
        if !has_interiors {
            return Ok(());
        }

        let exterior_poly = ring_as_polygon(self.exterior());
        for (interior_1_idx, interior_1) in self.interiors().iter().enumerate() {
            if interior_1.is_empty() {
                continue;
            }
            let ring_role_1 = RingRole::Interior(interior_1_idx);
            let interior_1_poly = ring_as_polygon(interior_1);
            let exterior_vs_interior = relate(&exterior_poly, &interior_1_poly);

            if !exterior_vs_interior.is_contains() {
                handle_validation_error(InvalidPolygon::InteriorRingNotContainedInExteriorRing(
                    ring_role_1,
                ))?;
            }
            if exterior_vs_interior.boundary_boundary_dim() == 1 {
                handle_validation_error(InvalidPolygon::IntersectingRingsOnALine(
                    RingRole::Exterior,
                    ring_role_1,
                ))?;
            }

            for (interior_2_idx, interior_2) in
                self.interiors().iter().enumerate().skip(interior_1_idx + 1)
            {
                if interior_2.is_empty() {
                    continue;
                }
                let ring_role_2 = RingRole::Interior(interior_2_idx);
                let im = relate(&interior_1_poly, &ring_as_polygon(interior_2));

                if im.interior_interior_dim() == 2 {
                    handle_validation_error(InvalidPolygon::IntersectingRingsOnAnArea(
                        ring_role_1,
                        ring_role_2,
                    ))?;
                }
                if im.boundary_boundary_dim() == 1 {
                    handle_validation_error(InvalidPolygon::IntersectingRingsOnALine(
                        ring_role_1,
                        ring_role_2,
                    ))?;
                }
            }
        }
        Ok(())
    }
}

/// View one ring as a standalone hole-free polygon, for [`relate`]-based ring
/// containment/overlap checks.
fn ring_as_polygon(ring: &LineString) -> Geometry {
    Geometry::Polygon(Polygon::new(ring.clone(), vec![]))
}

/// Why a [`MultiPoint`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidMultiPoint {
    /// Which member is invalid, and what was invalid about it.
    InvalidPoint(GeometryIndex, InvalidPoint),
}

impl fmt::Display for InvalidMultiPoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidMultiPoint::InvalidPoint(idx, err) => {
                write!(f, "point at index {} is invalid: {}", idx.0, err)
            }
        }
    }
}

impl std::error::Error for InvalidMultiPoint {}

impl Validation for MultiPoint {
    type Error = InvalidMultiPoint;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        for (i, point) in self.iter().enumerate() {
            point.visit_validation(Box::new(&mut |err| {
                handle_validation_error(InvalidMultiPoint::InvalidPoint(GeometryIndex(i), err))
            }))?;
        }
        Ok(())
    }
}

/// Why a [`MultiLineString`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidMultiLineString {
    /// Which member is invalid, and what was invalid about it.
    InvalidLineString(GeometryIndex, InvalidLineString),
}

impl fmt::Display for InvalidMultiLineString {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidMultiLineString::InvalidLineString(idx, err) => {
                write!(f, "line string at index {} is invalid: {}", idx.0, err)
            }
        }
    }
}

impl std::error::Error for InvalidMultiLineString {}

impl Validation for MultiLineString {
    type Error = InvalidMultiLineString;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        for (i, line_string) in self.iter().enumerate() {
            line_string.visit_validation(Box::new(&mut |err| {
                handle_validation_error(InvalidMultiLineString::InvalidLineString(
                    GeometryIndex(i),
                    err,
                ))
            }))?;
        }
        Ok(())
    }
}

/// Why a [`MultiPolygon`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidMultiPolygon {
    /// Which member is invalid, and what was invalid about it.
    InvalidPolygon(GeometryIndex, InvalidPolygon),
    /// Two members' interiors overlap (2-dimensional intersection).
    ElementsOverlaps(GeometryIndex, GeometryIndex),
    /// Two members touch along a line (1-dimensional intersection).
    ElementsTouchOnALine(GeometryIndex, GeometryIndex),
}

impl fmt::Display for InvalidMultiPolygon {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidMultiPolygon::InvalidPolygon(idx, err) => {
                write!(f, "polygon at index {} is invalid: {}", idx.0, err)
            }
            InvalidMultiPolygon::ElementsOverlaps(idx1, idx2) => {
                write!(f, "polygons at indices {} and {} overlap", idx1.0, idx2.0)
            }
            InvalidMultiPolygon::ElementsTouchOnALine(idx1, idx2) => {
                write!(
                    f,
                    "polygons at indices {} and {} touch on a line",
                    idx1.0, idx2.0
                )
            }
        }
    }
}

impl std::error::Error for InvalidMultiPolygon {}

impl Validation for MultiPolygon {
    type Error = InvalidMultiPolygon;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        for (i, polygon) in self.iter().enumerate() {
            polygon.visit_validation(Box::new(&mut |err| {
                handle_validation_error(InvalidMultiPolygon::InvalidPolygon(GeometryIndex(i), err))
            }))?;

            let polygon_geom = Geometry::Polygon(polygon.clone());
            for (j, polygon_2) in self.iter().enumerate().skip(i + 1) {
                let im = relate(&polygon_geom, &Geometry::Polygon(polygon_2.clone()));
                if im.interior_interior_dim() == 2 {
                    handle_validation_error(InvalidMultiPolygon::ElementsOverlaps(
                        GeometryIndex(i),
                        GeometryIndex(j),
                    ))?;
                }
                if im.boundary_boundary_dim() == 1 {
                    handle_validation_error(InvalidMultiPolygon::ElementsTouchOnALine(
                        GeometryIndex(i),
                        GeometryIndex(j),
                    ))?;
                }
            }
        }
        Ok(())
    }
}

/// Why a [`Rect`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidRect {
    /// A valid [`Rect`] must have finite corners; index 0 is `min`, index 1 `max`.
    NonFiniteCoord(CoordIndex),
}

impl fmt::Display for InvalidRect {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidRect::NonFiniteCoord(idx) => {
                let corner = if idx.0 == 0 { "min" } else { "max" };
                write!(f, "coordinate at rect's {corner} is non-finite")
            }
        }
    }
}

impl std::error::Error for InvalidRect {}

impl Validation for Rect {
    type Error = InvalidRect;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        if !coord_is_finite(&self.min()) {
            handle_validation_error(InvalidRect::NonFiniteCoord(CoordIndex(0)))?;
        }
        if !coord_is_finite(&self.max()) {
            handle_validation_error(InvalidRect::NonFiniteCoord(CoordIndex(1)))?;
        }
        Ok(())
    }
}

/// Why a [`Triangle`] is invalid.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidTriangle {
    /// A valid [`Triangle`] must have finite coordinates.
    NonFiniteCoord(CoordIndex),
    /// A valid [`Triangle`] must have distinct vertices.
    IdenticalCoords(CoordIndex, CoordIndex),
    /// A valid [`Triangle`] must have non-collinear vertices.
    CollinearCoords,
}

impl fmt::Display for InvalidTriangle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidTriangle::NonFiniteCoord(idx) => {
                write!(f, "coordinate at index {} is non-finite", idx.0)
            }
            InvalidTriangle::IdenticalCoords(idx1, idx2) => {
                write!(
                    f,
                    "coordinates at indices {} and {} are identical",
                    idx1.0, idx2.0
                )
            }
            InvalidTriangle::CollinearCoords => write!(f, "triangle has collinear coordinates"),
        }
    }
}

impl std::error::Error for InvalidTriangle {}

impl Validation for Triangle {
    type Error = InvalidTriangle;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        for (idx, v) in self.to_array().iter().enumerate() {
            if !coord_is_finite(v) {
                handle_validation_error(InvalidTriangle::NonFiniteCoord(CoordIndex(idx)))?;
            }
        }

        let mut identical = false;
        let [v1, v2, v3] = self.to_array();
        for (i, j, a, b) in [(0, 1, v1, v2), (0, 2, v1, v3), (1, 2, v2, v3)] {
            if a == b {
                handle_validation_error(InvalidTriangle::IdenticalCoords(
                    CoordIndex(i),
                    CoordIndex(j),
                ))?;
                identical = true;
            }
        }
        if !identical && orient(v1, v2, v3) == 0.0 {
            handle_validation_error(InvalidTriangle::CollinearCoords)?;
        }
        Ok(())
    }
}

/// Why a [`Geometry`] is invalid: the inner variant's reason, or for a
/// collection, which member failed and why.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidGeometry {
    InvalidPoint(InvalidPoint),
    InvalidLine(InvalidLine),
    InvalidLineString(InvalidLineString),
    InvalidPolygon(InvalidPolygon),
    InvalidMultiPoint(InvalidMultiPoint),
    InvalidMultiLineString(InvalidMultiLineString),
    InvalidMultiPolygon(InvalidMultiPolygon),
    InvalidRect(InvalidRect),
    InvalidTriangle(InvalidTriangle),
    /// Which member of a `GeometryCollection` is invalid, and why.
    InvalidGeometryCollection(GeometryIndex, Box<InvalidGeometry>),
}

impl fmt::Display for InvalidGeometry {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidGeometry::InvalidPoint(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidLine(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidLineString(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidPolygon(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidMultiPoint(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidMultiLineString(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidMultiPolygon(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidRect(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidTriangle(err) => write!(f, "{err}"),
            InvalidGeometry::InvalidGeometryCollection(idx, err) => {
                write!(f, "geometry at index {} is invalid: {}", idx.0, err)
            }
        }
    }
}

impl std::error::Error for InvalidGeometry {}

impl Validation for Geometry {
    type Error = InvalidGeometry;

    fn visit_validation<T>(
        &self,
        mut handle_validation_error: Box<dyn FnMut(Self::Error) -> Result<(), T> + '_>,
    ) -> Result<(), T> {
        match self {
            Geometry::Point(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidPoint(err))
            }))?,
            Geometry::Line(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidLine(err))
            }))?,
            Geometry::LineString(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidLineString(err))
            }))?,
            Geometry::Polygon(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidPolygon(err))
            }))?,
            Geometry::MultiPoint(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidMultiPoint(err))
            }))?,
            Geometry::MultiLineString(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidMultiLineString(err))
            }))?,
            Geometry::MultiPolygon(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidMultiPolygon(err))
            }))?,
            Geometry::Rect(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidRect(err))
            }))?,
            Geometry::Triangle(g) => g.visit_validation(Box::new(|err| {
                handle_validation_error(InvalidGeometry::InvalidTriangle(err))
            }))?,
            Geometry::GeometryCollection(gs) => {
                for (i, g) in gs.iter().enumerate() {
                    g.visit_validation(Box::new(&mut |err| {
                        handle_validation_error(InvalidGeometry::InvalidGeometryCollection(
                            GeometryIndex(i),
                            Box::new(err),
                        ))
                    }))?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ls(coords: &[(f64, f64)]) -> LineString {
        LineString::from(coords.to_vec())
    }

    #[test]
    fn simple_polygon_is_valid() {
        let polygon = Polygon::new(ls(&[(0.0, 0.0), (1.0, 1.0), (0.0, 1.0)]), vec![]);
        assert!(polygon.is_valid());
        assert!(polygon.validation_errors().is_empty());
    }

    #[test]
    fn bowtie_polygon_is_self_intersecting() {
        let polygon = Polygon::new(
            ls(&[(0.0, 0.0), (4.0, 0.0), (0.0, 2.0), (4.0, 2.0)]),
            vec![],
        );
        assert_eq!(
            polygon.check_validation().unwrap_err(),
            InvalidPolygon::SelfIntersection(RingRole::Exterior)
        );
    }

    #[test]
    fn spike_is_self_intersecting() {
        let polygon = Polygon::new(
            ls(&[
                (0.0, 0.0),
                (4.0, 0.0),
                (4.0, 4.0),
                (2.0, 4.0),
                (2.0, 6.0),
                (2.0, 4.0),
                (0.0, 4.0),
            ]),
            vec![],
        );
        assert_eq!(
            polygon.validation_errors(),
            vec![InvalidPolygon::SelfIntersection(RingRole::Exterior)]
        );
    }

    #[test]
    fn too_few_points_in_ring() {
        let polygon = Polygon::new(ls(&[(0.0, 0.0), (1.0, 1.0)]), vec![]);
        assert_eq!(
            polygon.validation_errors(),
            vec![InvalidPolygon::TooFewPointsInRing(RingRole::Exterior)]
        );
    }

    #[test]
    fn hole_outside_shell_not_contained() {
        let polygon = Polygon::new(
            ls(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
            vec![ls(&[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0)])],
        );
        assert_eq!(
            polygon.check_validation().unwrap_err(),
            InvalidPolygon::InteriorRingNotContainedInExteriorRing(RingRole::Interior(0))
        );
    }

    #[test]
    fn holes_touching_at_point_are_valid_on_a_line_invalid() {
        let shell = ls(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
        let point_touch = Polygon::new(
            shell.clone(),
            vec![
                ls(&[(1.0, 2.0), (2.0, 1.0), (3.0, 2.0), (2.0, 3.0)]),
                ls(&[(3.0, 2.0), (4.0, 1.0), (5.0, 2.0), (4.0, 3.0)]),
            ],
        );
        assert!(point_touch.is_valid());

        let line_touch = Polygon::new(
            shell,
            vec![
                ls(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]),
                ls(&[(3.0, 1.0), (5.0, 1.0), (5.0, 3.0), (3.0, 3.0)]),
            ],
        );
        assert_eq!(
            line_touch.check_validation().unwrap_err(),
            InvalidPolygon::IntersectingRingsOnALine(RingRole::Interior(0), RingRole::Interior(1))
        );
    }

    #[test]
    fn nested_holes_intersect_on_an_area() {
        let polygon = Polygon::new(
            ls(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
            vec![
                ls(&[(1.0, 1.0), (1.0, 9.0), (9.0, 9.0), (9.0, 1.0)]),
                ls(&[(2.0, 2.0), (2.0, 8.0), (8.0, 8.0), (8.0, 2.0)]),
            ],
        );
        assert_eq!(
            polygon.check_validation().unwrap_err(),
            InvalidPolygon::IntersectingRingsOnAnArea(RingRole::Interior(0), RingRole::Interior(1))
        );
    }

    #[test]
    fn nan_coordinate_is_reported_with_its_index() {
        let line = ls(&[(0.0, 0.0), (f64::NAN, 1.0), (2.0, 2.0)]);
        assert_eq!(
            line.validation_errors(),
            vec![InvalidLineString::NonFiniteCoord(CoordIndex(1))]
        );
    }

    #[test]
    fn multipolygon_overlap_and_touch() {
        let a = Polygon::new(
            ls(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
            vec![],
        );
        let overlapping = Polygon::new(
            ls(&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]),
            vec![],
        );
        let mp = MultiPolygon::new(vec![a.clone(), overlapping]);
        assert_eq!(
            mp.check_validation().unwrap_err(),
            InvalidMultiPolygon::ElementsOverlaps(GeometryIndex(0), GeometryIndex(1))
        );

        let edge_sharing = Polygon::new(
            ls(&[(4.0, 0.0), (8.0, 0.0), (8.0, 4.0), (4.0, 4.0)]),
            vec![],
        );
        let mp = MultiPolygon::new(vec![a.clone(), edge_sharing]);
        assert_eq!(
            mp.check_validation().unwrap_err(),
            InvalidMultiPolygon::ElementsTouchOnALine(GeometryIndex(0), GeometryIndex(1))
        );

        let disjoint = Polygon::new(
            ls(&[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0)]),
            vec![],
        );
        assert!(MultiPolygon::new(vec![a, disjoint]).is_valid());
    }

    #[test]
    fn geometry_collection_reports_member_index() {
        let g = Geometry::GeometryCollection(vec![
            Geometry::Point(Point::new(0.0, 0.0)),
            Geometry::LineString(ls(&[(0.0, 0.0), (0.0, 0.0)])),
        ]);
        assert_eq!(
            g.validation_errors(),
            vec![InvalidGeometry::InvalidGeometryCollection(
                GeometryIndex(1),
                Box::new(InvalidGeometry::InvalidLineString(
                    InvalidLineString::TooFewPoints
                )),
            )]
        );
    }
}
