use std::fmt;

use super::super::util::orient;
use super::*;

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
