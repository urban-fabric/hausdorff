use std::fmt;

use super::*;

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
