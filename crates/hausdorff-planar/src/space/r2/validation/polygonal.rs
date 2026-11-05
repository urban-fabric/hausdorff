use std::fmt;

use super::*;
use crate::space::r2::HasDimensions;

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
