//! Topological [`Dimensions`] of a geometry and the [`HasDimensions`] trait,
//! mirroring `geo`'s `dimensions` module (the `dimensions()` / `is_empty()`
//! surface the workspace uses).

use super::{
    Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect,
    Triangle,
};

/// The topological dimension of a geometry.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Dimensions {
    /// No coordinates.
    Empty,
    /// Points (0-D).
    ZeroDimensional,
    /// Curves (1-D).
    OneDimensional,
    /// Surfaces (2-D).
    TwoDimensional,
}

/// Report a geometry's [`Dimensions`] and emptiness.
pub trait HasDimensions {
    /// The geometry's topological dimension.
    fn dimensions(&self) -> Dimensions;

    /// True iff the geometry holds no coordinates.
    fn is_empty(&self) -> bool {
        self.dimensions() == Dimensions::Empty
    }
}

impl HasDimensions for Point {
    fn dimensions(&self) -> Dimensions {
        Dimensions::ZeroDimensional
    }
}

impl HasDimensions for MultiPoint {
    fn dimensions(&self) -> Dimensions {
        if self.is_empty() {
            Dimensions::Empty
        } else {
            Dimensions::ZeroDimensional
        }
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl HasDimensions for Line {
    fn dimensions(&self) -> Dimensions {
        if self.start == self.end {
            Dimensions::ZeroDimensional
        } else {
            Dimensions::OneDimensional
        }
    }
}

impl HasDimensions for LineString {
    fn dimensions(&self) -> Dimensions {
        if self.0.is_empty() {
            Dimensions::Empty
        } else if self.0.iter().all(|c| *c == self.0[0]) {
            Dimensions::ZeroDimensional
        } else {
            Dimensions::OneDimensional
        }
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl HasDimensions for MultiLineString {
    fn dimensions(&self) -> Dimensions {
        self.iter()
            .map(HasDimensions::dimensions)
            .max_by_key(rank)
            .unwrap_or(Dimensions::Empty)
    }

    fn is_empty(&self) -> bool {
        self.iter().all(HasDimensions::is_empty)
    }
}

impl HasDimensions for Polygon {
    fn dimensions(&self) -> Dimensions {
        if self.exterior().0.is_empty() {
            Dimensions::Empty
        } else {
            Dimensions::TwoDimensional
        }
    }

    fn is_empty(&self) -> bool {
        self.exterior().0.is_empty()
    }
}

impl HasDimensions for MultiPolygon {
    fn dimensions(&self) -> Dimensions {
        if self.is_empty() {
            Dimensions::Empty
        } else {
            Dimensions::TwoDimensional
        }
    }

    fn is_empty(&self) -> bool {
        self.iter().all(HasDimensions::is_empty)
    }
}

impl HasDimensions for Rect {
    fn dimensions(&self) -> Dimensions {
        Dimensions::TwoDimensional
    }
}

impl HasDimensions for Triangle {
    fn dimensions(&self) -> Dimensions {
        Dimensions::TwoDimensional
    }
}

impl HasDimensions for Geometry {
    fn dimensions(&self) -> Dimensions {
        match self {
            Geometry::Point(g) => g.dimensions(),
            Geometry::Line(g) => g.dimensions(),
            Geometry::LineString(g) => g.dimensions(),
            Geometry::Polygon(g) => g.dimensions(),
            Geometry::MultiPoint(g) => g.dimensions(),
            Geometry::MultiLineString(g) => g.dimensions(),
            Geometry::MultiPolygon(g) => g.dimensions(),
            Geometry::Rect(g) => g.dimensions(),
            Geometry::Triangle(g) => g.dimensions(),
            Geometry::GeometryCollection(gs) => gs
                .iter()
                .map(HasDimensions::dimensions)
                .max_by_key(rank)
                .unwrap_or(Dimensions::Empty),
        }
    }
}

fn rank(d: &Dimensions) -> u8 {
    match d {
        Dimensions::Empty => 0,
        Dimensions::ZeroDimensional => 1,
        Dimensions::OneDimensional => 2,
        Dimensions::TwoDimensional => 3,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planar::LineString;

    #[test]
    fn dimensions_by_type() {
        assert_eq!(Point::new(0.0, 0.0).dimensions(), Dimensions::ZeroDimensional);
        assert_eq!(
            LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]).dimensions(),
            Dimensions::OneDimensional
        );
        assert_eq!(LineString::default().dimensions(), Dimensions::Empty);
        assert!(LineString::default().is_empty());
        assert_eq!(
            Polygon::new(LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]), vec![])
                .dimensions(),
            Dimensions::TwoDimensional
        );
    }
}
