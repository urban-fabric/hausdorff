//! The [`Geometry`] sum type — any one of the planar geometries, for code that
//! handles heterogeneous shapes (e.g. a feature whose geometry varies at
//! runtime). Mirrors `geo_types::Geometry`.

use serde::{Deserialize, Serialize};

use super::{
    Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect, Triangle,
};

/// Any planar geometry.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub enum Geometry {
    Point(Point),
    Line(Line),
    LineString(LineString),
    Polygon(Polygon),
    MultiPoint(MultiPoint),
    MultiLineString(MultiLineString),
    MultiPolygon(MultiPolygon),
    Rect(Rect),
    Triangle(Triangle),
    GeometryCollection(Vec<Geometry>),
}

macro_rules! from_variant {
    ($variant:ident, $ty:ty) => {
        impl From<$ty> for Geometry {
            fn from(g: $ty) -> Self {
                Geometry::$variant(g)
            }
        }
    };
}

from_variant!(Point, Point);
from_variant!(Line, Line);
from_variant!(LineString, LineString);
from_variant!(Polygon, Polygon);
from_variant!(MultiPoint, MultiPoint);
from_variant!(MultiLineString, MultiLineString);
from_variant!(MultiPolygon, MultiPolygon);
from_variant!(Rect, Rect);
from_variant!(Triangle, Triangle);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_conversions() {
        let g: Geometry = Point::new(1.0, 2.0).into();
        assert!(matches!(g, Geometry::Point(_)));
        let g: Geometry = Polygon::default().into();
        assert!(matches!(g, Geometry::Polygon(_)));
    }
}
