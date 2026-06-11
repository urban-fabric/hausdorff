//! [`MapCoords`] — apply a coordinate-wise function to a geometry, producing a
//! geometry of the same shape. Mirrors `geo`'s `MapCoords`.
//!
//! This is the generic "transform every vertex" hook: translations, scalings,
//! axis swaps, and — via [`crate::earth::projection::Projection::project_geometry`]
//! — map-projection of whole geometries.

use super::{
    Coord, Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    Rect, Triangle,
};

/// Map a function over every coordinate of a geometry, returning a geometry of
/// the same type. The function is applied to each vertex independently;
/// topology (ring closure, part structure) is preserved.
pub trait MapCoords: Sized {
    /// This geometry with `f` applied to every coordinate.
    #[must_use]
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self;
}

impl MapCoords for Point {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        Point(f(self.0))
    }
}

impl MapCoords for Line {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        Line::new(f(self.start), f(self.end))
    }
}

impl MapCoords for LineString {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        LineString(self.0.iter().map(|&c| f(c)).collect())
    }
}

/// Rings stay rings: a closed ring maps vertex-by-vertex to a closed ring
/// (`Polygon::new` re-closes defensively).
impl MapCoords for Polygon {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        Polygon::new(
            self.exterior().map_coords(f),
            self.interiors().iter().map(|r| r.map_coords(f)).collect(),
        )
    }
}

impl MapCoords for MultiPoint {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        MultiPoint(self.0.iter().map(|p| p.map_coords(f)).collect())
    }
}

impl MapCoords for MultiLineString {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        MultiLineString(self.0.iter().map(|ls| ls.map_coords(f)).collect())
    }
}

impl MapCoords for MultiPolygon {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        MultiPolygon(self.0.iter().map(|p| p.map_coords(f)).collect())
    }
}

/// Only the two corners are mapped and the result re-normalized, so the result
/// is the axis-aligned rectangle spanning the mapped corners — exact for
/// rectilinear maps (translation, scaling, axis swap), an approximation for
/// anything that bends lines (e.g. a map projection). Convert
/// [`to_polygon`](Rect::to_polygon) first when the four edges must be mapped
/// faithfully.
impl MapCoords for Rect {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        Rect::new(f(self.min()), f(self.max()))
    }
}

impl MapCoords for Triangle {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        Triangle::new(f(self.0), f(self.1), f(self.2))
    }
}

impl MapCoords for Geometry {
    fn map_coords(&self, f: impl Fn(Coord) -> Coord + Copy) -> Self {
        match self {
            Geometry::Point(p) => Geometry::Point(p.map_coords(f)),
            Geometry::Line(l) => Geometry::Line(l.map_coords(f)),
            Geometry::LineString(ls) => Geometry::LineString(ls.map_coords(f)),
            Geometry::Polygon(p) => Geometry::Polygon(p.map_coords(f)),
            Geometry::MultiPoint(mp) => Geometry::MultiPoint(mp.map_coords(f)),
            Geometry::MultiLineString(mls) => Geometry::MultiLineString(mls.map_coords(f)),
            Geometry::MultiPolygon(mp) => Geometry::MultiPolygon(mp.map_coords(f)),
            Geometry::Rect(r) => Geometry::Rect(r.map_coords(f)),
            Geometry::Triangle(t) => Geometry::Triangle(t.map_coords(f)),
            Geometry::GeometryCollection(gs) => {
                Geometry::GeometryCollection(gs.iter().map(|g| g.map_coords(f)).collect())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shift(c: Coord) -> Coord {
        Coord::new(c.x + 1.0, c.y - 2.0)
    }

    #[test]
    fn polygon_rings_map_and_stay_closed() {
        let poly = Polygon::new(
            LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
            vec![LineString::from(vec![
                (1.0, 1.0),
                (2.0, 1.0),
                (2.0, 2.0),
                (1.0, 2.0),
            ])],
        );
        let mapped = poly.map_coords(shift);
        assert!(mapped.exterior().is_closed());
        assert_eq!(mapped.exterior().0[0], Coord::new(1.0, -2.0));
        assert_eq!(mapped.interiors().len(), 1);
        assert_eq!(mapped.interiors()[0].0[1], Coord::new(3.0, -1.0));
    }

    #[test]
    fn geometry_collection_recurses() {
        let g = Geometry::GeometryCollection(vec![
            Geometry::Point(Point::new(1.0, 1.0)),
            Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0)])),
        ]);
        let Geometry::GeometryCollection(parts) = g.map_coords(shift) else {
            panic!("collection should map to collection");
        };
        assert_eq!(parts[0], Geometry::Point(Point::new(2.0, -1.0)));
    }

    #[test]
    fn rect_renormalizes_swapped_corners() {
        let r = Rect::new(Coord::new(0.0, 0.0), Coord::new(2.0, 3.0));
        // A reflection swaps which corner is min/max; Rect::new re-normalizes.
        let mapped = r.map_coords(|c| Coord::new(-c.x, -c.y));
        assert_eq!(mapped.min(), Coord::new(-2.0, -3.0));
        assert_eq!(mapped.max(), Coord::new(0.0, 0.0));
    }
}
