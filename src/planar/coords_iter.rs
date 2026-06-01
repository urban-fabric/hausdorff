//! Iterate the coordinates of any planar geometry — mirror of `geo`'s
//! `CoordsIter` (`coords_iter` / `coords_count` / `exterior_coords_iter`).
//!
//! For simplicity the iterators are materialized (`Vec` → `IntoIter`); the
//! counts are computed structurally without allocating.

use super::{
    Coord, Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    Rect, Triangle,
};

/// Walk a geometry's coordinates.
pub trait CoordsIter {
    /// Every coordinate of the geometry, exterior and interior alike.
    fn coords_iter(&self) -> std::vec::IntoIter<Coord> {
        self.to_coords().into_iter()
    }

    /// Collect every coordinate into a `Vec` (the backing of [`coords_iter`]).
    ///
    /// [`coords_iter`]: CoordsIter::coords_iter
    fn to_coords(&self) -> Vec<Coord>;

    /// Number of coordinates, without allocating.
    fn coords_count(&self) -> usize;

    /// Only the exterior coordinates (for a polygon, its shell; for most other
    /// geometries, identical to [`coords_iter`](CoordsIter::coords_iter)).
    fn exterior_coords_iter(&self) -> std::vec::IntoIter<Coord> {
        self.to_coords().into_iter()
    }
}

impl CoordsIter for Coord {
    fn to_coords(&self) -> Vec<Coord> {
        vec![*self]
    }
    fn coords_count(&self) -> usize {
        1
    }
}

impl CoordsIter for Point {
    fn to_coords(&self) -> Vec<Coord> {
        vec![self.0]
    }
    fn coords_count(&self) -> usize {
        1
    }
}

impl CoordsIter for Line {
    fn to_coords(&self) -> Vec<Coord> {
        vec![self.start, self.end]
    }
    fn coords_count(&self) -> usize {
        2
    }
}

impl CoordsIter for Triangle {
    fn to_coords(&self) -> Vec<Coord> {
        vec![self.0, self.1, self.2]
    }
    fn coords_count(&self) -> usize {
        3
    }
}

impl CoordsIter for Rect {
    fn to_coords(&self) -> Vec<Coord> {
        self.to_polygon().to_coords()
    }
    fn coords_count(&self) -> usize {
        5
    }
}

impl CoordsIter for LineString {
    fn to_coords(&self) -> Vec<Coord> {
        self.0.clone()
    }
    fn coords_count(&self) -> usize {
        self.0.len()
    }
}

impl CoordsIter for MultiLineString {
    fn to_coords(&self) -> Vec<Coord> {
        self.iter().flat_map(|ls| ls.0.iter().copied()).collect()
    }
    fn coords_count(&self) -> usize {
        self.iter().map(LineString::coords_count).sum()
    }
}

impl CoordsIter for Polygon {
    fn to_coords(&self) -> Vec<Coord> {
        self.exterior()
            .0
            .iter()
            .chain(self.interiors().iter().flat_map(|r| r.0.iter()))
            .copied()
            .collect()
    }
    fn coords_count(&self) -> usize {
        self.exterior().coords_count()
            + self.interiors().iter().map(LineString::coords_count).sum::<usize>()
    }
    fn exterior_coords_iter(&self) -> std::vec::IntoIter<Coord> {
        self.exterior().0.clone().into_iter()
    }
}

impl CoordsIter for MultiPolygon {
    fn to_coords(&self) -> Vec<Coord> {
        self.iter().flat_map(Polygon::to_coords).collect()
    }
    fn coords_count(&self) -> usize {
        self.iter().map(Polygon::coords_count).sum()
    }
    fn exterior_coords_iter(&self) -> std::vec::IntoIter<Coord> {
        self.iter()
            .flat_map(|p| p.exterior().0.clone())
            .collect::<Vec<_>>()
            .into_iter()
    }
}

impl CoordsIter for MultiPoint {
    fn to_coords(&self) -> Vec<Coord> {
        self.iter().map(|p| p.0).collect()
    }
    fn coords_count(&self) -> usize {
        self.0.len()
    }
}

impl CoordsIter for Geometry {
    fn to_coords(&self) -> Vec<Coord> {
        match self {
            Geometry::Point(g) => g.to_coords(),
            Geometry::Line(g) => g.to_coords(),
            Geometry::LineString(g) => g.to_coords(),
            Geometry::Polygon(g) => g.to_coords(),
            Geometry::MultiPoint(g) => g.to_coords(),
            Geometry::MultiLineString(g) => g.to_coords(),
            Geometry::MultiPolygon(g) => g.to_coords(),
            Geometry::Rect(g) => g.to_coords(),
            Geometry::Triangle(g) => g.to_coords(),
            Geometry::GeometryCollection(gs) => gs.iter().flat_map(Geometry::to_coords).collect(),
        }
    }
    fn coords_count(&self) -> usize {
        match self {
            Geometry::Point(g) => g.coords_count(),
            Geometry::Line(g) => g.coords_count(),
            Geometry::LineString(g) => g.coords_count(),
            Geometry::Polygon(g) => g.coords_count(),
            Geometry::MultiPoint(g) => g.coords_count(),
            Geometry::MultiLineString(g) => g.coords_count(),
            Geometry::MultiPolygon(g) => g.coords_count(),
            Geometry::Rect(g) => g.coords_count(),
            Geometry::Triangle(g) => g.coords_count(),
            Geometry::GeometryCollection(gs) => gs.iter().map(Geometry::coords_count).sum(),
        }
    }
}
