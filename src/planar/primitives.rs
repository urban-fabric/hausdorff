//! Small fixed-vertex planar geometries: [`Line`] (a segment), [`Triangle`],
//! and the axis-aligned [`Rect`]. Each can expand to a [`Polygon`] for the
//! algorithms that only speak in rings.

use serde::{Deserialize, Serialize};

use super::{Coord, LineString, Point, Polygon};

/// A line segment from `start` to `end`.
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Line {
    pub start: Coord,
    pub end: Coord,
}

impl Line {
    /// The segment from `start` to `end`.
    #[must_use]
    pub fn new(start: impl Into<Coord>, end: impl Into<Coord>) -> Self {
        Self {
            start: start.into(),
            end: end.into(),
        }
    }

    /// `end - start`.
    #[inline]
    #[must_use]
    pub fn delta(&self) -> Coord {
        self.end - self.start
    }

    /// Squared Euclidean length — cheaper than [`length`](Self::length) when you
    /// only need to compare.
    #[inline]
    #[must_use]
    pub fn length_squared(&self) -> f64 {
        let d = self.delta();
        d.x * d.x + d.y * d.y
    }

    /// Euclidean length.
    #[inline]
    #[must_use]
    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }
}

/// A triangle defined by its three vertices.
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Triangle(pub Coord, pub Coord, pub Coord);

impl Triangle {
    /// The triangle with the given vertices.
    #[must_use]
    pub fn new(a: Coord, b: Coord, c: Coord) -> Self {
        Self(a, b, c)
    }

    /// The three vertices as an array.
    #[inline]
    #[must_use]
    pub fn to_array(self) -> [Coord; 3] {
        [self.0, self.1, self.2]
    }

    /// This triangle as a closed three-sided [`Polygon`].
    #[must_use]
    pub fn to_polygon(self) -> Polygon {
        Polygon::new(LineString::new(vec![self.0, self.1, self.2]), vec![])
    }
}

/// An axis-aligned rectangle, stored by its `min` (lowest x, y) and `max`
/// (highest x, y) corners. [`new`](Self::new) accepts the corners in any order
/// and normalizes them.
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Rect {
    min: Coord,
    max: Coord,
}

impl Rect {
    /// The rectangle spanning two opposite corners, given in any order.
    #[must_use]
    pub fn new(c1: impl Into<Coord>, c2: impl Into<Coord>) -> Self {
        let (a, b) = (c1.into(), c2.into());
        Self {
            min: Coord::new(a.x.min(b.x), a.y.min(b.y)),
            max: Coord::new(a.x.max(b.x), a.y.max(b.y)),
        }
    }

    /// The lower corner (minimum x and y).
    #[inline]
    #[must_use]
    pub fn min(self) -> Coord {
        self.min
    }

    /// The upper corner (maximum x and y).
    #[inline]
    #[must_use]
    pub fn max(self) -> Coord {
        self.max
    }

    /// Width (extent along x).
    #[inline]
    #[must_use]
    pub fn width(self) -> f64 {
        self.max.x - self.min.x
    }

    /// Height (extent along y).
    #[inline]
    #[must_use]
    pub fn height(self) -> f64 {
        self.max.y - self.min.y
    }

    /// The center point.
    #[inline]
    #[must_use]
    pub fn center(self) -> Point {
        Point::new(
            (self.min.x + self.max.x) * 0.5,
            (self.min.y + self.max.y) * 0.5,
        )
    }

    /// True iff `coord` lies within the rectangle (inclusive).
    #[inline]
    #[must_use]
    pub fn contains_coord(self, coord: Coord) -> bool {
        coord.x >= self.min.x
            && coord.x <= self.max.x
            && coord.y >= self.min.y
            && coord.y <= self.max.y
    }

    /// This rectangle as a closed, counter-clockwise [`Polygon`].
    #[must_use]
    pub fn to_polygon(self) -> Polygon {
        Polygon::new(
            LineString::new(vec![
                Coord::new(self.min.x, self.min.y),
                Coord::new(self.max.x, self.min.y),
                Coord::new(self.max.x, self.max.y),
                Coord::new(self.min.x, self.max.y),
            ]),
            vec![],
        )
    }
}
