//! Polylines: [`LineString`] (an ordered run of coordinates) and
//! [`MultiLineString`] (a collection of them).

use serde::{Deserialize, Serialize};

use crate::space::r2::{Coord, Line, Point};

/// An ordered sequence of [`Coord`]s. A *ring* is a `LineString` whose first and
/// last coordinate are equal ([`is_closed`](Self::is_closed)); polygon
/// boundaries are stored that way.
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct LineString(pub Vec<Coord>);

impl LineString {
    /// A line string over `coords`.
    #[must_use]
    pub fn new(coords: Vec<Coord>) -> Self {
        Self(coords)
    }

    /// The coordinates, by reference.
    pub fn coords(&self) -> impl Iterator<Item = &Coord> {
        self.0.iter()
    }

    /// The coordinates, owned/copied.
    pub fn coords_copied(&self) -> impl Iterator<Item = Coord> + '_ {
        self.0.iter().copied()
    }

    /// Number of coordinates.
    #[inline]
    #[must_use]
    pub fn coords_count(&self) -> usize {
        self.0.len()
    }

    /// True iff there are no coordinates.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// The coordinates as [`Point`]s.
    pub fn points(&self) -> impl Iterator<Item = Point> + '_ {
        self.0.iter().map(|&c| Point(c))
    }

    /// The segments between consecutive coordinates.
    pub fn lines(&self) -> impl Iterator<Item = Line> + '_ {
        self.0.windows(2).map(|w| Line::new(w[0], w[1]))
    }

    /// True iff the string is non-empty and its first and last coordinate match
    /// — i.e. it is a closed ring.
    #[must_use]
    pub fn is_closed(&self) -> bool {
        match (self.0.first(), self.0.last()) {
            (Some(f), Some(l)) => f == l,
            _ => false,
        }
    }

    /// Close the ring in place: if non-empty and not already closed, push a copy
    /// of the first coordinate onto the end.
    pub fn close(&mut self) {
        if !self.is_closed() && !self.0.is_empty() {
            self.0.push(self.0[0]);
        }
    }
}

impl From<Vec<Coord>> for LineString {
    fn from(coords: Vec<Coord>) -> Self {
        Self(coords)
    }
}

impl From<Vec<(f64, f64)>> for LineString {
    fn from(coords: Vec<(f64, f64)>) -> Self {
        Self(coords.into_iter().map(Coord::from).collect())
    }
}

impl From<Vec<Point>> for LineString {
    fn from(points: Vec<Point>) -> Self {
        Self(points.into_iter().map(|p| p.0).collect())
    }
}

impl FromIterator<Coord> for LineString {
    fn from_iter<I: IntoIterator<Item = Coord>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

/// A collection of [`LineString`]s.
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct MultiLineString(pub Vec<LineString>);

impl MultiLineString {
    /// A multi-line-string over `line_strings`.
    #[must_use]
    pub fn new(line_strings: Vec<LineString>) -> Self {
        Self(line_strings)
    }

    /// The component line strings.
    pub fn iter(&self) -> impl Iterator<Item = &LineString> {
        self.0.iter()
    }

    /// Number of component line strings.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// True iff there are no component line strings.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}
