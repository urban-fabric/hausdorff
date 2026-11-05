//! Area geometries: [`Polygon`] (an exterior ring with optional holes),
//! [`MultiPolygon`], and the point collection [`MultiPoint`].

use serde::{Deserialize, Serialize};

use hausdorff_traits::is_convex_ring;

use crate::space::r2::util::ring_vertices;
use crate::space::r2::{LineString, Point};

/// A polygon: one exterior ring plus zero or more interior rings (holes). Rings
/// are stored closed — [`new`](Self::new) closes any that are not.
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct Polygon {
    exterior: LineString,
    interiors: Vec<LineString>,
}

impl Polygon {
    /// A polygon with the given exterior ring and holes; all rings are closed.
    #[must_use]
    pub fn new(mut exterior: LineString, mut interiors: Vec<LineString>) -> Self {
        exterior.close();
        for ring in &mut interiors {
            ring.close();
        }
        Self {
            exterior,
            interiors,
        }
    }

    /// The exterior ring.
    #[inline]
    #[must_use]
    pub fn exterior(&self) -> &LineString {
        &self.exterior
    }

    /// The interior rings (holes).
    #[inline]
    #[must_use]
    pub fn interiors(&self) -> &[LineString] {
        &self.interiors
    }

    /// Mutate the exterior ring in place, re-closing it afterwards.
    pub fn exterior_mut(&mut self, f: impl FnOnce(&mut LineString)) {
        f(&mut self.exterior);
        self.exterior.close();
    }

    /// Decompose into `(exterior, interiors)`.
    #[must_use]
    pub fn into_inner(self) -> (LineString, Vec<LineString>) {
        (self.exterior, self.interiors)
    }

    /// True iff this polygon is **convex** — it has no holes and its exterior ring
    /// is a convex polygon (every vertex turns the same way). A convex polygon
    /// that *happens* to be one (e.g. a triangle or box read from data) is
    /// detected here so the convex fast paths apply; for a compile-time guarantee
    /// and the convex-only algorithms (O(log n) containment), upgrade with
    /// [`ConvexPolygon::from_polygon`](crate::space::r2::ConvexPolygon::from_polygon).
    ///
    /// The consistent-turn-sign test is the generic [`Orientation`] convexity
    /// classifier, proof-backed by `riemann`
    /// `Geometry/Orientation.OrientGeom.leftOfAll_convex` (planar
    /// `leftOfAllStrict_convex`): a consistently-wound ring's left-of-every-edge
    /// region is convex.
    ///
    /// [`Orientation`]: crate::algebra::Orientation
    #[must_use]
    pub fn is_convex(&self) -> bool {
        self.interiors.is_empty() && is_convex_ring(ring_vertices(&self.exterior.0))
    }
}

/// A collection of [`Polygon`]s.
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct MultiPolygon(pub Vec<Polygon>);

impl MultiPolygon {
    /// A multi-polygon over `polygons`.
    #[must_use]
    pub fn new(polygons: Vec<Polygon>) -> Self {
        Self(polygons)
    }

    /// The component polygons.
    pub fn iter(&self) -> impl Iterator<Item = &Polygon> {
        self.0.iter()
    }

    /// Number of component polygons.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// True iff there are no component polygons.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// A collection of [`Point`]s.
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct MultiPoint(pub Vec<Point>);

impl MultiPoint {
    /// A multi-point over `points`.
    #[must_use]
    pub fn new(points: Vec<Point>) -> Self {
        Self(points)
    }

    /// The component points.
    pub fn iter(&self) -> impl Iterator<Item = &Point> {
        self.0.iter()
    }

    /// Number of points.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// True iff there are no points.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl From<Vec<Point>> for MultiPoint {
    fn from(points: Vec<Point>) -> Self {
        Self(points)
    }
}

impl FromIterator<Point> for MultiPoint {
    fn from_iter<I: IntoIterator<Item = Point>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<'a> IntoIterator for &'a MultiPoint {
    type Item = &'a Point;
    type IntoIter = std::slice::Iter<'a, Point>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a MultiPolygon {
    type Item = &'a Polygon;
    type IntoIter = std::slice::Iter<'a, Polygon>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl From<Vec<(f64, f64)>> for MultiPoint {
    fn from(coords: Vec<(f64, f64)>) -> Self {
        Self(coords.into_iter().map(Point::from).collect())
    }
}
