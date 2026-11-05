//! Closest point on a geometry to a query point — mirror of `geo`'s
//! `ClosestPoint` trait and its `Closest` result enum.
//!
//! The answer distinguishes the point lying *on* the geometry
//! ([`Closest::Intersection`]), a unique nearest location
//! ([`Closest::SinglePoint`]), and the degenerate cases (empty or zero-length
//! geometry) where no single answer exists ([`Closest::Indeterminate`]).

use std::iter;

use crate::space::r2::util::point_in_ring;
use crate::space::r2::{
    Coord, Distance, Euclidean, Geometry, Line, LineString, MultiLineString, MultiPoint,
    MultiPolygon, Point, Polygon, Rect, Triangle,
};

/// The result of seeking the spot on a geometry closest to a point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Closest {
    /// The query point lies on the geometry itself.
    Intersection(Point),
    /// Exactly one place on the geometry is closest to the query point.
    SinglePoint(Point),
    /// Two or more (possibly infinitely many or undefined) places are equally
    /// close — e.g. the geometry is empty or a zero-length line.
    Indeterminate,
}

impl Closest {
    /// The better of two candidates relative to `p`: an intersection wins
    /// outright, `Indeterminate` loses to anything, and otherwise the nearer
    /// single point wins (`self` on ties).
    #[must_use]
    pub fn best_of_two(&self, other: &Self, p: Point) -> Self {
        let left = match *self {
            Closest::Indeterminate => return *other,
            Closest::Intersection(_) => return *self,
            Closest::SinglePoint(l) => l,
        };
        let right = match *other {
            Closest::Indeterminate => return *self,
            Closest::Intersection(_) => return *other,
            Closest::SinglePoint(r) => r,
        };
        if Euclidean.distance(left, p) <= Euclidean.distance(right, p) {
            *self
        } else {
            *other
        }
    }
}

/// Find the closest point between `self` and an input [`Point`].
pub trait ClosestPoint<Rhs = Point> {
    fn closest_point(&self, p: &Rhs) -> Closest;
}

impl<C: ClosestPoint> ClosestPoint for &C {
    fn closest_point(&self, p: &Point) -> Closest {
        (*self).closest_point(p)
    }
}

/// True iff `p` lies on the bounded segment `a`–`b` (collinear and within the
/// segment's bounding box).
fn on_segment(p: Coord, a: Coord, b: Coord) -> bool {
    let cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    cross == 0.0
        && p.x <= a.x.max(b.x)
        && p.x >= a.x.min(b.x)
        && p.y <= a.y.max(b.y)
        && p.y >= a.y.min(b.y)
}

/// True iff `c` lies on any segment of `ring`.
fn on_ring(ring: &LineString, c: Coord) -> bool {
    ring.lines().any(|l| on_segment(c, l.start, l.end))
}

/// Fold an iterator of sub-geometries down to the best `Closest`, short
/// circuiting on the first intersection (nothing can beat it).
fn closest_of<C, I>(iter: I, p: Point) -> Closest
where
    I: IntoIterator<Item = C>,
    C: ClosestPoint,
{
    let mut best = Closest::Indeterminate;
    for element in iter {
        best = element.closest_point(&p).best_of_two(&best, p);
        if matches!(best, Closest::Intersection(_)) {
            return best;
        }
    }
    best
}

impl ClosestPoint for Point {
    fn closest_point(&self, p: &Point) -> Closest {
        if self == p {
            Closest::Intersection(*self)
        } else {
            Closest::SinglePoint(*self)
        }
    }
}

impl ClosestPoint for Coord {
    fn closest_point(&self, p: &Point) -> Closest {
        Point(*self).closest_point(p)
    }
}

impl ClosestPoint for Line {
    fn closest_point(&self, p: &Point) -> Closest {
        if self.length() == 0.0 {
            return Closest::Indeterminate;
        }
        let direction = self.delta();
        let to_p = p.0 - self.start;
        let t = to_p.dot(direction) / direction.dot(direction);
        if t < 0.0 {
            return Closest::SinglePoint(Point(self.start));
        }
        if t > 1.0 {
            return Closest::SinglePoint(Point(self.end));
        }
        let c = Point(self.start + direction * t);
        if on_segment(p.0, self.start, self.end) {
            Closest::Intersection(c)
        } else {
            Closest::SinglePoint(c)
        }
    }
}

impl ClosestPoint for LineString {
    fn closest_point(&self, p: &Point) -> Closest {
        closest_of(self.lines(), *p)
    }
}

impl ClosestPoint for Polygon {
    fn closest_point(&self, p: &Point) -> Closest {
        let interior = point_in_ring(p.0, &self.exterior().0)
            && !self.interiors().iter().any(|h| point_in_ring(p.0, &h.0));
        let boundary = iter::once(self.exterior())
            .chain(self.interiors().iter())
            .any(|r| on_ring(r, p.0));
        if interior || boundary {
            return Closest::Intersection(*p);
        }
        let prospectives = self.interiors().iter().chain(iter::once(self.exterior()));
        closest_of(prospectives, *p)
    }
}

impl ClosestPoint for Triangle {
    fn closest_point(&self, p: &Point) -> Closest {
        let ring = [self.0, self.1, self.2, self.0];
        let edges = [
            Line::new(self.0, self.1),
            Line::new(self.1, self.2),
            Line::new(self.2, self.0),
        ];
        if point_in_ring(p.0, &ring) || edges.iter().any(|l| on_segment(p.0, l.start, l.end)) {
            return Closest::Intersection(*p);
        }
        closest_of(edges, *p)
    }
}

impl ClosestPoint for Rect {
    fn closest_point(&self, p: &Point) -> Closest {
        if self.contains_coord(p.0) {
            return Closest::Intersection(*p);
        }
        let (min, max) = (self.min(), self.max());
        let edges = [
            Line::new(Coord::new(max.x, min.y), Coord::new(max.x, max.y)),
            Line::new(Coord::new(max.x, max.y), Coord::new(min.x, max.y)),
            Line::new(Coord::new(min.x, max.y), Coord::new(min.x, min.y)),
            Line::new(Coord::new(min.x, min.y), Coord::new(max.x, min.y)),
        ];
        closest_of(edges, *p)
    }
}

impl ClosestPoint for MultiPoint {
    fn closest_point(&self, p: &Point) -> Closest {
        closest_of(self.iter(), *p)
    }
}

impl ClosestPoint for MultiLineString {
    fn closest_point(&self, p: &Point) -> Closest {
        closest_of(self.iter(), *p)
    }
}

impl ClosestPoint for MultiPolygon {
    fn closest_point(&self, p: &Point) -> Closest {
        closest_of(self.iter(), *p)
    }
}

impl ClosestPoint for Geometry {
    fn closest_point(&self, p: &Point) -> Closest {
        match self {
            Geometry::Point(g) => g.closest_point(p),
            Geometry::Line(g) => g.closest_point(p),
            Geometry::LineString(g) => g.closest_point(p),
            Geometry::Polygon(g) => g.closest_point(p),
            Geometry::MultiPoint(g) => g.closest_point(p),
            Geometry::MultiLineString(g) => g.closest_point(p),
            Geometry::MultiPolygon(g) => g.closest_point(p),
            Geometry::Rect(g) => g.closest_point(p),
            Geometry::Triangle(g) => g.closest_point(p),
            Geometry::GeometryCollection(gs) => closest_of(gs.iter(), *p),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn diagonal() -> Line {
        Line::new(Coord::new(0.0, 0.0), Coord::new(100.0, 100.0))
    }

    #[test]
    fn point_to_point() {
        let a = Point::new(1.0, 2.0);
        let b = Point::new(3.0, 4.0);
        assert_eq!(a.closest_point(&a), Closest::Intersection(a));
        assert_eq!(a.closest_point(&b), Closest::SinglePoint(a));
    }

    #[test]
    fn line_endpoints_interior_and_projection() {
        let line = diagonal();
        assert_eq!(
            line.closest_point(&Point::new(-10.0, 0.0)),
            Closest::SinglePoint(Point::new(0.0, 0.0))
        );
        assert_eq!(
            line.closest_point(&Point::new(200.0, 200.0)),
            Closest::SinglePoint(Point::new(100.0, 100.0))
        );
        assert_eq!(
            line.closest_point(&Point::new(50.0, 50.0)),
            Closest::Intersection(Point::new(50.0, 50.0))
        );
        assert_eq!(
            line.closest_point(&Point::new(0.0, 100.0)),
            Closest::SinglePoint(Point::new(50.0, 50.0))
        );
    }

    #[test]
    fn zero_length_line_is_indeterminate() {
        let line = Line::new(Coord::new(1.0, 1.0), Coord::new(1.0, 1.0));
        assert_eq!(
            line.closest_point(&Point::new(0.0, 0.0)),
            Closest::Indeterminate
        );
    }

    #[test]
    fn empty_geometries_are_indeterminate() {
        let p = Point::new(1.0, 2.0);
        assert_eq!(
            LineString::new(vec![]).closest_point(&p),
            Closest::Indeterminate
        );
        assert_eq!(
            MultiPoint::new(vec![]).closest_point(&p),
            Closest::Indeterminate
        );
        assert_eq!(
            MultiPolygon::new(vec![]).closest_point(&p),
            Closest::Indeterminate
        );
        assert_eq!(
            Geometry::GeometryCollection(vec![]).closest_point(&p),
            Closest::Indeterminate
        );
    }

    #[test]
    fn polygon_interior_and_hole() {
        let poly = Polygon::new(
            LineString::from(vec![(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
            vec![LineString::from(vec![
                (4.0, 4.0),
                (6.0, 4.0),
                (6.0, 6.0),
                (4.0, 6.0),
            ])],
        );
        let inside = Point::new(2.0, 2.0);
        assert_eq!(poly.closest_point(&inside), Closest::Intersection(inside));
        let in_hole = Point::new(5.0, 4.5);
        assert_eq!(
            poly.closest_point(&in_hole),
            Closest::SinglePoint(Point::new(5.0, 4.0))
        );
        let outside = Point::new(15.0, 5.0);
        assert_eq!(
            poly.closest_point(&outside),
            Closest::SinglePoint(Point::new(10.0, 5.0))
        );
    }

    #[test]
    fn rect_and_triangle() {
        let rect = Rect::new(Coord::new(0.0, 0.0), Coord::new(4.0, 2.0));
        assert_eq!(
            rect.closest_point(&Point::new(2.0, 1.0)),
            Closest::Intersection(Point::new(2.0, 1.0))
        );
        assert_eq!(
            rect.closest_point(&Point::new(2.0, 5.0)),
            Closest::SinglePoint(Point::new(2.0, 2.0))
        );
        let tri = Triangle::new(
            Coord::new(0.0, 0.0),
            Coord::new(4.0, 0.0),
            Coord::new(0.0, 4.0),
        );
        assert_eq!(
            tri.closest_point(&Point::new(1.0, 1.0)),
            Closest::Intersection(Point::new(1.0, 1.0))
        );
        assert_eq!(
            tri.closest_point(&Point::new(2.0, -3.0)),
            Closest::SinglePoint(Point::new(2.0, 0.0))
        );
    }
}
