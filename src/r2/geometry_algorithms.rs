//! Trait dispatch over the [`Geometry`] enum, so the planar algorithms run on a
//! heterogeneous geometry (not only the concrete types) — the native equivalent
//! of calling `geo`'s traits on a `geo::Geometry`.

use super::convex_hull::convex_hull_coords;
use super::{
    Area, BoundingRect, Centroid, ConvexHull, Coord, CoordsIter, Geometry, LineString, MultiPoint,
    Point, Polygon, Rect, Simplify,
};

impl Area for Geometry {
    fn signed_area(&self) -> f64 {
        match self {
            Geometry::Polygon(g) => g.signed_area(),
            Geometry::MultiPolygon(g) => g.signed_area(),
            Geometry::Rect(g) => g.signed_area(),
            Geometry::Triangle(g) => g.signed_area(),
            Geometry::GeometryCollection(gs) => gs.iter().map(Area::signed_area).sum(),
            // points and curves have no area
            Geometry::Point(_)
            | Geometry::Line(_)
            | Geometry::LineString(_)
            | Geometry::MultiPoint(_)
            | Geometry::MultiLineString(_) => 0.0,
        }
    }
}

impl BoundingRect for Geometry {
    type Output = Option<Rect>;

    fn bounding_rect(&self) -> Option<Rect> {
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        let mut any = false;
        for c in self.coords_iter() {
            if !c.is_finite() {
                continue;
            }
            any = true;
            min_x = min_x.min(c.x);
            min_y = min_y.min(c.y);
            max_x = max_x.max(c.x);
            max_y = max_y.max(c.y);
        }
        any.then(|| Rect::new(Coord::new(min_x, min_y), Coord::new(max_x, max_y)))
    }
}

impl ConvexHull for Geometry {
    fn convex_hull(&self) -> Polygon {
        let coords: Vec<Coord> = self.coords_iter().collect();
        Polygon::new(LineString::new(convex_hull_coords(&coords)), vec![])
    }
}

impl Centroid for Geometry {
    fn centroid(&self) -> Option<Point> {
        match self {
            Geometry::Point(g) => g.centroid(),
            Geometry::LineString(g) => g.centroid(),
            Geometry::Polygon(g) => g.centroid(),
            Geometry::MultiPolygon(g) => g.centroid(),
            Geometry::MultiPoint(g) => g.centroid(),
            Geometry::Rect(g) => Some(g.center()),
            Geometry::Triangle(g) => g.to_polygon().centroid(),
            Geometry::Line(l) => Some(Point((l.start + l.end) / 2.0)),
            Geometry::MultiLineString(mls) => {
                let pts: Vec<Point> = mls.iter().filter_map(Centroid::centroid).collect();
                MultiPoint::new(pts).centroid()
            }
            Geometry::GeometryCollection(gs) => {
                let pts: Vec<Point> = gs.iter().filter_map(Centroid::centroid).collect();
                MultiPoint::new(pts).centroid()
            }
        }
    }
}

impl Simplify for Geometry {
    fn simplify(&self, epsilon: f64) -> Self {
        match self {
            Geometry::LineString(g) => Geometry::LineString(g.simplify(epsilon)),
            Geometry::Polygon(g) => Geometry::Polygon(g.simplify(epsilon)),
            Geometry::MultiLineString(g) => Geometry::MultiLineString(g.simplify(epsilon)),
            Geometry::MultiPolygon(g) => Geometry::MultiPolygon(g.simplify(epsilon)),
            // nothing to simplify on points / lines / rects / triangles / collections
            other => other.clone(),
        }
    }
}
