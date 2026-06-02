//! Trait dispatch over the [`Geometry`] enum, so the planar algorithms run on a
//! heterogeneous geometry (not only the concrete types) — the native equivalent
//! of calling `geo`'s traits on a `geo::Geometry`.

use super::convex_hull::convex_hull_coords;
use super::{
    Area, BoundingRect, Centroid, ChamberlainDuquetteArea, ConvexHull, Coord, CoordsIter, Euclidean,
    Geometry, Length, LineString, MultiPoint, Point, Polygon, Rect, Simplify,
};

impl ChamberlainDuquetteArea for Geometry {
    fn chamberlain_duquette_signed_area(&self) -> f64 {
        match self {
            Geometry::Polygon(g) => g.chamberlain_duquette_signed_area(),
            Geometry::MultiPolygon(g) => g.chamberlain_duquette_signed_area(),
            Geometry::Rect(g) => g.to_polygon().chamberlain_duquette_signed_area(),
            Geometry::Triangle(g) => g.to_polygon().chamberlain_duquette_signed_area(),
            Geometry::GeometryCollection(gs) => {
                gs.iter().map(ChamberlainDuquetteArea::chamberlain_duquette_signed_area).sum()
            }
            _ => 0.0,
        }
    }
}

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

    fn unsigned_area(&self) -> f64 {
        // A collection's unsigned area is the sum of its parts' *unsigned* areas,
        // not `|Σ signed|` — members of opposite winding must not cancel.
        match self {
            Geometry::GeometryCollection(gs) => gs.iter().map(Area::unsigned_area).sum(),
            other => other.signed_area().abs(),
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
                // Length-weighted mean of the component centroids (a long line
                // pulls the centroid more than a short one).
                let mut acc = Coord::new(0.0, 0.0);
                let mut total = 0.0;
                for ls in mls.iter() {
                    if let Some(c) = ls.centroid() {
                        let len = Euclidean.length(ls);
                        acc = acc + c.0 * len;
                        total += len;
                    }
                }
                if total > 0.0 {
                    Some(Point(acc / total))
                } else {
                    // every component is zero-length: fall back to their mean
                    let pts: Vec<Point> = mls.iter().filter_map(Centroid::centroid).collect();
                    MultiPoint::new(pts).centroid()
                }
            }
            Geometry::GeometryCollection(_) => {
                // The mean of every coordinate in the collection.
                let mut sum = Coord::new(0.0, 0.0);
                let mut n = 0.0;
                for c in self.coords_iter() {
                    sum = sum + c;
                    n += 1.0;
                }
                (n > 0.0).then(|| Point(sum / n))
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
            // recurse into a collection, simplifying each member
            Geometry::GeometryCollection(gs) => {
                Geometry::GeometryCollection(gs.iter().map(|g| g.simplify(epsilon)).collect())
            }
            // nothing to simplify on points / lines / rects / triangles
            other => other.clone(),
        }
    }
}
