//! A representative point guaranteed to intersect the geometry, mirroring
//! `geo`'s `InteriorPoint` trait. Areal types use the midpoint of the widest
//! intersection of a horizontal scanline through the middle of the bounding
//! box (perturbed off vertices); lines use the non-endpoint vertex closest to
//! the centroid; points return themselves. Multi-geometries pick the best
//! member candidate (widest scanline segment for polygons, closest to the
//! collection centroid otherwise). Cross-checked against `geo` in
//! `tests/interior_point.rs`.

use std::cmp::{Ordering, Reverse};

use super::dimensions::{Dimensions, HasDimensions};
use super::relate::relate;
use super::{
    BoundingRect, Centroid, Coord, CoordsIter, Geometry, Line, LineString, MultiLineString,
    MultiPoint, MultiPolygon, Point, Polygon, Rect, Triangle,
};

/// A point guaranteed to be on the geometry: strictly interior for areal types
/// when possible, on the geometry otherwise, and reasonably central.
pub trait InteriorPoint {
    type Output;

    /// A representative point inside (or, for zero-area inputs, on) the geometry.
    fn interior_point(&self) -> Self::Output;
}

fn distance(a: Point, b: Point) -> f64 {
    (a.x() - b.x()).hypot(a.y() - b.y())
}

impl InteriorPoint for Point {
    type Output = Point;

    fn interior_point(&self) -> Self::Output {
        *self
    }
}

impl InteriorPoint for Line {
    type Output = Point;

    /// The start point: a computed midpoint is not guaranteed to lie exactly on
    /// the segment under floating-point rounding, the endpoint is.
    fn interior_point(&self) -> Self::Output {
        Point(self.start)
    }
}

impl InteriorPoint for LineString {
    type Output = Option<Point>;

    /// The non-endpoint vertex closest to the centroid, or the start point when
    /// there are no non-endpoint vertices.
    fn interior_point(&self) -> Self::Output {
        match self.0.len() {
            0 => None,
            1 | 2 => Some(Point(self.0[0])),
            _ => {
                let centroid = self.centroid()?;
                self.0[1..self.0.len() - 1]
                    .iter()
                    .map(|&c| {
                        let pt = Point(c);
                        (pt, distance(pt, centroid))
                    })
                    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less))
                    .map(|(pt, _)| pt)
            }
        }
    }
}

impl InteriorPoint for MultiLineString {
    type Output = Option<Point>;

    /// Of the members' interior points, the one closest to the centroid of the
    /// whole collection.
    fn interior_point(&self) -> Self::Output {
        let centroid = Geometry::MultiLineString(self.clone()).centroid()?;
        self.iter()
            .filter_map(|ls| ls.interior_point().map(|pt| (pt, distance(pt, centroid))))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less))
            .map(|(pt, _)| pt)
    }
}

/// X-coordinates where the horizontal line `y = y_scan` meets segment `a`–`b`:
/// both endpoints for a collinear horizontal segment, the touching endpoint or
/// the crossing point otherwise.
fn scanline_hits(a: Coord, b: Coord, y_scan: f64, xs: &mut Vec<f64>) {
    let (da, db) = (a.y - y_scan, b.y - y_scan);
    if da == 0.0 && db == 0.0 {
        xs.push(a.x);
        xs.push(b.x);
    } else if da == 0.0 {
        xs.push(a.x);
    } else if db == 0.0 {
        xs.push(b.x);
    } else if (da > 0.0) != (db > 0.0) {
        let t = da / (a.y - b.y);
        xs.push(a.x + t * (b.x - a.x));
    }
}

/// The widest-scanline interior point of `polygon` plus the width of its
/// scanline segment (`0` when the point is only on the boundary), so a
/// multi-polygon can rank candidates. `None` only for an empty polygon.
fn polygon_interior_point_with_segment_length(polygon: &Polygon) -> Option<(Point, f64)> {
    if polygon.exterior().0.len() == 1 {
        return Some((Point(polygon.exterior().0[0]), 0.0));
    }

    let bounds = polygon.bounding_rect()?;

    let mut y_mid = (bounds.min().y + bounds.max().y) / 2.0;
    if polygon.coords_iter().any(|c| c.y == y_mid) {
        let next_closest = polygon
            .coords_iter()
            .filter_map(|c| (c.y != y_mid).then(|| (c.y, (c.y - y_mid).abs())))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less));
        if let Some((closest, _)) = next_closest {
            y_mid = (y_mid + closest) / 2.0;
        }
    }

    let mut xs: Vec<f64> = Vec::new();
    let rings = std::iter::once(polygon.exterior()).chain(polygon.interiors().iter());
    for ring in rings {
        for w in ring.0.windows(2) {
            scanline_hits(w[0], w[1], y_mid, &mut xs);
        }
    }
    xs.sort_by(f64::total_cmp);

    let mut segments: Vec<(Point, f64)> = xs
        .windows(2)
        .map(|w| (Point::new((w[0] + w[1]) / 2.0, y_mid), w[1] - w[0]))
        .collect();
    segments.sort_by(|a, b| b.1.total_cmp(&a.1));

    let polygon_geom = Geometry::Polygon(polygon.clone());
    for (midpoint, segment_length) in segments {
        let relation = relate(&polygon_geom, &Geometry::Point(midpoint));
        if relation.is_intersects() {
            let length = if relation.is_contains() {
                segment_length
            } else {
                0.0
            };
            return Some((midpoint, length));
        }
    }

    polygon.coords_iter().next().map(|c| (Point(c), 0.0))
}

impl InteriorPoint for Polygon {
    type Output = Option<Point>;

    fn interior_point(&self) -> Self::Output {
        polygon_interior_point_with_segment_length(self).map(|(point, _)| point)
    }
}

impl InteriorPoint for MultiPolygon {
    type Output = Option<Point>;

    /// The member candidate whose scanline segment is widest.
    fn interior_point(&self) -> Self::Output {
        self.iter()
            .filter_map(polygon_interior_point_with_segment_length)
            .min_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Less))
            .map(|(point, _)| point)
    }
}

impl InteriorPoint for Rect {
    type Output = Point;

    fn interior_point(&self) -> Self::Output {
        self.center()
    }
}

impl InteriorPoint for Triangle {
    type Output = Point;

    fn interior_point(&self) -> Self::Output {
        Point((self.0 + self.1 + self.2) / 3.0)
    }
}

impl InteriorPoint for MultiPoint {
    type Output = Option<Point>;

    /// The member point closest to the centroid.
    fn interior_point(&self) -> Self::Output {
        let centroid = self.centroid()?;
        self.iter()
            .map(|&pt| (pt, distance(pt, centroid)))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less))
            .map(|(pt, _)| pt)
    }
}

fn dimensions_rank(d: Dimensions) -> u8 {
    match d {
        Dimensions::Empty => 0,
        Dimensions::ZeroDimensional => 1,
        Dimensions::OneDimensional => 2,
        Dimensions::TwoDimensional => 3,
    }
}

impl InteriorPoint for Geometry {
    type Output = Option<Point>;

    fn interior_point(&self) -> Self::Output {
        match self {
            Geometry::Point(g) => Some(g.interior_point()),
            Geometry::Line(g) => Some(g.interior_point()),
            Geometry::LineString(g) => g.interior_point(),
            Geometry::Polygon(g) => g.interior_point(),
            Geometry::MultiPoint(g) => g.interior_point(),
            Geometry::MultiLineString(g) => g.interior_point(),
            Geometry::MultiPolygon(g) => g.interior_point(),
            Geometry::Rect(g) => Some(g.interior_point()),
            Geometry::Triangle(g) => Some(g.interior_point()),
            Geometry::GeometryCollection(gs) => {
                let centroid = self.centroid()?;
                gs.iter()
                    .filter_map(|g| {
                        g.interior_point().map(|pt| {
                            (
                                pt,
                                (
                                    Reverse(dimensions_rank(g.dimensions())),
                                    distance(pt, centroid),
                                ),
                            )
                        })
                    })
                    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less))
                    .map(|(pt, _)| pt)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::{Contains, Intersects};
    use super::*;

    fn poly(rings: &[&[(f64, f64)]]) -> Polygon {
        Polygon::new(
            LineString::from(rings[0].to_vec()),
            rings[1..]
                .iter()
                .map(|r| LineString::from(r.to_vec()))
                .collect(),
        )
    }

    #[test]
    fn rhombus_interior_point_matches_geo_doc_example() {
        let polygon = poly(&[&[(-2.0, 1.0), (1.0, 3.0), (4.0, 1.0), (1.0, -1.0)]]);
        assert_eq!(polygon.interior_point(), Some(Point::new(1.0, 2.0)));
    }

    #[test]
    fn polygon_with_hole_avoids_the_hole() {
        let polygon = poly(&[
            &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
            &[(2.0, 2.0), (8.0, 2.0), (8.0, 8.0), (2.0, 8.0)],
        ]);
        let pt = polygon.interior_point().expect("non-empty polygon");
        assert!(polygon.contains(&pt));
    }

    #[test]
    fn line_string_uses_a_vertex() {
        let line = LineString::from(vec![(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (5.0, 5.0)]);
        let pt = line.interior_point().expect("non-empty line");
        assert!(line.0.contains(&pt.0));
        assert_ne!(pt.0, line.0[0]);
        assert_ne!(pt.0, line.0[line.0.len() - 1]);
    }

    #[test]
    fn empty_inputs_yield_none() {
        assert_eq!(LineString::new(vec![]).interior_point(), None);
        assert_eq!(MultiPoint::new(vec![]).interior_point(), None);
        assert_eq!(MultiPolygon::new(vec![]).interior_point(), None);
    }

    #[test]
    fn interior_point_intersects_its_geometry() {
        let shapes = vec![
            Geometry::Polygon(poly(&[&[(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (3.0, 3.0)]])),
            Geometry::MultiPolygon(MultiPolygon::new(vec![
                poly(&[&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]]),
                poly(&[&[(5.0, 5.0), (9.0, 5.0), (9.0, 9.0), (5.0, 9.0)]]),
            ])),
            Geometry::LineString(LineString::from(vec![(0.0, 0.0), (2.0, 3.0), (4.0, 0.0)])),
            Geometry::MultiPoint(MultiPoint::from(vec![(1.0, 1.0), (3.0, 3.0)])),
        ];
        for g in shapes {
            let pt = g.interior_point().expect("non-empty geometry");
            assert!(g.intersects(&Geometry::Point(pt)), "{g:?} missed {pt:?}");
        }
    }
}
