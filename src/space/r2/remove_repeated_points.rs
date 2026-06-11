//! Repeated-point removal, mirroring `geo`'s `RemoveRepeatedPoints` trait:
//! drop *consecutive* duplicate coordinates from line strings and polygon
//! rings, and *all* duplicate points (keeping first occurrences) from a
//! [`MultiPoint`]. Validated against `geo` in `tests/remove_repeated_points.rs`.

use super::{
    Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect,
    Triangle,
};

/// Remove repeated points from a [`MultiPoint`] and repeated consecutive
/// coordinates from [`LineString`], [`Polygon`], [`MultiLineString`] and
/// [`MultiPolygon`].
///
/// For [`Geometry::GeometryCollection`] each member is cleaned individually.
/// [`Point`], [`Line`], [`Rect`] and [`Triangle`] are returned unchanged.
pub trait RemoveRepeatedPoints {
    /// A copy of the geometry with (consecutive) repeated points removed.
    fn remove_repeated_points(&self) -> Self;
    /// Remove (consecutive) repeated points in place.
    fn remove_repeated_points_mut(&mut self);
}

impl RemoveRepeatedPoints for MultiPoint {
    fn remove_repeated_points(&self) -> Self {
        let mut points: Vec<Point> = Vec::new();
        for p in self.iter() {
            if !points.contains(p) {
                points.push(*p);
            }
        }
        MultiPoint(points)
    }

    fn remove_repeated_points_mut(&mut self) {
        *self = self.remove_repeated_points();
    }
}

impl RemoveRepeatedPoints for LineString {
    fn remove_repeated_points(&self) -> Self {
        let mut coords = self.0.clone();
        coords.dedup();
        LineString(coords)
    }

    fn remove_repeated_points_mut(&mut self) {
        self.0.dedup();
    }
}

impl RemoveRepeatedPoints for Polygon {
    fn remove_repeated_points(&self) -> Self {
        Polygon::new(
            self.exterior().remove_repeated_points(),
            self.interiors()
                .iter()
                .map(RemoveRepeatedPoints::remove_repeated_points)
                .collect(),
        )
    }

    fn remove_repeated_points_mut(&mut self) {
        let (mut exterior, mut interiors) = std::mem::take(self).into_inner();
        exterior.remove_repeated_points_mut();
        for ring in &mut interiors {
            ring.remove_repeated_points_mut();
        }
        *self = Polygon::new(exterior, interiors);
    }
}

impl RemoveRepeatedPoints for MultiLineString {
    fn remove_repeated_points(&self) -> Self {
        MultiLineString::new(
            self.iter()
                .map(RemoveRepeatedPoints::remove_repeated_points)
                .collect(),
        )
    }

    fn remove_repeated_points_mut(&mut self) {
        for ls in &mut self.0 {
            ls.remove_repeated_points_mut();
        }
    }
}

impl RemoveRepeatedPoints for MultiPolygon {
    fn remove_repeated_points(&self) -> Self {
        MultiPolygon::new(
            self.iter()
                .map(RemoveRepeatedPoints::remove_repeated_points)
                .collect(),
        )
    }

    fn remove_repeated_points_mut(&mut self) {
        for p in &mut self.0 {
            p.remove_repeated_points_mut();
        }
    }
}

macro_rules! noop_remove_repeated_points {
    ($ty:ident) => {
        /// Single fixed-vertex geometry: nothing to deduplicate, so
        /// `remove_repeated_points` is a clone and the `_mut` form a no-op.
        impl RemoveRepeatedPoints for $ty {
            fn remove_repeated_points(&self) -> Self {
                *self
            }

            fn remove_repeated_points_mut(&mut self) {}
        }
    };
}

noop_remove_repeated_points!(Point);
noop_remove_repeated_points!(Line);
noop_remove_repeated_points!(Rect);
noop_remove_repeated_points!(Triangle);

impl RemoveRepeatedPoints for Geometry {
    fn remove_repeated_points(&self) -> Self {
        match self {
            Geometry::Point(g) => Geometry::Point(g.remove_repeated_points()),
            Geometry::Line(g) => Geometry::Line(g.remove_repeated_points()),
            Geometry::LineString(g) => Geometry::LineString(g.remove_repeated_points()),
            Geometry::Polygon(g) => Geometry::Polygon(g.remove_repeated_points()),
            Geometry::MultiPoint(g) => Geometry::MultiPoint(g.remove_repeated_points()),
            Geometry::MultiLineString(g) => Geometry::MultiLineString(g.remove_repeated_points()),
            Geometry::MultiPolygon(g) => Geometry::MultiPolygon(g.remove_repeated_points()),
            Geometry::Rect(g) => Geometry::Rect(g.remove_repeated_points()),
            Geometry::Triangle(g) => Geometry::Triangle(g.remove_repeated_points()),
            Geometry::GeometryCollection(gs) => Geometry::GeometryCollection(
                gs.iter()
                    .map(RemoveRepeatedPoints::remove_repeated_points)
                    .collect(),
            ),
        }
    }

    fn remove_repeated_points_mut(&mut self) {
        match self {
            Geometry::Point(g) => g.remove_repeated_points_mut(),
            Geometry::Line(g) => g.remove_repeated_points_mut(),
            Geometry::LineString(g) => g.remove_repeated_points_mut(),
            Geometry::Polygon(g) => g.remove_repeated_points_mut(),
            Geometry::MultiPoint(g) => g.remove_repeated_points_mut(),
            Geometry::MultiLineString(g) => g.remove_repeated_points_mut(),
            Geometry::MultiPolygon(g) => g.remove_repeated_points_mut(),
            Geometry::Rect(g) => g.remove_repeated_points_mut(),
            Geometry::Triangle(g) => g.remove_repeated_points_mut(),
            Geometry::GeometryCollection(gs) => {
                for g in gs {
                    g.remove_repeated_points_mut();
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::Coord;
    use super::*;

    fn ls(coords: &[(f64, f64)]) -> LineString {
        LineString::from(coords.to_vec())
    }

    #[test]
    fn line_string_consecutive_dedup_only() {
        let input = ls(&[
            (0.0, 0.0),
            (1.0, 1.0),
            (1.0, 1.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (0.0, 0.0),
        ]);
        let expected = ls(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (0.0, 0.0)]);
        assert_eq!(input.remove_repeated_points(), expected);

        let mut in_place = input;
        in_place.remove_repeated_points_mut();
        assert_eq!(in_place, expected);
    }

    #[test]
    fn multi_point_full_dedup_keeps_first_occurrence_order() {
        let input = MultiPoint::from(vec![
            (0.0, 0.0),
            (1.0, 1.0),
            (1.0, 1.0),
            (2.0, 2.0),
            (0.0, 0.0),
        ]);
        let expected = MultiPoint::from(vec![(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)]);
        assert_eq!(input.remove_repeated_points(), expected);
    }

    #[test]
    fn polygon_rings_are_deduped_and_stay_closed() {
        let poly = Polygon::new(
            ls(&[
                (0.0, 0.0),
                (4.0, 0.0),
                (4.0, 0.0),
                (4.0, 4.0),
                (0.0, 4.0),
                (0.0, 4.0),
            ]),
            vec![ls(&[(1.0, 1.0), (2.0, 1.0), (2.0, 1.0), (2.0, 2.0)])],
        );
        let out = poly.remove_repeated_points();
        assert!(out.exterior().is_closed());
        assert_eq!(out.exterior().0.len(), 5);
        assert_eq!(out.interiors()[0].0.len(), 4);
    }

    #[test]
    fn geometry_collection_recurses() {
        let g = Geometry::GeometryCollection(vec![Geometry::LineString(ls(&[
            (0.0, 0.0),
            (0.0, 0.0),
            (1.0, 1.0),
        ]))]);
        let cleaned = g.remove_repeated_points();
        match cleaned {
            Geometry::GeometryCollection(gs) => match &gs[0] {
                Geometry::LineString(out) => {
                    assert_eq!(out.0, vec![Coord::new(0.0, 0.0), Coord::new(1.0, 1.0)]);
                }
                other => panic!("unexpected {other:?}"),
            },
            other => panic!("unexpected {other:?}"),
        }
    }
}
