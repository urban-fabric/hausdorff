//! [`Densify`] — a copy of a geometry with vertices inserted so no segment is
//! longer than a cap, mirroring `geo`'s `Densify`. Available on any metric
//! space implementing [`Distance`] + [`InterpolatePoint`] (today [`Euclidean`]
//! and [`Rhumb`](super::Rhumb)); segment lengths and the inserted points are
//! both measured in that metric.
//!
//! [`Euclidean`]: super::Euclidean

use super::line_measures::{Distance, InterpolatePoint};
use super::{Line, LineString, MultiLineString, MultiPolygon, Point, Polygon, Rect, Triangle};

/// A copy of the geometry with points inserted as necessary so that no two
/// successive vertices are more than `max_segment_length` apart under the
/// metric (`geo`'s `Densify`). `max_segment_length` must be greater than `0`.
///
/// `Output` differs from the input only where densifying changes the natural
/// shape: a [`Line`] becomes a [`LineString`], a [`Rect`]/[`Triangle`] a
/// [`Polygon`].
pub trait Densify<G> {
    /// The densified shape (the input type, widened where necessary).
    type Output;

    fn densify(&self, geometry: &G, max_segment_length: f64) -> Self::Output;
}

/// Pushes the interior points subdividing `start`–`end` into runs of at most
/// `max_segment_length` onto `points` (endpoints excluded). Subdivision is at
/// equal ratios of the segment's metric length, matching `geo`.
fn densify_between<M>(
    metric: &M,
    start: Point,
    end: Point,
    points: &mut Vec<Point>,
    max_segment_length: f64,
) where
    M: Distance<Point, Point> + InterpolatePoint,
{
    assert!(
        max_segment_length > 0.0,
        "max_segment_length must be greater than 0"
    );
    let num_segments = (metric.distance(start, end) / max_segment_length).ceil() as u64;
    let frac = 1.0 / num_segments as f64;
    for segment_num in 1..num_segments {
        let ratio = frac * segment_num as f64;
        points.push(metric.point_at_ratio_between(start, end, ratio));
    }
}

impl<M> Densify<Line> for M
where
    M: Distance<Point, Point> + InterpolatePoint,
{
    type Output = LineString;

    fn densify(&self, line: &Line, max_segment_length: f64) -> LineString {
        let (start, end) = (Point(line.start), Point(line.end));
        let mut points = vec![start];
        densify_between(self, start, end, &mut points, max_segment_length);
        points.push(end);
        LineString::from(points)
    }
}

impl<M> Densify<LineString> for M
where
    M: Distance<Point, Point> + InterpolatePoint,
{
    type Output = LineString;

    fn densify(&self, ls: &LineString, max_segment_length: f64) -> LineString {
        if ls.is_empty() {
            return LineString::new(vec![]);
        }

        let mut points = vec![];
        for line in ls.lines() {
            points.push(Point(line.start));
            densify_between(
                self,
                Point(line.start),
                Point(line.end),
                &mut points,
                max_segment_length,
            );
        }
        points.push(Point(*ls.0.last().expect("non-empty by the guard above")));
        LineString::from(points)
    }
}

impl<M> Densify<MultiLineString> for M
where
    M: Distance<Point, Point> + InterpolatePoint,
{
    type Output = MultiLineString;

    fn densify(&self, mls: &MultiLineString, max_segment_length: f64) -> MultiLineString {
        MultiLineString::new(
            mls.iter()
                .map(|ls| self.densify(ls, max_segment_length))
                .collect(),
        )
    }
}

impl<M> Densify<Polygon> for M
where
    M: Distance<Point, Point> + InterpolatePoint,
{
    type Output = Polygon;

    fn densify(&self, polygon: &Polygon, max_segment_length: f64) -> Polygon {
        Polygon::new(
            self.densify(polygon.exterior(), max_segment_length),
            polygon
                .interiors()
                .iter()
                .map(|hole| self.densify(hole, max_segment_length))
                .collect(),
        )
    }
}

impl<M> Densify<MultiPolygon> for M
where
    M: Distance<Point, Point> + InterpolatePoint,
{
    type Output = MultiPolygon;

    fn densify(&self, mp: &MultiPolygon, max_segment_length: f64) -> MultiPolygon {
        MultiPolygon::new(
            mp.iter()
                .map(|polygon| self.densify(polygon, max_segment_length))
                .collect(),
        )
    }
}

impl<M> Densify<Rect> for M
where
    M: Distance<Point, Point> + InterpolatePoint,
{
    type Output = Polygon;

    fn densify(&self, rect: &Rect, max_segment_length: f64) -> Polygon {
        self.densify(&rect.to_polygon(), max_segment_length)
    }
}

impl<M> Densify<Triangle> for M
where
    M: Distance<Point, Point> + InterpolatePoint,
{
    type Output = Polygon;

    fn densify(&self, triangle: &Triangle, max_segment_length: f64) -> Polygon {
        self.densify(&triangle.to_polygon(), max_segment_length)
    }
}
