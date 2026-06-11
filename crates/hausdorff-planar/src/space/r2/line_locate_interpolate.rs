//! Linear referencing along a line — mirrors of `geo`'s `LineInterpolatePoint`
//! (fraction of total length → point) and `LineLocatePoint` (point → fraction
//! of total length of the closest location).
//!
//! Fractions are measured along the line's Euclidean length in `[0, 1]`;
//! out-of-range fractions clamp to the endpoints and `NaN` inputs (or
//! non-finite coordinates) yield `None`.

use super::{Coord, Euclidean, Length, Line, LineString, Point};

/// The point that lies a given fraction along a line.
///
/// Fractions below `0` (including negative infinity) return the starting
/// point; fractions above `1` (including infinity) return the ending point.
/// `NaN` fractions and non-finite coordinates return `None`.
pub trait LineInterpolatePoint {
    /// `Option<Point>` — `None` for degenerate inputs.
    type Output;

    fn line_interpolate_point(&self, fraction: f64) -> Self::Output;
}

impl LineInterpolatePoint for Line {
    type Output = Option<Point>;

    fn line_interpolate_point(&self, fraction: f64) -> Option<Point> {
        if (0.0..=1.0).contains(&fraction) {
            let diff = self.end - self.start;
            let r = self.start + diff * fraction;
            return r.is_finite().then_some(Point(r));
        }
        if fraction < 0.0 {
            return self.line_interpolate_point(0.0);
        }
        if fraction > 1.0 {
            return self.line_interpolate_point(1.0);
        }
        None
    }
}

impl LineInterpolatePoint for LineString {
    type Output = Option<Point>;

    fn line_interpolate_point(&self, fraction: f64) -> Option<Point> {
        if (0.0..=1.0).contains(&fraction) {
            let total_length = Euclidean.length(self);
            let fractional_length = total_length * fraction;
            let mut cum_length = 0.0;
            for segment in self.lines() {
                let length = segment.length();
                if cum_length + length >= fractional_length {
                    let segment_fraction = (fractional_length - cum_length) / length;
                    return segment.line_interpolate_point(segment_fraction);
                }
                cum_length += length;
            }
            return None;
        }
        if fraction < 0.0 {
            return self.line_interpolate_point(0.0);
        }
        if fraction > 1.0 {
            return self.line_interpolate_point(1.0);
        }
        None
    }
}

/// The fraction of a line's total length locating the closest point on the
/// line to the given point.
///
/// A zero-length line locates everything at fraction `0`. Non-finite inputs
/// (the point or any line coordinate) return `None`.
pub trait LineLocatePoint<Rhs = Point> {
    /// `Option<f64>` — `None` for degenerate inputs.
    type Output;

    fn line_locate_point(&self, p: &Rhs) -> Self::Output;
}

impl LineLocatePoint for Line {
    type Output = Option<f64>;

    fn line_locate_point(&self, p: &Point) -> Option<f64> {
        let sp = p.0 - self.start;
        let v = self.delta();
        let v_sq = v.dot(v);
        if v_sq == 0.0 {
            return Some(0.0);
        }
        let l = v.dot(sp) / v_sq;
        l.is_finite().then(|| l.clamp(0.0, 1.0))
    }
}

/// Distance from `p` to the bounded segment `a`–`b` (the projection clamped to
/// the segment).
fn point_segment_distance(p: Coord, a: Coord, b: Coord) -> f64 {
    let ab = b - a;
    let len2 = ab.dot(ab);
    let proj = if len2 == 0.0 {
        a
    } else {
        a + ab * ((p - a).dot(ab) / len2).clamp(0.0, 1.0)
    };
    let d = p - proj;
    d.dot(d).sqrt()
}

impl LineLocatePoint for LineString {
    type Output = Option<f64>;

    fn line_locate_point(&self, p: &Point) -> Option<f64> {
        let total_length = Euclidean.length(self);
        if total_length == 0.0 {
            return Some(0.0);
        }
        let mut cum_length = 0.0;
        let mut closest_dist_to_point = f64::INFINITY;
        let mut fraction = 0.0;
        for segment in self.lines() {
            let segment_distance_to_point = point_segment_distance(p.0, segment.start, segment.end);
            let segment_length = segment.length();
            let segment_fraction = segment.line_locate_point(p)?;
            if segment_distance_to_point < closest_dist_to_point {
                closest_dist_to_point = segment_distance_to_point;
                fraction = (cum_length + segment_fraction * segment_length) / total_length;
            }
            cum_length += segment_length;
        }
        Some(fraction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn elbow() -> LineString {
        LineString::from(vec![(-1.0, 0.0), (0.0, 0.0), (0.0, 1.0)])
    }

    #[test]
    fn interpolate_line() {
        let line = Line::new(Coord::new(-1.0, 0.0), Coord::new(1.0, 0.0));
        assert_eq!(
            line.line_interpolate_point(0.0),
            Some(Point::new(-1.0, 0.0))
        );
        assert_eq!(line.line_interpolate_point(0.5), Some(Point::new(0.0, 0.0)));
        assert_eq!(line.line_interpolate_point(1.0), Some(Point::new(1.0, 0.0)));
        assert_eq!(
            line.line_interpolate_point(-2.0),
            Some(Point::new(-1.0, 0.0))
        );
        assert_eq!(line.line_interpolate_point(2.0), Some(Point::new(1.0, 0.0)));
        assert_eq!(line.line_interpolate_point(f64::NAN), None);
    }

    #[test]
    fn interpolate_line_string() {
        let ls = elbow();
        assert_eq!(ls.line_interpolate_point(-1.0), Some(Point::new(-1.0, 0.0)));
        assert_eq!(ls.line_interpolate_point(0.25), Some(Point::new(-0.5, 0.0)));
        assert_eq!(ls.line_interpolate_point(0.5), Some(Point::new(0.0, 0.0)));
        assert_eq!(ls.line_interpolate_point(0.75), Some(Point::new(0.0, 0.5)));
        assert_eq!(ls.line_interpolate_point(2.0), Some(Point::new(0.0, 1.0)));
        assert_eq!(ls.line_interpolate_point(f64::NAN), None);
        assert_eq!(LineString::new(vec![]).line_interpolate_point(0.5), None);
    }

    #[test]
    fn locate_line() {
        let line = Line::new(Coord::new(0.0, 0.0), Coord::new(2.0, 0.0));
        assert_eq!(line.line_locate_point(&Point::new(0.5, 1.0)), Some(0.25));
        assert_eq!(line.line_locate_point(&Point::new(-5.0, 0.0)), Some(0.0));
        assert_eq!(line.line_locate_point(&Point::new(5.0, 0.0)), Some(1.0));
        let degenerate = Line::new(Coord::new(1.0, 1.0), Coord::new(1.0, 1.0));
        assert_eq!(
            degenerate.line_locate_point(&Point::new(2.0, 2.0)),
            Some(0.0)
        );
        assert_eq!(line.line_locate_point(&Point::new(f64::NAN, 0.0)), None);
    }

    #[test]
    fn locate_line_string() {
        let ls = elbow();
        assert_eq!(ls.line_locate_point(&Point::new(-1.0, 0.0)), Some(0.0));
        assert_eq!(ls.line_locate_point(&Point::new(-0.5, 0.0)), Some(0.25));
        assert_eq!(ls.line_locate_point(&Point::new(0.0, 0.0)), Some(0.5));
        assert_eq!(ls.line_locate_point(&Point::new(0.0, 0.5)), Some(0.75));
        assert_eq!(ls.line_locate_point(&Point::new(0.0, 1.0)), Some(1.0));
        assert_eq!(
            LineString::new(vec![]).line_locate_point(&Point::new(0.0, 0.0)),
            Some(0.0)
        );
    }

    #[test]
    fn round_trip_locate_interpolate() {
        let ls = LineString::from(vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)]);
        for i in 0..=20 {
            let f = i as f64 / 20.0;
            let p = ls.line_interpolate_point(f).unwrap();
            let back = ls.line_locate_point(&p).unwrap();
            assert!((back - f).abs() < 1e-12, "round trip {f} -> {back}");
        }
    }
}
