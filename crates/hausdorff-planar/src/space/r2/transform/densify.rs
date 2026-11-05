//! [`Densify`] — insert vertices along each segment of a geometry so no
//! segment exceeds a maximum length *under a chosen metric*, mirroring `geo`'s
//! `Densify`. The inserted vertices lie **on the metric's own path**:
//! straight lines for [`Euclidean`], great circles for [`Haversine`], true
//! WGS84 geodesics for [`Geodesic`] (via [`hausdorff_base::earth::GeodesicLine`], so a
//! densified segment is sampled from one precomputed line, not n separate
//! solves), and constant-bearing loxodromes for [`Rhumb`].
//!
//! Densifying before `Projection::project_geometry` (in `hausdorff-earth`)
//! is how a long edge gets carried into a projection faithfully: the projected
//! polyline then follows the projected curve instead of cutting straight
//! between distant endpoints.

use hausdorff_base::earth::{
    Bearing as CompassBearing, Ellipsoid, GeoPoint, destination_m, haversine_m, rhumb_distance_m,
    rhumb_point_at_ratio,
};

use crate::space::r2::{
    Coord, Geometry, Line, LineString, MultiLineString, MultiPolygon, Polygon, Rect, Triangle,
};
use crate::space::r2::{Distance, Euclidean, Geodesic, Haversine, Rhumb};

/// Densify a geometry under this metric so no segment is longer than
/// `max_segment_length` (in the metric's unit — metres for [`Haversine`] /
/// [`Geodesic`] / [`Rhumb`], the plane's unit for [`Euclidean`]).
///
/// `Output` differs from the input only where densifying changes the natural
/// shape: a [`Line`] becomes a [`LineString`], a [`Rect`] / [`Triangle`] a
/// [`Polygon`] (the fixed-vertex types cannot hold inserted vertices).
///
/// Proof-backed by `riemann` `Hausdorff/Densify`: inserted samples lie on the
/// segment as the convex combination `(1−i/n)a + (i/n)b` (`densPoint_eq_combo`),
/// are evenly spaced (`densPoint_step`), and preserve length exactly — the pieces
/// sum to `‖b−a‖` (`densify_length`); the `⌈‖b−a‖/maxLen⌉` count keeps every piece
/// `≤ maxLen` (`pieces_satisfies_cap`).
pub trait Densify<G> {
    /// The densified shape (the input type, widened where necessary).
    type Output;

    /// `geometry` with extra vertices inserted along each over-long segment.
    /// Original vertices are always preserved.
    ///
    /// # Panics
    /// Panics if `max_segment_length` is not finite and positive.
    #[must_use]
    fn densify(&self, geometry: &G, max_segment_length: f64) -> Self::Output;
}

/// The per-metric kernel: append `a` plus the interior vertices of the `a → b`
/// segment (not `b` itself) to `out`, splitting so no piece exceeds `max_len`.
trait SegmentFill {
    fn fill_segment(&self, a: Coord, b: Coord, max_len: f64, out: &mut Vec<Coord>);
}

/// Number of pieces the segment must split into (≥ 1; degenerate or NaN
/// distances collapse to a single piece).
fn pieces(length: f64, max_len: f64) -> usize {
    if length.is_finite() && length > max_len {
        (length / max_len).ceil() as usize
    } else {
        1
    }
}

impl SegmentFill for Euclidean {
    fn fill_segment(&self, a: Coord, b: Coord, max_len: f64, out: &mut Vec<Coord>) {
        out.push(a);
        let n = pieces(self.distance(a, b), max_len);
        for i in 1..n {
            let t = i as f64 / n as f64;
            out.push(a + (b - a) * t);
        }
    }
}

impl SegmentFill for Haversine {
    fn fill_segment(&self, a: Coord, b: Coord, max_len: f64, out: &mut Vec<Coord>) {
        out.push(a);
        let ga = GeoPoint::from_lon_lat(a.x, a.y);
        let gb = GeoPoint::from_lon_lat(b.x, b.y);
        let length = haversine_m(ga, gb);
        let n = pieces(length, max_len);
        if n > 1 {
            // Walking the initial bearing traces the great circle through
            // a and b, so every sample stays on the spherical path.
            let bearing = CompassBearing::initial(ga, gb);
            for i in 1..n {
                let p = destination_m(ga, bearing, length * i as f64 / n as f64);
                out.push(Coord::new(p.lon, p.lat));
            }
        }
    }
}

impl SegmentFill for Geodesic {
    fn fill_segment(&self, a: Coord, b: Coord, max_len: f64, out: &mut Vec<Coord>) {
        out.push(a);
        let ga = GeoPoint::from_lon_lat(a.x, a.y);
        let gb = GeoPoint::from_lon_lat(b.x, b.y);
        let inv = Ellipsoid::WGS84.inverse(ga, gb);
        let n = pieces(inv.distance_m, max_len);
        if n > 1 {
            let line = Ellipsoid::WGS84.line(ga, inv.initial_bearing);
            for i in 1..n {
                let p = line.point_m(inv.distance_m * i as f64 / n as f64);
                out.push(Coord::new(p.lon, p.lat));
            }
        }
    }
}

impl SegmentFill for Rhumb {
    fn fill_segment(&self, a: Coord, b: Coord, max_len: f64, out: &mut Vec<Coord>) {
        out.push(a);
        let ga = GeoPoint::from_lon_lat(a.x, a.y);
        let gb = GeoPoint::from_lon_lat(b.x, b.y);
        let n = pieces(rhumb_distance_m(ga, gb), max_len);
        for i in 1..n {
            let p = rhumb_point_at_ratio(ga, gb, i as f64 / n as f64);
            out.push(Coord::new(p.lon, p.lat));
        }
    }
}

fn densify_line_string<M: SegmentFill>(m: &M, ls: &LineString, max_len: f64) -> LineString {
    let coords = &ls.0;
    let Some(&last) = coords.last() else {
        return LineString::new(Vec::new());
    };
    let mut out = Vec::with_capacity(coords.len());
    for w in coords.windows(2) {
        m.fill_segment(w[0], w[1], max_len, &mut out);
    }
    out.push(last);
    LineString::new(out)
}

fn assert_max_len(max_segment_length: f64) {
    assert!(
        max_segment_length.is_finite() && max_segment_length > 0.0,
        "densify: max_segment_length must be finite and positive, got {max_segment_length}"
    );
}

impl<M: SegmentFill> Densify<LineString> for M {
    type Output = LineString;

    fn densify(&self, ls: &LineString, max_segment_length: f64) -> LineString {
        assert_max_len(max_segment_length);
        densify_line_string(self, ls, max_segment_length)
    }
}

impl<M: SegmentFill> Densify<Line> for M {
    type Output = LineString;

    fn densify(&self, line: &Line, max_segment_length: f64) -> LineString {
        assert_max_len(max_segment_length);
        densify_line_string(
            self,
            &LineString::new(vec![line.start, line.end]),
            max_segment_length,
        )
    }
}

impl<M: SegmentFill> Densify<MultiLineString> for M {
    type Output = MultiLineString;

    fn densify(&self, mls: &MultiLineString, max_segment_length: f64) -> MultiLineString {
        assert_max_len(max_segment_length);
        MultiLineString(
            mls.0
                .iter()
                .map(|ls| densify_line_string(self, ls, max_segment_length))
                .collect(),
        )
    }
}

impl<M: SegmentFill> Densify<Polygon> for M {
    type Output = Polygon;

    fn densify(&self, poly: &Polygon, max_segment_length: f64) -> Polygon {
        assert_max_len(max_segment_length);
        Polygon::new(
            densify_line_string(self, poly.exterior(), max_segment_length),
            poly.interiors()
                .iter()
                .map(|r| densify_line_string(self, r, max_segment_length))
                .collect(),
        )
    }
}

impl<M: SegmentFill> Densify<MultiPolygon> for M {
    type Output = MultiPolygon;

    fn densify(&self, mp: &MultiPolygon, max_segment_length: f64) -> MultiPolygon {
        assert_max_len(max_segment_length);
        MultiPolygon(
            mp.0.iter()
                .map(|p| self.densify(p, max_segment_length))
                .collect(),
        )
    }
}

impl<M: SegmentFill> Densify<Rect> for M {
    type Output = Polygon;

    fn densify(&self, rect: &Rect, max_segment_length: f64) -> Polygon {
        self.densify(&rect.to_polygon(), max_segment_length)
    }
}

impl<M: SegmentFill> Densify<Triangle> for M {
    type Output = Polygon;

    fn densify(&self, triangle: &Triangle, max_segment_length: f64) -> Polygon {
        self.densify(&triangle.to_polygon(), max_segment_length)
    }
}

/// Linear and areal parts densify; puntal parts pass through. A [`Line`],
/// [`Rect`], or [`Triangle`] variant becomes the densified [`LineString`] /
/// [`Polygon`] it describes (the fixed-vertex types cannot hold inserted
/// vertices).
impl<M: SegmentFill> Densify<Geometry> for M {
    type Output = Geometry;

    fn densify(&self, g: &Geometry, max_segment_length: f64) -> Geometry {
        assert_max_len(max_segment_length);
        match g {
            Geometry::Point(_) | Geometry::MultiPoint(_) => g.clone(),
            Geometry::Line(l) => Geometry::LineString(densify_line_string(
                self,
                &LineString::new(vec![l.start, l.end]),
                max_segment_length,
            )),
            Geometry::LineString(ls) => {
                Geometry::LineString(densify_line_string(self, ls, max_segment_length))
            }
            Geometry::MultiLineString(mls) => {
                Geometry::MultiLineString(self.densify(mls, max_segment_length))
            }
            Geometry::Polygon(p) => Geometry::Polygon(self.densify(p, max_segment_length)),
            Geometry::MultiPolygon(mp) => {
                Geometry::MultiPolygon(self.densify(mp, max_segment_length))
            }
            Geometry::Rect(r) => {
                Geometry::Polygon(self.densify(&r.to_polygon(), max_segment_length))
            }
            Geometry::Triangle(t) => {
                Geometry::Polygon(self.densify(&t.to_polygon(), max_segment_length))
            }
            Geometry::GeometryCollection(gs) => Geometry::GeometryCollection(
                gs.iter()
                    .map(|sub| self.densify(sub, max_segment_length))
                    .collect(),
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::r2::Length;

    #[test]
    fn euclidean_densify_splits_evenly() {
        let ls = LineString::from(vec![(0.0, 0.0), (10.0, 0.0)]);
        let dense = Euclidean.densify(&ls, 3.0);
        // 10/3 → 4 pieces → 5 vertices, evenly spaced.
        assert_eq!(dense.0.len(), 5);
        assert_eq!(dense.0[2], Coord::new(5.0, 0.0));
        assert_eq!(Euclidean.length(&dense), 10.0);
    }

    #[test]
    fn metric_densify_preserves_length_and_caps_segments() {
        // Paris → Moscow, one ~2500 km segment. The inserted vertices lie on
        // the metric's own path, so the total length is unchanged and every
        // piece respects the cap.
        let ls = LineString::from(vec![(2.35, 48.85), (37.62, 55.75)]);

        let hav = Haversine.densify(&ls, 100_000.0);
        assert!(hav.0.len() > ls.0.len(), "no vertices inserted");
        assert!(hav.0.first() == ls.0.first() && hav.0.last() == ls.0.last());
        assert!((Haversine.length(&hav) - Haversine.length(&ls)).abs() < 1e-6);
        for w in hav.0.windows(2) {
            assert!(Haversine.distance(w[0], w[1]) <= 100_000.0 + 1e-6);
        }

        let geo = Geodesic.densify(&ls, 100_000.0);
        assert!(geo.0.len() > ls.0.len(), "no vertices inserted");
        assert!(geo.0.first() == ls.0.first() && geo.0.last() == ls.0.last());
        assert!((Geodesic.length(&geo) - Geodesic.length(&ls)).abs() < 1e-6);
        for w in geo.0.windows(2) {
            assert!(Geodesic.distance(w[0], w[1]) <= 100_000.0 + 1e-6);
        }

        let rhumb = Rhumb.densify(&ls, 100_000.0);
        assert!(rhumb.0.len() > ls.0.len(), "no vertices inserted");
        assert!(rhumb.0.first() == ls.0.first() && rhumb.0.last() == ls.0.last());
        assert!((Rhumb.length(&rhumb) - Rhumb.length(&ls)).abs() < 1e-6);
        for w in rhumb.0.windows(2) {
            assert!(Rhumb.distance(w[0], w[1]) <= 100_000.0 + 1e-6);
        }
    }

    #[test]
    fn line_densifies_into_a_line_string() {
        let line = Line::new(Coord::new(0.0, 0.0), Coord::new(10.0, 0.0));
        let dense: LineString = Euclidean.densify(&line, 3.0);
        assert_eq!(dense.0.len(), 5);
        assert_eq!(dense.0.first(), Some(&line.start));
        assert_eq!(dense.0.last(), Some(&line.end));
    }

    #[test]
    fn polygon_densify_keeps_rings_closed() {
        let poly = Polygon::new(
            LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
            vec![],
        );
        let dense: Polygon = Haversine.densify(&poly, 20_000.0);
        assert!(dense.exterior().is_closed());
        assert!(dense.exterior().0.len() > poly.exterior().0.len());
    }

    #[test]
    #[should_panic(expected = "max_segment_length")]
    fn zero_max_length_panics() {
        let ls = LineString::from(vec![(0.0, 0.0), (1.0, 0.0)]);
        let _ = Euclidean.densify(&ls, 0.0);
    }
}
