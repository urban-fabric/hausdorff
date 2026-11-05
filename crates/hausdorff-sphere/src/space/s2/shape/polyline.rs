//! `s2::Polyline` — a chain of great-circle edges on the unit sphere, after
//! golang/geo's `s2.Polyline`. A polyline is one-dimensional: as a [`Region`]
//! it never *contains* a point or a cell, but its [`Region::intersects_cell`]
//! is exact enough for covering — a cell intersects the polyline iff a vertex
//! lies in the cell or an edge meets the cell boundary — so a
//! [`crate::space::s2::RegionCoverer`] produces boundary covers of line geometries with no
//! false negatives.

use hausdorff_base::earth::GeoPoint;
use hausdorff_base::space::s1::Angle;

use super::region::{Region, bounding_cap_of, cell_loop};
use crate::space::s2::{Cap, Rect, S2CellId, S2Point, arc};
use crate::space::s2::{Crossing, crossing_sign};

/// A polyline on the unit sphere: consecutive vertices joined by great-circle
/// edges.
#[derive(Clone, Debug, PartialEq)]
pub struct Polyline {
    vertices: Vec<S2Point>,
}

impl Polyline {
    /// A polyline over `vertices` (an empty or single-vertex polyline is
    /// allowed; the latter behaves as a degenerate point).
    #[must_use]
    pub fn new(vertices: Vec<S2Point>) -> Self {
        Self { vertices }
    }

    /// A polyline through the given WGS84 coordinates.
    #[must_use]
    pub fn from_geo_points(points: &[GeoPoint]) -> Self {
        Self::new(points.iter().map(|&p| S2Point::from_geo_point(p)).collect())
    }

    /// The vertices in order.
    #[must_use]
    pub fn vertices(&self) -> &[S2Point] {
        &self.vertices
    }

    /// Number of vertices.
    #[must_use]
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Total length — the sum of the great-circle edge lengths — as an
    /// [`Angle`]; multiply by a planet radius for metres.
    ///
    /// Proof-backed by `riemann` `Hausdorff/SpherePolyline`: the length is
    /// non-negative (`spherePolyLen_nonneg`), a two-vertex polyline measures the
    /// endpoint great-circle distance (`spherePolyLen_pair`), and the **chord
    /// bound** holds — the direct great-circle distance between the first and last
    /// vertices is at most the total length (`sphereDist_le_polyLen`, the
    /// spherical triangle inequality chained). Same fold as `Loop::perimeter`.
    #[must_use]
    pub fn length(&self) -> Angle {
        let total: f64 = self
            .edges()
            .map(|(a, b)| a.vector().angle(b.vector()).radians())
            .sum();
        Angle::from_radians(total)
    }

    /// The directed edges `(vᵢ, vᵢ₊₁)` (empty for fewer than two vertices).
    fn edges(&self) -> impl Iterator<Item = (S2Point, S2Point)> + '_ {
        self.vertices.windows(2).map(|w| (w[0], w[1]))
    }

    /// A conservative lat/lng rectangle containing every great-circle edge.
    #[must_use]
    pub fn rect_bound(&self) -> Rect {
        match self.vertices.as_slice() {
            [] => Rect::empty(),
            [p] => {
                let g = p.to_geo_point();
                Rect::from_degrees(g.lat, g.lon, g.lat, g.lon)
            }
            _ => self
                .edges()
                .map(|(a, b)| arc::rect_bound(a, b))
                .fold(Rect::empty(), |acc, r| acc.union(&r)),
        }
    }
}

impl Region for Polyline {
    /// Always false: a polyline has no interior.
    fn contains_point(&self, _p: S2Point) -> bool {
        false
    }

    fn bounding_cap(&self) -> Cap {
        bounding_cap_of(&self.vertices)
    }

    fn bounding_rect(&self) -> Rect {
        self.rect_bound()
    }

    /// True iff the polyline may intersect `cell`: a vertex lies inside the
    /// cell, or an edge crosses *or touches* a cell edge
    /// ([`Crossing::Maybe`] counts, so boundary-coincident cells are kept —
    /// the conservative contract of [`Region`]).
    fn intersects_cell(&self, cell: S2CellId) -> bool {
        if self.vertices.is_empty() {
            return false;
        }
        let cl = cell_loop(cell);
        if self.vertices.iter().any(|&v| cl.contains_point(v)) {
            return true;
        }
        let cv = cl.vertices();
        let n = cv.len();
        for (a, b) in self.edges() {
            for i in 0..n {
                if crossing_sign(a, b, cv[i], cv[(i + 1) % n]) != Crossing::DoNotCross {
                    return true;
                }
            }
        }
        false
    }

    /// Always false: a one-dimensional region cannot contain a cell.
    fn contains_cell(&self, _cell: S2CellId) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::s2::{Region, RegionCoverer};

    fn pt(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    /// Points along the great circle from `a` to `b` (spherical interpolation),
    /// endpoints excluded.
    fn sample_edge(a: S2Point, b: S2Point, steps: u32) -> Vec<S2Point> {
        let theta = a.vector().angle(b.vector()).radians();
        (1..steps)
            .map(|k| {
                let t = f64::from(k) / f64::from(steps);
                let v = a
                    .vector()
                    .mul(((1.0 - t) * theta).sin())
                    .add(b.vector().mul((t * theta).sin()))
                    .normalize();
                S2Point::from_vector(v)
            })
            .collect()
    }

    #[test]
    fn intersects_cell_positive_and_negative() {
        let line = Polyline::new(vec![pt(0.0, 0.0), pt(0.0, 10.0)]);
        let through = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(0.0, 5.0), 8);
        assert!(
            line.intersects_cell(through),
            "edge passes through the cell"
        );
        let endpoint = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(0.0, 0.0), 8);
        assert!(line.intersects_cell(endpoint), "vertex lies in the cell");
        let far = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(40.0, 90.0), 8);
        assert!(!line.intersects_cell(far));
        assert!(!line.contains_cell(through));
        assert!(!line.contains_point(pt(0.0, 5.0)));
    }

    #[test]
    fn bounding_cap_contains_vertices_and_length_adds_up() {
        let line = Polyline::new(vec![pt(0.0, 0.0), pt(0.0, 90.0), pt(90.0, 0.0)]);
        let cap = line.bounding_cap();
        for &v in line.vertices() {
            assert!(cap.contains_point(v));
        }
        assert!(
            (line.length().radians() - std::f64::consts::PI).abs() < 1e-12,
            "two quarter arcs"
        );
        assert_eq!(line.num_vertices(), 3);
    }

    #[test]
    fn rect_bound_contains_sampled_edges() {
        let line = Polyline::new(vec![
            pt(45.0, 0.0),
            pt(45.0, 90.0),
            pt(10.0, 170.0),
            pt(10.0, -170.0),
        ]);
        let bound = line.rect_bound();
        assert_eq!(Region::bounding_rect(&line), bound);
        for w in line.vertices().windows(2) {
            for s in sample_edge(w[0], w[1], 200) {
                assert!(
                    bound.contains_lat_lng(s.to_geo_point()),
                    "{:?} escaped {bound:?}",
                    s.to_geo_point()
                );
            }
        }
    }

    #[test]
    fn fixed_level_cover_has_no_false_negatives_along_edges() {
        let line = Polyline::from_geo_points(&[
            GeoPoint::from_lat_lon(50.0, 9.5),
            GeoPoint::from_lat_lon(50.6, 10.4),
            GeoPoint::from_lat_lon(49.8, 11.2),
        ]);
        let level = 11;
        let cover = RegionCoverer::at_level(level).covering(&line);
        assert!(!cover.is_empty());
        for &c in &cover {
            assert_eq!(c.level(), level);
            assert!(line.intersects_cell(c));
        }
        for w in line.vertices().windows(2) {
            for s in sample_edge(w[0], w[1], 400) {
                let leaf = S2CellId::from_point(s);
                assert!(
                    cover.iter().any(|c| c.contains(leaf)),
                    "sampled edge point not covered"
                );
            }
        }
    }

    #[test]
    fn adaptive_cover_follows_the_line() {
        let line = Polyline::from_geo_points(&[
            GeoPoint::from_lat_lon(0.0, 0.0),
            GeoPoint::from_lat_lon(5.0, 20.0),
        ]);
        let cover = RegionCoverer {
            min_level: 0,
            max_level: 12,
            max_cells: 24,
        }
        .adaptive_covering(&line);
        assert!(!cover.is_empty());
        assert!(cover.len() <= 24);
        for s in sample_edge(line.vertices()[0], line.vertices()[1], 200) {
            assert!(cover.contains_point(s), "sampled point escaped the cover");
        }
    }
}
