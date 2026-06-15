//! Planar [`ConvexHull`] via Andrew's monotone chain, mirroring `geo`'s
//! `ConvexHull` trait. Returns a closed, counter-clockwise [`Polygon`].
//!
//! **Minimality is proof-backed** by `riemann` `ExactArith/ConvexHull`: the hull
//! is the *smallest* convex set containing the input vertices
//! (`inHull_minimal` — any convex region containing the points contains the whole
//! hull; `inHull_subset_leftOfAll` at the half-plane region). That universal
//! property is the correct statement: the naïve "all left turns ⟹ convex" is
//! false without simplicity (a pentagram turns left at every vertex yet is not
//! convex), which is why the chain algorithm discards reflex/collinear vertices.

use super::{Coord, LineString, MultiPoint, MultiPolygon, Polygon};

/// The convex hull of a geometry, as a polygon that covers it.
pub trait ConvexHull {
    fn convex_hull(&self) -> Polygon;
}

/// Counter-clockwise convex hull of a coordinate set, returned as a closed ring
/// (first coordinate repeated at the end). Non-finite coordinates are dropped.
pub(crate) fn convex_hull_coords(points: &[Coord]) -> Vec<Coord> {
    let mut pts: Vec<Coord> = points.iter().copied().filter(|c| c.is_finite()).collect();
    pts.sort_unstable_by(|a, b| a.x.total_cmp(&b.x).then_with(|| a.y.total_cmp(&b.y)));
    pts.dedup();

    if pts.len() < 3 {
        if let Some(&first) = pts.first() {
            pts.push(first);
        }
        return pts;
    }

    let cross = |o: Coord, a: Coord, b: Coord| (a - o).x * (b - o).y - (a - o).y * (b - o).x;

    let mut lower: Vec<Coord> = Vec::new();
    for &p in &pts {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 0.0 {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<Coord> = Vec::new();
    for &p in pts.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 0.0 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    if let Some(&first) = lower.first() {
        lower.push(first);
    }
    lower
}

fn hull_polygon(points: &[Coord]) -> Polygon {
    Polygon::new(LineString::new(convex_hull_coords(points)), vec![])
}

impl ConvexHull for MultiPoint {
    fn convex_hull(&self) -> Polygon {
        let coords: Vec<Coord> = self.0.iter().map(|p| p.0).collect();
        hull_polygon(&coords)
    }
}

impl ConvexHull for LineString {
    fn convex_hull(&self) -> Polygon {
        hull_polygon(&self.0)
    }
}

impl ConvexHull for Polygon {
    fn convex_hull(&self) -> Polygon {
        hull_polygon(&self.exterior().0)
    }
}

impl ConvexHull for MultiPolygon {
    fn convex_hull(&self) -> Polygon {
        let coords: Vec<Coord> = self
            .iter()
            .flat_map(|p| p.exterior().0.iter().copied())
            .collect();
        hull_polygon(&coords)
    }
}
