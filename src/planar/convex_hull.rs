//! Planar [`ConvexHull`] via Andrew's monotone chain, mirroring `geo`'s
//! `ConvexHull` trait. Returns a closed, counter-clockwise [`Polygon`].

use super::{Coord, LineString, MultiPoint, MultiPolygon, Polygon};

/// The convex hull of a geometry, as a polygon that covers it.
pub trait ConvexHull {
    fn convex_hull(&self) -> Polygon;
}

/// Counter-clockwise convex hull of a coordinate set, returned as a closed ring
/// (first coordinate repeated at the end). Non-finite coordinates are dropped.
pub(crate) fn convex_hull_coords(points: &[Coord]) -> Vec<Coord> {
    let mut pts: Vec<Coord> = points.iter().copied().filter(|c| c.is_finite()).collect();
    pts.sort_by(|a, b| {
        a.x.partial_cmp(&b.x)
            .unwrap()
            .then(a.y.partial_cmp(&b.y).unwrap())
    });
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planar::area::Area;

    #[test]
    fn hull_of_square_with_interior_point() {
        // Four corners plus a point inside; the hull is just the square.
        let mp = MultiPoint::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (2.0, 2.0), // interior — should be dropped
        ]);
        let hull = mp.convex_hull();
        assert_eq!(hull.unsigned_area(), 16.0);
        assert!(hull.exterior().is_closed());
        // 4 distinct corners + closing point
        assert_eq!(hull.exterior().coords_count(), 5);
    }

    #[test]
    fn hull_is_counter_clockwise() {
        let mp = MultiPoint::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        let hull = mp.convex_hull();
        assert!(hull.signed_area() > 0.0); // CCW
    }
}
