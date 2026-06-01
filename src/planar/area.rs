//! Planar [`Area`] — signed and unsigned shoelace area, mirroring `geo`'s
//! `Area` trait. Areas are in the (unspecified) planar unit²; for metric
//! geographic area use the spherical helpers after projecting, not this.

use super::util::signed_area2;
use super::{LineString, MultiPolygon, Point, Polygon, Rect, Triangle};

/// Signed and unsigned planar area.
pub trait Area {
    /// Signed area: positive when the exterior ring is counter-clockwise,
    /// negative when clockwise. Holes reduce the magnitude.
    fn signed_area(&self) -> f64;

    /// Absolute area — `signed_area().abs()`.
    fn unsigned_area(&self) -> f64 {
        self.signed_area().abs()
    }
}

impl Area for Polygon {
    fn signed_area(&self) -> f64 {
        let ext = signed_area2(&self.exterior().0) / 2.0;
        let holes: f64 = self
            .interiors()
            .iter()
            .map(|h| (signed_area2(&h.0) / 2.0).abs())
            .sum();
        let magnitude = (ext.abs() - holes).max(0.0);
        if ext < 0.0 { -magnitude } else { magnitude }
    }
}

impl Area for MultiPolygon {
    fn signed_area(&self) -> f64 {
        self.iter().map(Area::signed_area).sum()
    }
}

impl Area for Rect {
    fn signed_area(&self) -> f64 {
        self.width() * self.height()
    }
}

impl Area for Triangle {
    fn signed_area(&self) -> f64 {
        signed_area2(&[self.0, self.1, self.2]) / 2.0
    }
}

impl Area for LineString {
    fn signed_area(&self) -> f64 {
        0.0
    }
}

impl Area for Point {
    fn signed_area(&self) -> f64 {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planar::Coord;

    fn poly(coords: &[(f64, f64)]) -> Polygon {
        Polygon::new(
            LineString::new(coords.iter().map(|&(x, y)| Coord::new(x, y)).collect()),
            vec![],
        )
    }

    #[test]
    fn unit_square_area() {
        let p = poly(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        assert_eq!(p.unsigned_area(), 1.0);
        assert_eq!(p.signed_area(), 1.0); // CCW -> positive
    }

    #[test]
    fn clockwise_is_negative_signed() {
        let p = poly(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]);
        assert_eq!(p.signed_area(), -1.0);
        assert_eq!(p.unsigned_area(), 1.0);
    }

    #[test]
    fn hole_subtracts_area() {
        let p = Polygon::new(
            LineString::new(
                [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
                    .iter()
                    .map(|&(x, y)| Coord::new(x, y))
                    .collect(),
            ),
            vec![LineString::new(
                [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]
                    .iter()
                    .map(|&(x, y)| Coord::new(x, y))
                    .collect(),
            )],
        );
        assert_eq!(p.unsigned_area(), 16.0 - 1.0);
    }

    #[test]
    fn rect_and_triangle() {
        assert_eq!(Rect::new((0.0, 0.0), (2.0, 3.0)).unsigned_area(), 6.0);
        let t = Triangle::new(Coord::new(0.0, 0.0), Coord::new(4.0, 0.0), Coord::new(0.0, 3.0));
        assert_eq!(t.unsigned_area(), 6.0);
    }
}
