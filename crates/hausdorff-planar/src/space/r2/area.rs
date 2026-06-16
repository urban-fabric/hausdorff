//! Planar [`Area`] — signed and unsigned shoelace area, mirroring `geo`'s
//! `Area` trait. Areas are in the (unspecified) planar unit²; for metric
//! geographic area use the spherical helpers after projecting, not this.

use super::util::signed_area2;
use super::{LineString, MultiPolygon, Point, Polygon, Rect, Triangle};

/// Signed and unsigned planar area.
///
/// Proof-backed by `riemann` `ExactArith/Shoelace` and `ExactArith/PlanarArea`:
/// the shoelace sum is the triangle area (`shoelace_triangle`,
/// `triangle_shoelace_eq_det`), translation-invariant (`shoelace_translation` /
/// `triSignedArea_translation`) and sign-flips on vertex swap
/// (`triSignedArea_swap`) — so signed area's sign is the orientation.
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
