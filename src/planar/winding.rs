//! Ring [`WindingOrder`] detection ([`Winding`]) and polygon reorientation
//! ([`Orient`]), mirroring `geo`'s `Winding` / `Orient` traits.

use super::util::signed_area2;
use super::{LineString, MultiPolygon, Polygon};

/// How a ring is wound.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum WindingOrder {
    Clockwise,
    CounterClockwise,
}

/// Determine the winding order of a ring.
pub trait Winding {
    /// The winding order, or `None` if the ring has no signed area (fewer than
    /// three distinct vertices, or collinear/degenerate).
    fn winding_order(&self) -> Option<WindingOrder>;

    /// True iff wound clockwise.
    fn is_cw(&self) -> bool {
        self.winding_order() == Some(WindingOrder::Clockwise)
    }

    /// True iff wound counter-clockwise.
    fn is_ccw(&self) -> bool {
        self.winding_order() == Some(WindingOrder::CounterClockwise)
    }
}

impl Winding for LineString {
    fn winding_order(&self) -> Option<WindingOrder> {
        let a2 = signed_area2(&self.0);
        if a2 > 0.0 {
            Some(WindingOrder::CounterClockwise)
        } else if a2 < 0.0 {
            Some(WindingOrder::Clockwise)
        } else {
            None
        }
    }
}

/// Target orientation for [`Orient::orient`].
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Direction {
    /// Exterior counter-clockwise, interiors clockwise (the OGC / "right-hand
    /// rule" convention, matching `geo`'s `Direction::Default`).
    Default,
    /// The reverse: exterior clockwise, interiors counter-clockwise.
    Reversed,
}

/// Reorient a polygon's rings to a chosen [`Direction`].
pub trait Orient {
    fn orient(&self, direction: Direction) -> Self;
}

fn ring_to_winding(ring: &LineString, target: WindingOrder) -> LineString {
    match ring.winding_order() {
        Some(w) if w != target => {
            let mut v = ring.0.clone();
            v.reverse(); // reversing a closed ring keeps it closed
            LineString::new(v)
        }
        _ => ring.clone(),
    }
}

impl Orient for Polygon {
    fn orient(&self, direction: Direction) -> Polygon {
        let (ext, int) = match direction {
            Direction::Default => (WindingOrder::CounterClockwise, WindingOrder::Clockwise),
            Direction::Reversed => (WindingOrder::Clockwise, WindingOrder::CounterClockwise),
        };
        Polygon::new(
            ring_to_winding(self.exterior(), ext),
            self.interiors()
                .iter()
                .map(|h| ring_to_winding(h, int))
                .collect(),
        )
    }
}

impl Orient for MultiPolygon {
    fn orient(&self, direction: Direction) -> MultiPolygon {
        MultiPolygon::new(self.iter().map(|p| p.orient(direction)).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planar::Coord;

    fn ls(coords: &[(f64, f64)]) -> LineString {
        LineString::new(coords.iter().map(|&(x, y)| Coord::new(x, y)).collect())
    }

    #[test]
    fn winding_detection() {
        // CCW square
        let ccw = ls(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]);
        assert_eq!(ccw.winding_order(), Some(WindingOrder::CounterClockwise));
        assert!(ccw.is_ccw());

        let mut rev = ccw.0.clone();
        rev.reverse();
        assert_eq!(LineString::new(rev).winding_order(), Some(WindingOrder::Clockwise));
    }

    #[test]
    fn orient_default_makes_exterior_ccw() {
        // Start with a CW exterior.
        let cw = Polygon::new(
            ls(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]),
            vec![],
        );
        assert!(cw.exterior().is_cw());
        let oriented = cw.orient(Direction::Default);
        assert!(oriented.exterior().is_ccw());
    }
}
