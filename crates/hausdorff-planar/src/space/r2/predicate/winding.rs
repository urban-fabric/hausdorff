//! Ring [`WindingOrder`] detection ([`Winding`]) and polygon reorientation
//! ([`Orient`]), mirroring `geo`'s `Winding` / `Orient` traits.

use crate::space::r2::util::signed_area2;
use crate::space::r2::{Geometry, LineString, MultiPolygon, Polygon};

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
            v.reverse();
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

impl Orient for Geometry {
    /// Reorient the areal members (`Polygon`/`MultiPolygon`, and recursively a
    /// `GeometryCollection`); other variants have no ring winding and are
    /// returned unchanged.
    fn orient(&self, direction: Direction) -> Geometry {
        match self {
            Geometry::Polygon(p) => Geometry::Polygon(p.orient(direction)),
            Geometry::MultiPolygon(mp) => Geometry::MultiPolygon(mp.orient(direction)),
            Geometry::GeometryCollection(gs) => {
                Geometry::GeometryCollection(gs.iter().map(|g| g.orient(direction)).collect())
            }
            other => other.clone(),
        }
    }
}
