//! Planar containment predicates ([`Contains`]), mirroring `geo`'s `Contains`
//! trait for the point-in-area tests the workspace uses.
//!
//! Containment is *interior* membership via an even-odd ray cast: a polygon
//! contains a point when it falls inside the exterior ring and outside every
//! hole. Results for points exactly on a boundary are not guaranteed.

use super::util::point_in_ring;
use super::{Coord, Geometry, MultiPolygon, Point, Polygon, Rect};

/// Whether `self` contains `rhs`.
pub trait Contains<Rhs = Self> {
    fn contains(&self, rhs: &Rhs) -> bool;
}

impl Contains<Coord> for Polygon {
    fn contains(&self, rhs: &Coord) -> bool {
        point_in_ring(*rhs, &self.exterior().0)
            && !self.interiors().iter().any(|h| point_in_ring(*rhs, &h.0))
    }
}

impl Contains<Point> for Polygon {
    fn contains(&self, rhs: &Point) -> bool {
        self.contains(&rhs.0)
    }
}

impl Contains<Coord> for MultiPolygon {
    fn contains(&self, rhs: &Coord) -> bool {
        self.iter().any(|p| p.contains(rhs))
    }
}

impl Contains<Point> for MultiPolygon {
    fn contains(&self, rhs: &Point) -> bool {
        self.contains(&rhs.0)
    }
}

impl Contains<Coord> for Rect {
    fn contains(&self, rhs: &Coord) -> bool {
        rhs.x > self.min().x && rhs.x < self.max().x && rhs.y > self.min().y && rhs.y < self.max().y
    }
}

impl Contains<Point> for Rect {
    fn contains(&self, rhs: &Point) -> bool {
        self.contains(&rhs.0)
    }
}

/// Areal interior containment lifted to the [`Geometry`] enum: a coordinate is
/// contained iff the geometry is areal (`Polygon`/`MultiPolygon`/`Rect`) and the
/// coordinate lies in its interior. Non-areal geometries contain nothing under
/// this predicate — use [`relate`](super::relate) for full DE-9IM containment
/// (e.g. a line containing a point on it).
impl Contains<Coord> for Geometry {
    fn contains(&self, rhs: &Coord) -> bool {
        match self {
            Geometry::Polygon(p) => p.contains(rhs),
            Geometry::MultiPolygon(mp) => mp.contains(rhs),
            Geometry::Rect(r) => r.contains(rhs),
            _ => false,
        }
    }
}

impl Contains<Point> for Geometry {
    fn contains(&self, rhs: &Point) -> bool {
        self.contains(&rhs.0)
    }
}
