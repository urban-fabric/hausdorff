//! Axis-aligned [`BoundingRect`] of a planar geometry, mirroring `geo`'s
//! `BoundingRect` trait (with its `Output` distinction: always-present for a
//! single coordinate, `Option` for collections that may be empty).

use crate::space::r2::{
    Coord, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect,
    Triangle,
};

/// The bounding rectangle of a geometry.
pub trait BoundingRect {
    /// `Rect` for geometries that always have one (a point); `Option<Rect>` for
    /// collections that may be empty.
    type Output;

    fn bounding_rect(&self) -> Self::Output;
}

/// The axis-aligned bounding box of a coordinate set: `[min x, max x] × [min y,
/// max y]`, or `None` if empty.
///
/// Proof-backed by `riemann` `ExactArith/BoundingRect`: the running-`min`/`max`
/// fold is **sound** — every input point lies in the box (`boundingRect_contains`)
/// — and **minimal** — any box containing all the points encloses it
/// (`boundingRect_minimal`), i.e. it is the smallest enclosing box (and so the
/// R-tree `leaf ⊆ node` invariant feeding `AABB.prune_sound`).
fn rect_of(coords: impl IntoIterator<Item = Coord>) -> Option<Rect> {
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut any = false;
    for c in coords {
        if !c.is_finite() {
            continue;
        }
        any = true;
        min_x = min_x.min(c.x);
        min_y = min_y.min(c.y);
        max_x = max_x.max(c.x);
        max_y = max_y.max(c.y);
    }
    any.then(|| Rect::new(Coord::new(min_x, min_y), Coord::new(max_x, max_y)))
}

impl BoundingRect for Coord {
    type Output = Rect;
    fn bounding_rect(&self) -> Rect {
        Rect::new(*self, *self)
    }
}

impl BoundingRect for Point {
    type Output = Rect;
    fn bounding_rect(&self) -> Rect {
        Rect::new(self.0, self.0)
    }
}

impl BoundingRect for Rect {
    type Output = Rect;
    fn bounding_rect(&self) -> Rect {
        *self
    }
}

impl BoundingRect for Line {
    type Output = Rect;
    fn bounding_rect(&self) -> Rect {
        Rect::new(self.start, self.end)
    }
}

impl BoundingRect for Triangle {
    type Output = Rect;
    fn bounding_rect(&self) -> Rect {
        Rect::new(
            Coord::new(
                self.0.x.min(self.1.x).min(self.2.x),
                self.0.y.min(self.1.y).min(self.2.y),
            ),
            Coord::new(
                self.0.x.max(self.1.x).max(self.2.x),
                self.0.y.max(self.1.y).max(self.2.y),
            ),
        )
    }
}

impl BoundingRect for LineString {
    type Output = Option<Rect>;
    fn bounding_rect(&self) -> Option<Rect> {
        rect_of(self.coords_copied())
    }
}

impl BoundingRect for MultiLineString {
    type Output = Option<Rect>;
    fn bounding_rect(&self) -> Option<Rect> {
        rect_of(self.iter().flat_map(|ls| ls.coords_copied()))
    }
}

impl BoundingRect for Polygon {
    type Output = Option<Rect>;
    fn bounding_rect(&self) -> Option<Rect> {
        rect_of(self.exterior().coords_copied())
    }
}

impl BoundingRect for MultiPolygon {
    type Output = Option<Rect>;
    fn bounding_rect(&self) -> Option<Rect> {
        rect_of(self.iter().flat_map(|p| p.exterior().coords_copied()))
    }
}

impl BoundingRect for MultiPoint {
    type Output = Option<Rect>;
    fn bounding_rect(&self) -> Option<Rect> {
        rect_of(self.iter().map(|p| p.0))
    }
}
