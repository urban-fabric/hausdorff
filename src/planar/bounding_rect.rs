//! Axis-aligned [`BoundingRect`] of a planar geometry, mirroring `geo`'s
//! `BoundingRect` trait (with its `Output` distinction: always-present for a
//! single coordinate, `Option` for collections that may be empty).

use super::{
    Coord, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect,
};

/// The bounding rectangle of a geometry.
pub trait BoundingRect {
    /// `Rect` for geometries that always have one (a point); `Option<Rect>` for
    /// collections that may be empty.
    type Output;

    fn bounding_rect(&self) -> Self::Output;
}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polygon_bounding_rect() {
        let p = Polygon::new(
            LineString::from(vec![(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (1.0, 5.0)]),
            vec![],
        );
        let r = p.bounding_rect().unwrap();
        assert_eq!(r.min(), Coord::new(0.0, 0.0));
        assert_eq!(r.max(), Coord::new(3.0, 5.0));
    }

    #[test]
    fn point_bounding_rect_is_unconditional() {
        let r: Rect = Point::new(2.0, 4.0).bounding_rect();
        assert_eq!(r.min(), Coord::new(2.0, 4.0));
        assert_eq!(r.max(), Coord::new(2.0, 4.0));
    }

    #[test]
    fn empty_line_string_has_no_rect() {
        assert!(LineString::default().bounding_rect().is_none());
    }
}
