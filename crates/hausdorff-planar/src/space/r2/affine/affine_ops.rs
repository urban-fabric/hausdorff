use super::*;

/// Apply an [`AffineTransform`] to a geometry, immutably or in place.
pub trait AffineOps {
    /// Apply `transform`, outputting a new geometry.
    #[must_use]
    fn affine_transform(&self, transform: &AffineTransform) -> Self;

    /// Apply `transform` to mutate `self`.
    fn affine_transform_mut(&mut self, transform: &AffineTransform);
}

impl AffineOps for Point {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        Point(transform.apply(self.0))
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        self.0 = transform.apply(self.0);
    }
}

impl AffineOps for Line {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        Line::new(transform.apply(self.start), transform.apply(self.end))
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        self.start = transform.apply(self.start);
        self.end = transform.apply(self.end);
    }
}

impl AffineOps for LineString {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        LineString::new(self.0.iter().map(|&c| transform.apply(c)).collect())
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        for c in &mut self.0 {
            *c = transform.apply(*c);
        }
    }
}

impl AffineOps for Polygon {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        Polygon::new(
            self.exterior().affine_transform(transform),
            self.interiors()
                .iter()
                .map(|r| r.affine_transform(transform))
                .collect(),
        )
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        *self = self.affine_transform(transform);
    }
}

impl AffineOps for MultiPoint {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        MultiPoint::new(self.iter().map(|p| p.affine_transform(transform)).collect())
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        for p in &mut self.0 {
            p.affine_transform_mut(transform);
        }
    }
}

impl AffineOps for MultiLineString {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        MultiLineString::new(
            self.iter()
                .map(|ls| ls.affine_transform(transform))
                .collect(),
        )
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        for ls in &mut self.0 {
            ls.affine_transform_mut(transform);
        }
    }
}

impl AffineOps for MultiPolygon {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        MultiPolygon::new(self.iter().map(|p| p.affine_transform(transform)).collect())
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        for p in &mut self.0 {
            p.affine_transform_mut(transform);
        }
    }
}

impl AffineOps for Rect {
    /// The transformed corners renormalize to an axis-aligned rectangle, the
    /// same min/max-corner mapping as `geo`'s `Rect` under `MapCoords`.
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        Rect::new(transform.apply(self.min()), transform.apply(self.max()))
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        *self = self.affine_transform(transform);
    }
}

impl AffineOps for Triangle {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        Triangle::new(
            transform.apply(self.0),
            transform.apply(self.1),
            transform.apply(self.2),
        )
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        *self = self.affine_transform(transform);
    }
}

impl AffineOps for Geometry {
    fn affine_transform(&self, transform: &AffineTransform) -> Self {
        match self {
            Geometry::Point(g) => Geometry::Point(g.affine_transform(transform)),
            Geometry::Line(g) => Geometry::Line(g.affine_transform(transform)),
            Geometry::LineString(g) => Geometry::LineString(g.affine_transform(transform)),
            Geometry::Polygon(g) => Geometry::Polygon(g.affine_transform(transform)),
            Geometry::MultiPoint(g) => Geometry::MultiPoint(g.affine_transform(transform)),
            Geometry::MultiLineString(g) => {
                Geometry::MultiLineString(g.affine_transform(transform))
            }
            Geometry::MultiPolygon(g) => Geometry::MultiPolygon(g.affine_transform(transform)),
            Geometry::Rect(g) => Geometry::Rect(g.affine_transform(transform)),
            Geometry::Triangle(g) => Geometry::Triangle(g.affine_transform(transform)),
            Geometry::GeometryCollection(gs) => Geometry::GeometryCollection(
                gs.iter().map(|g| g.affine_transform(transform)).collect(),
            ),
        }
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform) {
        *self = self.affine_transform(transform);
    }
}
