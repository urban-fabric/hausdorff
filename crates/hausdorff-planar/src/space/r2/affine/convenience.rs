use super::*;

/// Translate a geometry along its axes by the given offsets.
pub trait Translate {
    #[must_use]
    fn translate(&self, x_offset: f64, y_offset: f64) -> Self;

    /// Mutable version of [`translate`](Self::translate).
    fn translate_mut(&mut self, x_offset: f64, y_offset: f64);
}

impl<G: AffineOps> Translate for G {
    fn translate(&self, x_offset: f64, y_offset: f64) -> Self {
        self.affine_transform(&AffineTransform::translate(x_offset, y_offset))
    }

    fn translate_mut(&mut self, x_offset: f64, y_offset: f64) {
        self.affine_transform_mut(&AffineTransform::translate(x_offset, y_offset));
    }
}

/// Rotate a geometry around a point by an angle, in degrees. Positive angles
/// are counter-clockwise.
pub trait Rotate {
    /// Rotate around the geometry's [centroid](Centroid); empty geometries are
    /// returned unchanged.
    #[must_use]
    fn rotate_around_centroid(&self, degrees: f64) -> Self;

    /// Mutable version of [`rotate_around_centroid`](Self::rotate_around_centroid).
    fn rotate_around_centroid_mut(&mut self, degrees: f64);

    /// Rotate around the centre of the geometry's [bounding
    /// rect](BoundingRect); empty geometries are returned unchanged.
    #[must_use]
    fn rotate_around_center(&self, degrees: f64) -> Self;

    /// Mutable version of [`rotate_around_center`](Self::rotate_around_center).
    fn rotate_around_center_mut(&mut self, degrees: f64);

    /// Rotate around an arbitrary point.
    #[must_use]
    fn rotate_around_point(&self, degrees: f64, point: Point) -> Self;

    /// Mutable version of [`rotate_around_point`](Self::rotate_around_point).
    fn rotate_around_point_mut(&mut self, degrees: f64, point: Point);
}

impl<G> Rotate for G
where
    G: Clone + AffineOps + Centroid + BoundingRect,
    <G as BoundingRect>::Output: Into<Option<Rect>>,
{
    fn rotate_around_centroid(&self, degrees: f64) -> Self {
        match self.centroid() {
            Some(point) => self.rotate_around_point(degrees, point),
            None => self.clone(),
        }
    }

    fn rotate_around_centroid_mut(&mut self, degrees: f64) {
        if let Some(point) = self.centroid() {
            self.rotate_around_point_mut(degrees, point);
        }
    }

    fn rotate_around_center(&self, degrees: f64) -> Self {
        match self.bounding_rect().into() {
            Some(rect) => self.rotate_around_point(degrees, rect.center()),
            None => self.clone(),
        }
    }

    fn rotate_around_center_mut(&mut self, degrees: f64) {
        if let Some(rect) = self.bounding_rect().into() {
            self.rotate_around_point_mut(degrees, rect.center());
        }
    }

    fn rotate_around_point(&self, degrees: f64, point: Point) -> Self {
        self.affine_transform(&AffineTransform::rotate(degrees, point))
    }

    fn rotate_around_point_mut(&mut self, degrees: f64, point: Point) {
        self.affine_transform_mut(&AffineTransform::rotate(degrees, point));
    }
}

/// Scale a geometry up or down by factors, by default from its bounding-box
/// centre.
pub trait Scale {
    /// Scale uniformly from the bounding-box centre; empty geometries are
    /// returned unchanged.
    #[must_use]
    fn scale(&self, scale_factor: f64) -> Self;

    /// Mutable version of [`scale`](Self::scale).
    fn scale_mut(&mut self, scale_factor: f64);

    /// Scale from the bounding-box centre with separate `x` / `y` factors;
    /// empty geometries are returned unchanged.
    #[must_use]
    fn scale_xy(&self, x_factor: f64, y_factor: f64) -> Self;

    /// Mutable version of [`scale_xy`](Self::scale_xy).
    fn scale_xy_mut(&mut self, x_factor: f64, y_factor: f64);

    /// Scale around an arbitrary point of `origin`.
    #[must_use]
    fn scale_around_point(&self, x_factor: f64, y_factor: f64, origin: impl Into<Coord>) -> Self;

    /// Mutable version of [`scale_around_point`](Self::scale_around_point).
    fn scale_around_point_mut(&mut self, x_factor: f64, y_factor: f64, origin: impl Into<Coord>);
}

impl<G> Scale for G
where
    G: Clone + AffineOps + BoundingRect,
    <G as BoundingRect>::Output: Into<Option<Rect>>,
{
    fn scale(&self, scale_factor: f64) -> Self {
        self.scale_xy(scale_factor, scale_factor)
    }

    fn scale_mut(&mut self, scale_factor: f64) {
        self.scale_xy_mut(scale_factor, scale_factor);
    }

    fn scale_xy(&self, x_factor: f64, y_factor: f64) -> Self {
        match self.bounding_rect().into() {
            Some(rect) => self.scale_around_point(x_factor, y_factor, rect.center()),
            None => self.clone(),
        }
    }

    fn scale_xy_mut(&mut self, x_factor: f64, y_factor: f64) {
        if let Some(rect) = self.bounding_rect().into() {
            self.scale_around_point_mut(x_factor, y_factor, rect.center());
        }
    }

    fn scale_around_point(&self, x_factor: f64, y_factor: f64, origin: impl Into<Coord>) -> Self {
        self.affine_transform(&AffineTransform::scale(x_factor, y_factor, origin))
    }

    fn scale_around_point_mut(&mut self, x_factor: f64, y_factor: f64, origin: impl Into<Coord>) {
        self.affine_transform_mut(&AffineTransform::scale(x_factor, y_factor, origin));
    }
}

/// Skew (shear) a geometry by angles in degrees along the `x` and `y`
/// dimensions, by default from its bounding-box centre.
pub trait Skew {
    /// Shear by a uniform angle along both dimensions; empty geometries are
    /// returned unchanged.
    #[must_use]
    fn skew(&self, degrees: f64) -> Self;

    /// Mutable version of [`skew`](Self::skew).
    fn skew_mut(&mut self, degrees: f64);

    /// Shear by separate angles along the `x` and `y` dimensions; empty
    /// geometries are returned unchanged.
    #[must_use]
    fn skew_xy(&self, degrees_x: f64, degrees_y: f64) -> Self;

    /// Mutable version of [`skew_xy`](Self::skew_xy).
    fn skew_xy_mut(&mut self, degrees_x: f64, degrees_y: f64);

    /// Shear around an arbitrary point of `origin`.
    #[must_use]
    fn skew_around_point(&self, degrees_x: f64, degrees_y: f64, origin: impl Into<Coord>) -> Self;

    /// Mutable version of [`skew_around_point`](Self::skew_around_point).
    fn skew_around_point_mut(&mut self, degrees_x: f64, degrees_y: f64, origin: impl Into<Coord>);
}

impl<G> Skew for G
where
    G: Clone + AffineOps + BoundingRect,
    <G as BoundingRect>::Output: Into<Option<Rect>>,
{
    fn skew(&self, degrees: f64) -> Self {
        self.skew_xy(degrees, degrees)
    }

    fn skew_mut(&mut self, degrees: f64) {
        self.skew_xy_mut(degrees, degrees);
    }

    fn skew_xy(&self, degrees_x: f64, degrees_y: f64) -> Self {
        match self.bounding_rect().into() {
            Some(rect) => self.skew_around_point(degrees_x, degrees_y, rect.center()),
            None => self.clone(),
        }
    }

    fn skew_xy_mut(&mut self, degrees_x: f64, degrees_y: f64) {
        if let Some(rect) = self.bounding_rect().into() {
            self.skew_around_point_mut(degrees_x, degrees_y, rect.center());
        }
    }

    fn skew_around_point(&self, degrees_x: f64, degrees_y: f64, origin: impl Into<Coord>) -> Self {
        self.affine_transform(&AffineTransform::skew(degrees_x, degrees_y, origin))
    }

    fn skew_around_point_mut(&mut self, degrees_x: f64, degrees_y: f64, origin: impl Into<Coord>) {
        self.affine_transform_mut(&AffineTransform::skew(degrees_x, degrees_y, origin));
    }
}
