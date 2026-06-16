//! Planar affine transformations — mirror of `geo`'s [`AffineTransform`]
//! matrix plus the convenience traits [`AffineOps`], [`Translate`], [`Rotate`],
//! [`Scale`], and [`Skew`].
//!
//! An [`AffineTransform`] is a row-major 3×3 matrix whose last row is fixed at
//! `[0, 0, 1]`. Transforms **compose**, so an arbitrary chain of operations
//! collapses to a single matrix applied once per coordinate. Constructors use
//! the present tense (`translate`, `scale`, `rotate`, `skew`); the past
//! participle (`translated`, `scaled`, …) *adds* the operation to an existing
//! transform. All angles are in degrees; positive rotation angles are
//! counter-clockwise.

use std::fmt;

use super::{
    BoundingRect, Centroid, Coord, Geometry, Line, LineString, MultiLineString, MultiPoint,
    MultiPolygon, Point, Polygon, Rect, Triangle,
};

/// A general affine transformation matrix mapping `(x, y)` to
/// `(a·x + b·y + xoff, d·x + e·y + yoff)` — row-major, with the rows ordered
/// `[[a, b, xoff], [d, e, yoff], [0, 0, 1]]`.
///
/// Proof-backed by `riemann` `ExactArith/Affine`: affine maps scale the
/// orientation determinant by `det` (`triDet_affMap`, with `triDet_affMap_pos`
/// when `det > 0`) — so a positive-determinant transform preserves orientation —
/// and preserve collinearity (`collinear_affMap`); rotations preserve it exactly
/// (`triDet_rotation`).
#[derive(Copy, Clone, PartialEq)]
pub struct AffineTransform([[f64; 3]; 3]);

impl Default for AffineTransform {
    fn default() -> Self {
        Self::identity()
    }
}

impl AffineTransform {
    /// A custom transform from the six free matrix parameters
    /// `[[a, b, xoff], [d, e, yoff]]`.
    #[must_use]
    pub fn new(a: f64, b: f64, xoff: f64, d: f64, e: f64, yoff: f64) -> Self {
        Self([[a, b, xoff], [d, e, yoff], [0.0, 0.0, 1.0]])
    }

    /// The identity matrix: applying it is a no-op.
    #[must_use]
    pub fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    }

    /// Whether this transform is the [identity](Self::identity).
    #[must_use]
    pub fn is_identity(&self) -> bool {
        self == &Self::identity()
    }

    /// Compose this transform with `other`, so `other` runs *after* `self`.
    /// This is a **cumulative** operation; the new transform is *added* to the
    /// existing one.
    #[must_use]
    pub fn compose(&self, other: &Self) -> Self {
        let mut out = [[0.0; 3]; 3];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = (0..3).map(|k| other.0[i][k] * self.0[k][j]).sum();
            }
        }
        Self(out)
    }

    /// Compose this transform with an arbitrary number of others, applied in
    /// order after `self`.
    #[must_use]
    pub fn compose_many(&self, transforms: &[Self]) -> Self {
        self.compose(
            &transforms
                .iter()
                .fold(Self::identity(), |acc, transform| acc.compose(transform)),
        )
    }

    /// Apply the transform to a coordinate.
    #[must_use]
    pub fn apply(&self, coord: Coord) -> Coord {
        Coord::new(
            self.0[0][0] * coord.x + self.0[0][1] * coord.y + self.0[0][2],
            self.0[1][0] * coord.x + self.0[1][1] * coord.y + self.0[1][2],
        )
    }

    /// **Create** a translation by offsets along the `x` and `y` dimensions.
    #[must_use]
    pub fn translate(xoff: f64, yoff: f64) -> Self {
        Self::new(1.0, 0.0, xoff, 0.0, 1.0, yoff)
    }

    /// **Add** a translation to this transform (cumulative).
    #[must_use]
    pub fn translated(mut self, xoff: f64, yoff: f64) -> Self {
        self.0 = self.compose(&Self::translate(xoff, yoff)).0;
        self
    }

    /// **Create** a scaling by factors along the `x` and `y` dimensions about
    /// `origin` (usually the geometry's bounding-box centre). Negative factors
    /// mirror or reflect.
    #[must_use]
    pub fn scale(xfact: f64, yfact: f64, origin: impl Into<Coord>) -> Self {
        let (x0, y0) = origin.into().x_y();
        let xoff = x0 - (x0 * xfact);
        let yoff = y0 - (y0 * yfact);
        Self::new(xfact, 0.0, xoff, 0.0, yfact, yoff)
    }

    /// **Add** a scaling about `origin` to this transform (cumulative).
    #[must_use]
    pub fn scaled(mut self, xfact: f64, yfact: f64, origin: impl Into<Coord>) -> Self {
        self.0 = self.compose(&Self::scale(xfact, yfact, origin)).0;
        self
    }

    /// **Create** a rotation by `degrees` (counter-clockwise when positive)
    /// about `origin`.
    #[must_use]
    pub fn rotate(degrees: f64, origin: impl Into<Coord>) -> Self {
        let (sin_theta, cos_theta) = degrees.to_radians().sin_cos();
        let (x0, y0) = origin.into().x_y();
        let xoff = x0 - (x0 * cos_theta) + (y0 * sin_theta);
        let yoff = y0 - (x0 * sin_theta) - (y0 * cos_theta);
        Self::new(cos_theta, -sin_theta, xoff, sin_theta, cos_theta, yoff)
    }

    /// **Add** a rotation about `origin` to this transform (cumulative).
    #[must_use]
    pub fn rotated(mut self, degrees: f64, origin: impl Into<Coord>) -> Self {
        self.0 = self.compose(&Self::rotate(degrees, origin)).0;
        self
    }

    /// **Create** a shear by angles (in degrees) along the `x` (`xs`) and `y`
    /// (`ys`) dimensions about `origin` (usually the geometry's bounding-box
    /// centre).
    #[must_use]
    pub fn skew(xs: f64, ys: f64, origin: impl Into<Coord>) -> Self {
        let (x0, y0) = origin.into().x_y();
        let mut tanx = xs.to_radians().tan();
        let mut tany = ys.to_radians().tan();
        if tanx.abs() < 2.5e-16 {
            tanx = 0.0;
        }
        if tany.abs() < 2.5e-16 {
            tany = 0.0;
        }
        let xoff = -y0 * tanx;
        let yoff = -x0 * tany;
        Self::new(1.0, tanx, xoff, tany, 1.0, yoff)
    }

    /// **Add** a shear about `origin` to this transform (cumulative).
    #[must_use]
    pub fn skewed(mut self, xs: f64, ys: f64, origin: impl Into<Coord>) -> Self {
        self.0 = self.compose(&Self::skew(xs, ys, origin)).0;
        self
    }

    /// The inverse transform, or `None` when the matrix is singular. Composing
    /// a transform with its inverse yields the [identity](Self::identity).
    #[must_use]
    pub fn inverse(&self) -> Option<Self> {
        let (a, b, xoff) = (self.0[0][0], self.0[0][1], self.0[0][2]);
        let (d, e, yoff) = (self.0[1][0], self.0[1][1], self.0[1][2]);
        let determinant = a * e - b * d;
        if determinant == 0.0 {
            return None;
        }
        let inv_det = 1.0 / determinant;
        Some(Self::new(
            e * inv_det,
            -b * inv_det,
            (b * yoff - e * xoff) * inv_det,
            -d * inv_det,
            a * inv_det,
            (d * xoff - a * yoff) * inv_det,
        ))
    }

    /// Matrix parameter `a` (row 0, column 0).
    #[must_use]
    pub fn a(&self) -> f64 {
        self.0[0][0]
    }

    /// Matrix parameter `b` (row 0, column 1).
    #[must_use]
    pub fn b(&self) -> f64 {
        self.0[0][1]
    }

    /// Matrix parameter `xoff` (row 0, column 2).
    #[must_use]
    pub fn xoff(&self) -> f64 {
        self.0[0][2]
    }

    /// Matrix parameter `d` (row 1, column 0).
    #[must_use]
    pub fn d(&self) -> f64 {
        self.0[1][0]
    }

    /// Matrix parameter `e` (row 1, column 1).
    #[must_use]
    pub fn e(&self) -> f64 {
        self.0[1][1]
    }

    /// Matrix parameter `yoff` (row 1, column 2).
    #[must_use]
    pub fn yoff(&self) -> f64 {
        self.0[1][2]
    }

    /// Whether this transform is **rigid** (a Euclidean isometry — distance- and
    /// angle-preserving): true iff the linear part `[[a, b], [d, e]]` has
    /// orthonormal columns — `a² + d² ≈ 1`, `b² + e² ≈ 1`, `a·b + d·e ≈ 0`
    /// (within a `1e-9` tolerance). Orthonormal columns force `|det| = 1`, so
    /// this is the full rigid / distance-preserving test (rotations, reflections,
    /// translations and their composites pass; non-uniform scales and shears do
    /// not).
    ///
    /// Proof-backed by `riemann` `Hausdorff/AffineIsometry.act_isometry`
    /// (orthonormal-column maps preserve squared distance) and `rotation_isometry`
    /// (rotations are a special case).
    #[must_use]
    pub fn is_rigid(&self) -> bool {
        const TOL: f64 = 1e-9;
        let (a, b) = (self.0[0][0], self.0[0][1]);
        let (d, e) = (self.0[1][0], self.0[1][1]);
        (a * a + d * d - 1.0).abs() <= TOL
            && (b * b + e * e - 1.0).abs() <= TOL
            && (a * b + d * e).abs() <= TOL
    }
}

/// Affine transforms form a **group** under composition. `a * b` is
/// `a.compose(&b)` — applying `a` *then* `b`, the same order as
/// [`compose`](AffineTransform::compose). With [`identity`](AffineTransform::identity)
/// the unit and [`inverse`](AffineTransform::inverse) the (partial) inverse, this
/// is the group structure.
///
/// Proof-backed by `riemann` `Hausdorff/AffineMap`: `act_comp` (compose is
/// function composition — `(a*b).apply(p) == b.apply(a.apply(p))`), `comp_assoc`
/// (associativity), `id_comp`/`comp_id` (identity), and `comp_inv`/`inv_comp`
/// (inverses).
impl std::ops::Mul for AffineTransform {
    type Output = AffineTransform;

    /// Compose, applying `self` then `rhs` (same order as
    /// [`compose`](AffineTransform::compose)).
    fn mul(self, rhs: AffineTransform) -> AffineTransform {
        self.compose(&rhs)
    }
}

/// By-reference composition, so chains need not move/copy the matrices. Same
/// `self`-then-`rhs` order as [`AffineTransform::compose`].
impl std::ops::Mul<&AffineTransform> for &AffineTransform {
    type Output = AffineTransform;

    fn mul(self, rhs: &AffineTransform) -> AffineTransform {
        self.compose(rhs)
    }
}

/// In-place composition: `a *= b` makes `a` apply its old self then `b` (same
/// order as [`AffineTransform::compose`]).
impl std::ops::MulAssign for AffineTransform {
    fn mul_assign(&mut self, rhs: AffineTransform) {
        *self = self.compose(&rhs);
    }
}

impl fmt::Debug for AffineTransform {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AffineTransform")
            .field("a", &self.0[0][0])
            .field("b", &self.0[0][1])
            .field("xoff", &self.0[0][2])
            .field("d", &self.0[1][0])
            .field("e", &self.0[1][1])
            .field("yoff", &self.0[1][2])
            .finish()
    }
}

impl From<[f64; 6]> for AffineTransform {
    fn from(arr: [f64; 6]) -> Self {
        Self::new(arr[0], arr[1], arr[2], arr[3], arr[4], arr[5])
    }
}

impl From<(f64, f64, f64, f64, f64, f64)> for AffineTransform {
    fn from(tup: (f64, f64, f64, f64, f64, f64)) -> Self {
        Self::new(tup.0, tup.1, tup.2, tup.3, tup.4, tup.5)
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    fn assert_coords_eq(a: &LineString, b: &[(f64, f64)]) {
        assert_eq!(a.0.len(), b.len());
        for (c, &(x, y)) in a.0.iter().zip(b) {
            assert!(
                (c.x - x).abs() < EPS && (c.y - y).abs() < EPS,
                "({}, {}) vs ({x}, {y})",
                c.x,
                c.y
            );
        }
    }

    #[test]
    fn identity_and_is_identity() {
        let t = AffineTransform::identity();
        assert!(t.is_identity());
        assert!(!t.translated(1.0, 2.0).is_identity());
        assert!(t.translated(1.0, 2.0).translated(-1.0, -2.0).is_identity());
        assert_eq!(t.apply(Coord::new(3.0, -4.0)), Coord::new(3.0, -4.0));
    }

    #[test]
    fn compose_runs_in_order() {
        let scale = AffineTransform::scale(2.0, 2.0, Coord::new(0.0, 0.0));
        let translate = AffineTransform::translate(10.0, 0.0);
        let scale_then_translate = scale.compose(&translate);
        assert_eq!(
            scale_then_translate.apply(Coord::new(1.0, 1.0)),
            Coord::new(12.0, 2.0)
        );
        let translate_then_scale = translate.compose(&scale);
        assert_eq!(
            translate_then_scale.apply(Coord::new(1.0, 1.0)),
            Coord::new(22.0, 2.0)
        );
        let many = AffineTransform::identity().compose_many(&[scale, translate]);
        assert_eq!(many.apply(Coord::new(1.0, 1.0)), Coord::new(12.0, 2.0));
    }

    #[test]
    fn inverse_round_trips() {
        let t = AffineTransform::translate(3.0, -7.0)
            .rotated(31.0, Coord::new(1.0, 2.0))
            .scaled(2.0, 0.5, Coord::new(-3.0, 4.0));
        let inv = t.inverse().unwrap();
        let round = t.compose(&inv);
        let c = Coord::new(5.5, -2.25);
        let rc = round.apply(c);
        assert!((rc.x - c.x).abs() < EPS && (rc.y - c.y).abs() < EPS);
        assert!(
            AffineTransform::scale(0.0, 1.0, Coord::new(0.0, 0.0))
                .inverse()
                .is_none()
        );
    }

    #[test]
    fn translate_line_string() {
        let ls = LineString::from(vec![(0.0, 0.0), (5.0, 5.0), (10.0, 10.0)]);
        let translated = ls.translate(1.5, 3.5);
        assert_coords_eq(&translated, &[(1.5, 3.5), (6.5, 8.5), (11.5, 13.5)]);
        let mut in_place = ls;
        in_place.translate_mut(1.5, 3.5);
        assert_eq!(in_place, translated);
    }

    #[test]
    fn rotate_point_about_itself_is_fixed() {
        let p = Point::new(1.0, 5.0);
        assert_eq!(p.rotate_around_centroid(30.0), p);
        let rotated = p.rotate_around_point(90.0, Point::new(0.0, 0.0));
        assert!((rotated.x() + 5.0).abs() < EPS && (rotated.y() - 1.0).abs() < EPS);
    }

    #[test]
    fn scale_about_bounding_box_centre() {
        let ls = LineString::from(vec![(0.0, 0.0), (10.0, 10.0)]);
        assert_coords_eq(&ls.scale(2.0), &[(-5.0, -5.0), (15.0, 15.0)]);
        assert_coords_eq(&ls.scale_xy(2.0, 4.0), &[(-5.0, -15.0), (15.0, 25.0)]);
        assert_coords_eq(
            &ls.scale_around_point(2.0, 4.0, Coord::new(100.0, 100.0)),
            &[(-100.0, -300.0), (-80.0, -260.0)],
        );
    }

    #[test]
    fn skew_zero_is_identity() {
        let poly = Polygon::new(
            LineString::from(vec![(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
            vec![],
        );
        assert_eq!(poly.skew(0.0), poly);
        let skewed = poly.skew(30.0);
        assert!((skewed.exterior().0[0].x - -2.886751345948129).abs() < 1e-9);
    }

    #[test]
    fn mul_is_compose_and_a_group() {
        let a = AffineTransform::scale(2.0, 3.0, Coord::new(0.0, 0.0));
        let b = AffineTransform::translate(10.0, -4.0);
        let p = Coord::new(1.5, -2.0);
        // `a * b` applies a then b — identical to `b.apply(a.apply(p))`.
        assert_eq!((a * b).apply(p), b.apply(a.apply(p)));
        assert_eq!((a * b).apply(p), a.compose(&b).apply(p));
        // By-reference Mul agrees (called explicitly: the values are `Copy`, so
        // the `&a * &b` sugar would trip `clippy::op_ref`).
        assert_eq!(std::ops::Mul::mul(&a, &b).apply(p), (a * b).apply(p));
        // Identity is the unit.
        assert_eq!((a * AffineTransform::identity()).apply(p), a.apply(p));
        assert_eq!((AffineTransform::identity() * a).apply(p), a.apply(p));
        // MulAssign matches Mul.
        let mut acc = a;
        acc *= b;
        assert_eq!(acc.apply(p), (a * b).apply(p));
        // Inverse round-trips to the identity.
        let inv = a.inverse().unwrap();
        let round = (a * inv).apply(p);
        assert!((round.x - p.x).abs() < EPS && (round.y - p.y).abs() < EPS);
    }

    #[test]
    fn is_rigid_distinguishes_isometries() {
        // A rotation is rigid (orthonormal columns).
        assert!(AffineTransform::rotate(37.0, Coord::new(2.0, -1.0)).is_rigid());
        // A pure translation is rigid.
        assert!(AffineTransform::translate(5.0, -3.0).is_rigid());
        // The identity is rigid.
        assert!(AffineTransform::identity().is_rigid());
        // A reflection (negative uniform scale) is rigid (|det| = 1).
        assert!(AffineTransform::scale(-1.0, 1.0, Coord::new(0.0, 0.0)).is_rigid());
        // A non-uniform scale (x by 2) is NOT rigid.
        assert!(!AffineTransform::scale(2.0, 1.0, Coord::new(0.0, 0.0)).is_rigid());
        // A uniform scale ≠ 1 is not rigid (columns not unit length).
        assert!(!AffineTransform::scale(2.0, 2.0, Coord::new(0.0, 0.0)).is_rigid());
        // A shear is not rigid (columns not orthogonal).
        assert!(!AffineTransform::skew(30.0, 0.0, Coord::new(0.0, 0.0)).is_rigid());
        // A rotation composed with a translation stays rigid.
        let rt = AffineTransform::rotate(20.0, Coord::new(0.0, 0.0))
            * AffineTransform::translate(7.0, 2.0);
        assert!(rt.is_rigid());
    }

    #[test]
    fn empty_geometries_pass_through() {
        let empty = LineString::new(vec![]);
        assert_eq!(empty.scale(2.0), empty);
        assert_eq!(empty.skew(30.0), empty);
        assert_eq!(empty.rotate_around_centroid(45.0), empty);
        assert_eq!(empty.rotate_around_center(45.0), empty);
    }
}
