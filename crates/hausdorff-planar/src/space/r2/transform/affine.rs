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

use crate::space::r2::{
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

mod affine_ops;
pub use affine_ops::AffineOps;

mod convenience;
pub use convenience::{Rotate, Scale, Skew, Translate};

#[cfg(test)]
mod tests;
