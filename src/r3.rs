//! `r3` — three-dimensional Euclidean geometry: the [`Vector`] type, mirroring
//! golang/geo's `r3` package.
//!
//! `r3::Vector` is the substrate of the spherical world: an [`crate::s2::S2Point`]
//! is a *unit* `r3::Vector`, and spherical measures/predicates are computed from
//! vector dot/cross products rather than lat/lon, which is what makes them free
//! of pole and antimeridian singularities.

use crate::s1::Angle;

/// A vector in ℝ³.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// The three axes of ℝ³.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Vector {
    /// The vector `(x, y, z)`.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Dot (scalar) product.
    #[inline]
    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross (vector) product.
    #[inline]
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Vector sum.
    #[inline]
    #[must_use]
    pub fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }

    /// Vector difference.
    #[inline]
    #[must_use]
    pub fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }

    /// Scalar multiple.
    #[inline]
    #[must_use]
    pub fn mul(self, m: f64) -> Self {
        Self::new(self.x * m, self.y * m, self.z * m)
    }

    /// Squared Euclidean norm (`self · self`).
    #[inline]
    #[must_use]
    pub fn norm2(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm(self) -> f64 {
        self.norm2().sqrt()
    }

    /// A unit vector in the same direction, or the zero vector if `self` is zero
    /// (matching `r3.Vector.Normalize`).
    #[must_use]
    pub fn normalize(self) -> Self {
        let n2 = self.norm2();
        if n2 == 0.0 {
            return Self::new(0.0, 0.0, 0.0);
        }
        self.mul(1.0 / n2.sqrt())
    }

    /// Componentwise absolute value.
    #[must_use]
    pub fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// True iff the vector is (near) unit length.
    #[must_use]
    pub fn is_unit(self) -> bool {
        (self.norm2() - 1.0).abs() <= 1e-14
    }

    /// The angle between `self` and `o`, robust for both near-parallel and
    /// near-antiparallel vectors (`atan2(|a×b|, a·b)`).
    #[must_use]
    pub fn angle(self, o: Self) -> Angle {
        Angle::from_radians(self.cross(o).norm().atan2(self.dot(o)))
    }

    /// The axis of the largest-magnitude component. Ties resolve toward `Z` then
    /// `Y` (strict `>`), exactly as `r3.Vector.LargestComponent`.
    #[must_use]
    pub fn largest_component(self) -> Axis {
        let (tx, ty, tz) = (self.x.abs(), self.y.abs(), self.z.abs());
        if tx > ty {
            if tx > tz { Axis::X } else { Axis::Z }
        } else if ty > tz {
            Axis::Y
        } else {
            Axis::Z
        }
    }

    /// A unit vector orthogonal to `self`, chosen deterministically (built from
    /// the smallest-magnitude axis for numerical stability), as `r3.Vector.Ortho`.
    #[must_use]
    pub fn ortho(self) -> Self {
        let (tx, ty, tz) = (self.x.abs(), self.y.abs(), self.z.abs());
        // pick the axis with the smallest component
        let other = if tx < ty && tx < tz {
            Self::new(1.0, 0.0, 0.0)
        } else if ty < tz {
            Self::new(0.0, 1.0, 0.0)
        } else {
            Self::new(0.0, 0.0, 1.0)
        };
        self.cross(other).normalize()
    }
}
