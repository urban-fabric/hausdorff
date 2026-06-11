//! Vector-space algebra shared by the coordinate types, after `spire`'s
//! `VectorSpace` / `InnerProductSpace` / `NormedVectorSpace` hierarchy.
//!
//! [`crate::space::r2::Coord`] (the plane) and [`hausdorff_base::space::r3::Vector`] (the sphere's ℝ³
//! substrate) already carry the operations inherently; these traits give them a
//! shared abstraction so generic code (centroids, orthogonalization, generic
//! norms) can be written once. The scalar field is fixed to `f64`.

use crate::space::r2::Coord;
use hausdorff_base::space::r3::Vector;

/// An `f64`-scalar vector space: addition and scaling with an identity.
pub trait VectorSpace: Copy {
    /// The additive identity (origin / zero vector).
    fn zero() -> Self;
    /// Vector addition.
    fn add(self, rhs: Self) -> Self;
    /// Vector subtraction.
    fn sub(self, rhs: Self) -> Self;
    /// Scalar multiple.
    fn scale(self, factor: f64) -> Self;
}

/// A vector space equipped with an inner (dot) product.
pub trait InnerProduct: VectorSpace {
    /// The inner (dot) product.
    fn dot(self, rhs: Self) -> f64;
}

/// An inner-product space with the induced Euclidean norm.
pub trait Norm: InnerProduct {
    /// Squared norm (`self · self`) — avoids the `sqrt` when only comparing.
    fn norm_squared(self) -> f64 {
        self.dot(self)
    }
    /// Euclidean norm.
    fn norm(self) -> f64 {
        self.norm_squared().sqrt()
    }
    /// A unit vector in the same direction, or [`zero`](VectorSpace::zero) if
    /// `self` is the zero vector.
    fn normalized(self) -> Self {
        let n2 = self.norm_squared();
        if n2 == 0.0 {
            Self::zero()
        } else {
            self.scale(1.0 / n2.sqrt())
        }
    }
}

impl VectorSpace for Vector {
    fn zero() -> Self {
        Vector::new(0.0, 0.0, 0.0)
    }
    fn add(self, rhs: Self) -> Self {
        Vector::add(self, rhs)
    }
    fn sub(self, rhs: Self) -> Self {
        Vector::sub(self, rhs)
    }
    fn scale(self, factor: f64) -> Self {
        self.mul(factor)
    }
}

impl InnerProduct for Vector {
    fn dot(self, rhs: Self) -> f64 {
        Vector::dot(self, rhs)
    }
}

impl Norm for Vector {}

impl VectorSpace for Coord {
    fn zero() -> Self {
        Coord::new(0.0, 0.0)
    }
    fn add(self, rhs: Self) -> Self {
        self + rhs
    }
    fn sub(self, rhs: Self) -> Self {
        self - rhs
    }
    fn scale(self, factor: f64) -> Self {
        self * factor
    }
}

impl InnerProduct for Coord {
    fn dot(self, rhs: Self) -> f64 {
        Coord::dot(self, rhs)
    }
}

impl Norm for Coord {}
