//! Vector-space algebra shared by the coordinate types, after `spire`'s
//! `VectorSpace` / `InnerProductSpace` / `NormedVectorSpace` hierarchy.
//!
//! `Coord` (the plane, in `hausdorff-planar`) and `Vector` (the sphere's ℝ³
//! substrate, in `hausdorff-base`) already carry the operations inherently;
//! these traits give them a shared abstraction so generic code (centroids,
//! orthogonalization, generic norms) can be written once. The scalar field is the
//! associated [`Scalar`](VectorSpace::Scalar): any [`LazyReal`] (`f64` by default,
//! `f32`, a lazy tensor value, …).

use euclid::numeric::LazyReal;

/// A vector space over a [`LazyReal`] scalar: addition and scaling with an identity.
pub trait VectorSpace: Clone {
    /// The scalar field this vector space is over.
    type Scalar: LazyReal;
    /// The additive identity (origin / zero vector).
    fn zero() -> Self;
    /// Vector addition.
    fn add(self, rhs: Self) -> Self;
    /// Vector subtraction.
    fn sub(self, rhs: Self) -> Self;
    /// Scalar multiple.
    fn scale(self, factor: Self::Scalar) -> Self;
}

/// A vector space equipped with an inner (dot) product.
pub trait InnerProduct: VectorSpace {
    /// The inner (dot) product.
    fn dot(self, rhs: Self) -> Self::Scalar;
}

/// An inner-product space with the induced Euclidean norm.
pub trait Norm: InnerProduct {
    /// Squared norm (`self · self`) — avoids the `sqrt` when only comparing.
    fn norm_squared(self) -> Self::Scalar {
        self.clone().dot(self)
    }
    /// Euclidean norm.
    fn norm(self) -> Self::Scalar {
        self.norm_squared().sqrt()
    }
    /// A unit vector in the same direction, or [`zero`](VectorSpace::zero) if
    /// `self` is the zero vector.
    fn normalized(self) -> Self
    where
        Self::Scalar: PartialEq,
    {
        let n2 = self.clone().norm_squared();
        if n2 == Self::Scalar::zero() {
            Self::zero()
        } else {
            self.scale(n2.sqrt().recip())
        }
    }
}
