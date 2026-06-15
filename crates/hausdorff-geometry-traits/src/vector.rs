//! Vector-space algebra shared by the coordinate types, after `spire`'s
//! `VectorSpace` / `InnerProductSpace` / `NormedVectorSpace` hierarchy.
//!
//! `Coord` (the plane, in `hausdorff-planar`) and `Vector` (the sphere's ℝ³
//! substrate, in `hausdorff-base`) already carry the operations inherently;
//! these traits give them a shared abstraction so generic code (centroids,
//! orthogonalization, generic norms) can be written once. The scalar field is
//! fixed to `f64`.

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
