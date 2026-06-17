//! Planar [`VectorSpace`] / [`InnerProduct`] / [`Norm`] instance for [`Coord`].
//! The trait hierarchy lives upstream in [`hausdorff_geometry_traits`]; this just
//! routes it onto `Coord`'s inherent arithmetic so generic code (centroids,
//! orthogonalization, generic norms) runs in the plane.

use hausdorff_geometry_traits::{InnerProduct, Norm, VectorSpace};

use crate::space::r2::Coord;

impl VectorSpace for Coord {
    type Scalar = f64;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coord_vector_space_and_norm() {
        let v = Coord::new(3.0, 4.0);
        assert_eq!(<Coord as VectorSpace>::zero(), Coord::new(0.0, 0.0));
        assert_eq!(VectorSpace::add(v, Coord::new(1.0, 1.0)), Coord::new(4.0, 5.0));
        assert!((Norm::norm(v) - 5.0).abs() < 1e-12);
        assert!((Norm::norm_squared(v) - 25.0).abs() < 1e-12);
        // `normalized` is a unit vector; the zero vector maps to zero.
        assert!((Norm::norm(Norm::normalized(v)) - 1.0).abs() < 1e-12);
        assert_eq!(Norm::normalized(Coord::new(0.0, 0.0)), Coord::new(0.0, 0.0));
    }
}
