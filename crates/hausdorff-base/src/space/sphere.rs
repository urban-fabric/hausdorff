//! `sphere` — the **unit hypersphere `Sⁿ⁻¹`** as a space, plus the [`Sphere`] manifold
//! trait that the concrete spheres a distribution can be defined over all implement.
//!
//! A probability distribution *over a sphere* (von Mises–Fisher and kin) needs the
//! sphere as an object: its ambient dimension, the inner product (cosine of the geodesic
//! angle), and the boundary to and from ambient ℝⁿ coordinates. [`Sphere`] captures
//! exactly that, implemented by the fixed `S2Point` (in `hausdorff-sphere`) for the
//! 2-sphere and by the dimension-general [`UnitVector`] here for `Sⁿ⁻¹` — so the
//! distribution is written once and runs on the geographic 2-sphere *and* on the
//! high-dimensional hypersphere that L2-normalized embeddings live on.

/// A point on a unit sphere, viewed as the **manifold a distribution lives on**. The
/// minimal surface a directional distribution needs: the ambient dimension `n`, the inner
/// product, and conversion to/from ambient ℝⁿ coordinates (the rest — densities,
/// sampling, fitting — is generic over this trait).
pub trait Sphere: Clone {
    /// The ambient dimension `n` (the sphere is `Sⁿ⁻¹`).
    fn ambient_dim(&self) -> usize;

    /// The inner product `self · other` — the cosine of the geodesic (central) angle.
    fn dot(&self, other: &Self) -> f64;

    /// Write this point's ambient **unit** coordinates into `out` (length `ambient_dim`).
    fn write_ambient(&self, out: &mut [f64]);

    /// The point in the direction of the ambient vector `coords` (its normalization onto
    /// the sphere). `coords` is the space's ambient dimension and must be non-zero.
    fn from_ambient(coords: &[f64]) -> Self;
}

/// A point on the unit hypersphere `Sⁿ⁻¹` — a unit vector in ℝⁿ, dimension-general
/// (`n = ambient_dim`, fixed by construction). The space high-dimensional directional
/// data lives on (e.g. L2-normalized embeddings), the general counterpart of the fixed
/// `S2Point`.
#[derive(Clone, Debug, PartialEq)]
pub struct UnitVector(Box<[f64]>);

impl UnitVector {
    /// The unit vector in the **direction of** `coords` (its normalization). `coords` must
    /// be non-zero (a zero direction yields NaNs, as for any normalization).
    #[must_use]
    pub fn from_coords(coords: &[f64]) -> Self {
        let norm = coords.iter().map(|c| c * c).sum::<f64>().sqrt();
        Self(coords.iter().map(|c| c / norm).collect())
    }

    /// The ambient unit coordinates.
    #[must_use]
    pub fn coords(&self) -> &[f64] {
        &self.0
    }

    /// The ambient dimension `n`.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.0.len()
    }

    /// The inner product `self · other` (cosine of the geodesic angle).
    #[must_use]
    pub fn dot(&self, other: &Self) -> f64 {
        self.0.iter().zip(other.0.iter()).map(|(a, b)| a * b).sum()
    }
}

impl Sphere for UnitVector {
    fn ambient_dim(&self) -> usize {
        self.0.len()
    }
    fn dot(&self, other: &Self) -> f64 {
        UnitVector::dot(self, other)
    }
    fn write_ambient(&self, out: &mut [f64]) {
        out.copy_from_slice(&self.0);
    }
    fn from_ambient(coords: &[f64]) -> Self {
        Self::from_coords(coords)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_coords_normalizes() {
        let u = UnitVector::from_coords(&[3.0, 0.0, 4.0, 0.0]);
        assert_eq!(u.dim(), 4);
        assert!((u.dot(&u) - 1.0).abs() < 1e-12);
        assert!((u.coords()[0] - 0.6).abs() < 1e-12);
        assert!((u.coords()[2] - 0.8).abs() < 1e-12);
    }

    #[test]
    fn dot_is_cosine_and_orthogonality_is_zero() {
        let e0 = UnitVector::from_coords(&[1.0, 0.0, 0.0, 0.0]);
        let e1 = UnitVector::from_coords(&[0.0, 1.0, 0.0, 0.0]);
        assert!((e0.dot(&e0) - 1.0).abs() < 1e-12);
        assert!(e0.dot(&e1).abs() < 1e-12);
    }

    #[test]
    fn sphere_trait_round_trips_ambient_coordinates() {
        let u = UnitVector::from_coords(&[1.0, 2.0, -2.0]);
        let mut out = vec![0.0; u.ambient_dim()];
        u.write_ambient(&mut out);
        let v = <UnitVector as Sphere>::from_ambient(&out);
        assert!((u.dot(&v) - 1.0).abs() < 1e-12);
    }
}
