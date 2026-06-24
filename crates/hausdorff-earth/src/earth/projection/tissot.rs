//! [`Tissot`]'s indicatrix — the local distortion ellipse of a map projection.
//!
//! Map distortion *is* an ellipse: an infinitesimal circle on the ground projects
//! to an ellipse, and its shape is exactly the singular-value decomposition of the
//! projection's 2×2 Jacobian `J`. So this reuses the symmetric-2×2 eigenstructure
//! already proven for [`ConfidenceEllipse`] (`riemann`
//! `Hausdorff/Ellipse.{sym2_eigenvalue, sym2_eig_prod_eq_det}`) — the indicatrix
//! is the confidence ellipse of the Gram matrix `JᵀJ`.

use hausdorff_planar::space::r2::{ConfidenceEllipse, Coord};

/// **Tissot's indicatrix** — the distortion ellipse an infinitesimal ground circle
/// maps to under a projection, the standard summary of local map distortion.
///
/// It is built from the projection's 2×2 Jacobian `J` measured against true ground
/// distance (projected metres per unit east/north ground metre, so an undistorted
/// map has `J = I` and a unit-circle indicatrix). The indicatrix semi-axes are the
/// **singular values** of `J` — `√` the eigenvalues of the symmetric Gram matrix
/// `JᵀJ` — so it is exactly the [`ConfidenceEllipse`] of `JᵀJ`, reusing that
/// proven closed form (`riemann` `Hausdorff/Ellipse.sym2_eigenvalue`).
///
/// From it: the principal scale factors (max/min linear scale), the **area scale**
/// `a·b = |det J|`, and the **maximum angular distortion** `ω`. A projection is
/// **conformal** at the point iff the indicatrix is a circle (`a = b`) and
/// **equal-area** iff the area scale is `1`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tissot {
    indicatrix: ConfidenceEllipse,
}

impl Tissot {
    /// The indicatrix of the projection Jacobian
    /// `j = [[∂x/∂E, ∂x/∂N], [∂y/∂E, ∂y/∂N]]` (projected metres per unit ground
    /// distance east/north). The Gram matrix `JᵀJ` is symmetric positive-(semi)definite;
    /// the indicatrix is its [`ConfidenceEllipse`] (`k = 1`), whose semi-axes are
    /// the singular values of `J`.
    #[must_use]
    pub fn from_jacobian(j: [[f64; 2]; 2]) -> Self {
        let [[a, b], [c, d]] = j;
        // JᵀJ — the metric by which the unit ground circle is measured.
        let m11 = a * a + c * c;
        let m12 = a * b + c * d;
        let m22 = b * b + d * d;
        Self {
            indicatrix: ConfidenceEllipse::from_covariance(
                Coord::new(0.0, 0.0),
                m11,
                m12,
                m22,
                1.0,
            ),
        }
    }

    /// The principal scale factors `(a, b)` — the singular values of `J`, the
    /// **maximum and minimum linear scale** at the point (the indicatrix
    /// semi-axes). For a conformal projection `a = b`.
    #[must_use]
    pub fn principal_scales(self) -> (f64, f64) {
        (self.indicatrix.semi_major(), self.indicatrix.semi_minor())
    }

    /// The **area scale** `a·b = |det J|` (projected area ÷ true area) — equal to
    /// `√det(JᵀJ)`, the product of the indicatrix eigenvalues (`riemann`
    /// `Hausdorff/Ellipse.sym2_eig_prod_eq_det`). `1` for an equal-area map.
    #[must_use]
    pub fn area_scale(self) -> f64 {
        self.indicatrix.det().max(0.0).sqrt()
    }

    /// The **maximum angular distortion** `ω` (radians): `sin(ω/2) = (a − b)/(a + b)`,
    /// the largest error a projected angle can suffer. Zero iff conformal.
    #[must_use]
    pub fn max_angular_distortion(self) -> f64 {
        let (a, b) = self.principal_scales();
        if a + b == 0.0 {
            return 0.0;
        }
        2.0 * ((a - b) / (a + b)).asin()
    }

    /// True iff the projection is **conformal** at the point — the indicatrix is a
    /// circle (`a = b`) to within the relative tolerance `tol`.
    #[must_use]
    pub fn is_conformal(self, tol: f64) -> bool {
        let (a, b) = self.principal_scales();
        (a - b).abs() <= tol * a.max(b)
    }

    /// True iff the projection is **equal-area** at the point — the area scale is
    /// `1` to within the relative tolerance `tol`.
    #[must_use]
    pub fn is_equal_area(self, tol: f64) -> bool {
        (self.area_scale() - 1.0).abs() <= tol
    }

    /// The distortion ellipse itself — semi-axes and orientation, ready to draw
    /// (its [`ConfidenceEllipse::axes`] are the principal directions).
    #[must_use]
    pub fn indicatrix(self) -> ConfidenceEllipse {
        self.indicatrix
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    /// An isotropic Jacobian `k·I` (a conformal projection's local map, e.g. Web
    /// Mercator's `sec φ` scale) gives a *circular* indicatrix: equal principal
    /// scales `k`, area scale `k²`, zero angular distortion, conformal.
    #[test]
    fn isotropic_jacobian_is_conformal_circle() {
        let k = 2.0; // e.g. Web Mercator point scale at 60°N (sec 60° = 2).
        let t = Tissot::from_jacobian([[k, 0.0], [0.0, k]]);
        let (a, b) = t.principal_scales();
        assert!(close(a, k, 1e-12) && close(b, k, 1e-12));
        assert!(close(t.area_scale(), k * k, 1e-12));
        assert!(close(t.max_angular_distortion(), 0.0, 1e-12));
        assert!(t.is_conformal(1e-9));
        assert!(!t.is_equal_area(1e-9));
    }

    /// A rotation composed with a uniform scale is still conformal — the
    /// indicatrix only sees `JᵀJ`, which is rotation-invariant.
    #[test]
    fn rotation_plus_scale_stays_conformal() {
        let (k, theta) = (1.7_f64, 0.6_f64);
        let (s, c) = theta.sin_cos();
        let t = Tissot::from_jacobian([[k * c, -k * s], [k * s, k * c]]);
        let (a, b) = t.principal_scales();
        assert!(close(a, k, 1e-12) && close(b, k, 1e-12));
        assert!(t.is_conformal(1e-9));
    }

    /// An anisotropic Jacobian distorts: different principal scales, nonzero
    /// angular distortion, not conformal. With `det J = 1` it is still equal-area.
    #[test]
    fn anisotropic_jacobian_distorts_but_can_preserve_area() {
        // J = diag(2, 1/2): area scale = 1 (equal-area), but a ≠ b (not conformal).
        let t = Tissot::from_jacobian([[2.0, 0.0], [0.0, 0.5]]);
        let (a, b) = t.principal_scales();
        assert!(close(a, 2.0, 1e-12) && close(b, 0.5, 1e-12));
        assert!(close(t.area_scale(), 1.0, 1e-12));
        assert!(t.is_equal_area(1e-9));
        assert!(!t.is_conformal(1e-9));
        // ω: sin(ω/2) = (2 − 0.5)/(2 + 0.5) = 0.6 ⇒ ω = 2·asin(0.6).
        assert!(close(
            t.max_angular_distortion(),
            2.0 * 0.6_f64.asin(),
            1e-12
        ));
    }
}
