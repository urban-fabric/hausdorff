//! The ellipsoid **latitude system** and surface-of-revolution radii — the
//! conversions that make this a *geodetic* (oblate-spheroid) kernel rather than
//! a spherical one.
//!
//! On a sphere there is one latitude; on an oblate spheroid the single notion
//! splits into three, and ECEF / geodesics / projections each need a specific
//! one:
//!
//! - **geodetic** `φ` — the angle the surface *normal* makes with the equator
//!   (what GPS reports);
//! - **geocentric** `ψ` — the angle of the *position vector*, `tan ψ = (1−e²)tan φ`;
//! - **reduced / parametric** `β` — the latitude on the geodesic *auxiliary
//!   sphere* (Karney/Vincenty) and the variable Bowring's ECEF inverse iterates,
//!   `tan β = (b/a) tan φ = (1−f) tan φ`.
//!
//! These were previously re-derived inline — the `N` in [`Ellipsoid::geodetic_to_ecef`],
//! Bowring's `(1−f)·tanφ` in [`Ellipsoid::ecef_to_geodetic`], the `up` normal in
//! [`Ellipsoid::enu_frame`]. This module makes them first-class and reusable, and
//! adds the meridian radius of curvature `M` (previously absent).
//!
//! Every conversion here is machine-checked in `riemann`
//! (`RiemannLean/Hausdorff/{Geodetic,Spheroid}.lean`):
//!
//! - geodetic→reduced `tan β = (b/a) tan φ` — `reduced_tan_of_geodetic`;
//! - geodetic→geocentric `tan ψ = (1−e²) tan φ` — `geocentric_tan_of_geodetic`;
//! - the surface normal `(cosφcosλ, cosφsinλ, sinφ)` *is* the geodetic up
//!   (so `φ` is the normal angle — the defining property) — `ecef_normal_is_geodetic`,
//!   `Spheroid.surface_normal`;
//! - the spheroid is a surface of revolution whose first fundamental form is
//!   diagonal (meridians ⟂ parallels), with `G = (N cosφ)²` and `E = M²` — the
//!   prime-vertical and meridian radii — `Spheroid.{parallel_normSq, meridian_normSq,
//!   meridian_perp_parallel}`.

use crate::earth::ellipsoid::Ellipsoid;
use crate::space::surface::SurfaceOfRevolution;

/// The oblate spheroid is a [`SurfaceOfRevolution`]: its meridian arc speed is the
/// meridian radius of curvature `M = √E` and its parallel radius is `N cosφ = √G`
/// (`riemann` `Hausdorff/Spheroid`). The sphere is `Ellipsoid::SPHERE`.
impl SurfaceOfRevolution for Ellipsoid {
    fn meridian_speed(&self, lat_deg: f64) -> f64 {
        self.meridian_radius(lat_deg)
    }

    fn parallel_radius(&self, lat_deg: f64) -> f64 {
        self.prime_vertical_radius(lat_deg) * lat_deg.to_radians().cos()
    }

    fn meridian_domain_deg(&self) -> (f64, f64) {
        (-90.0, 90.0)
    }
}

impl Ellipsoid {
    /// **Reduced (parametric) latitude** `β` of the geodetic latitude `lat_deg`:
    /// `tan β = (b/a) tan φ = (1−f) tan φ`. This is the latitude on the geodesic
    /// *auxiliary sphere* (Karney/Vincenty) and the value Bowring's
    /// [`Ellipsoid::ecef_to_geodetic`] iterates. On a sphere `β = φ`.
    ///
    /// Proof: `riemann` `Hausdorff/Geodetic.reduced_tan_of_geodetic`.
    #[must_use]
    pub fn reduced_latitude(self, lat_deg: f64) -> f64 {
        let (s, c) = lat_deg.to_radians().sin_cos();
        ((1.0 - self.f) * s).atan2(c).to_degrees()
    }

    /// Geodetic latitude `φ` from the reduced latitude `β`: `tan φ = (a/b) tan β`
    /// — the inverse of [`Ellipsoid::reduced_latitude`].
    #[must_use]
    pub fn geodetic_from_reduced(self, beta_deg: f64) -> f64 {
        let (s, c) = beta_deg.to_radians().sin_cos();
        (s / (1.0 - self.f)).atan2(c).to_degrees()
    }

    /// **Geocentric latitude** `ψ` of the geodetic latitude `lat_deg`: the angle
    /// of the position vector from the centre, `tan ψ = (1−e²) tan φ`. The
    /// spheroid squashes it toward the equator (`|ψ| ≤ |φ|`); on a sphere `ψ = φ`.
    ///
    /// Proof: `riemann` `Hausdorff/Geodetic.geocentric_tan_of_geodetic`.
    #[must_use]
    pub fn geocentric_latitude(self, lat_deg: f64) -> f64 {
        let (s, c) = lat_deg.to_radians().sin_cos();
        ((1.0 - self.e2()) * s).atan2(c).to_degrees()
    }

    /// Geodetic latitude `φ` from the geocentric latitude `ψ`:
    /// `tan φ = tan ψ / (1−e²)` — the inverse of [`Ellipsoid::geocentric_latitude`].
    #[must_use]
    pub fn geodetic_from_geocentric(self, psi_deg: f64) -> f64 {
        let (s, c) = psi_deg.to_radians().sin_cos();
        (s / (1.0 - self.e2())).atan2(c).to_degrees()
    }

    /// **Prime-vertical radius of curvature** `N = a / √(1 − e² sin²φ)`: the
    /// length of the surface normal from the point to the polar axis. The `N` of
    /// geodetic→ECEF (`X = N cosφ`); the circle of latitude has radius `N cosφ`,
    /// so `G = (N cosφ)²` is the parallel coefficient of the first fundamental
    /// form. `N ≥ a`, with equality at the equator.
    ///
    /// Proof: `N²(1 − e² sin²φ) = a²` underlies
    /// `Hausdorff/Ellipsoid.ecef_on_ellipsoid` and
    /// `Hausdorff/Geodetic.ecef_normal_is_geodetic`; `G = (N cosφ)²` is
    /// `Hausdorff/Spheroid.parallel_normSq`.
    #[must_use]
    pub fn prime_vertical_radius(self, lat_deg: f64) -> f64 {
        let s = lat_deg.to_radians().sin();
        self.a / (1.0 - self.e2() * s * s).sqrt()
    }

    /// **Meridian radius of curvature** `M = a(1 − e²) / (1 − e² sin²φ)^{3/2}`:
    /// the radius of curvature in the north–south (meridian) direction, and the
    /// meridian arc-speed `E = M²` (the other diagonal coefficient of the first
    /// fundamental form). With `N` these are the two principal radii of the
    /// spheroid; on a sphere both collapse to `a`.
    ///
    /// Structure: `Hausdorff/Spheroid.{meridian_normSq, parallel_normSq}` — the
    /// first fundamental form is diagonal (`Hausdorff/Spheroid.meridian_perp_parallel`:
    /// meridians ⟂ parallels), so `M`, `N cosφ` fully describe the local metric.
    #[must_use]
    pub fn meridian_radius(self, lat_deg: f64) -> f64 {
        let s = lat_deg.to_radians().sin();
        let w2 = 1.0 - self.e2() * s * s;
        self.a * (1.0 - self.e2()) / (w2 * w2.sqrt())
    }

    /// The **surface-area element** `dA / (dφ dλ) = M · N cosφ` (m² per radian²)
    /// — the integrand of the spheroid's surface area `∬ M N cosφ dφ dλ`. It is
    /// `‖rφ × rλ‖`, the area element of a parametrized surface; the diagonal
    /// first fundamental form makes it `√(E·G) = M · (N cosφ)` (the two principal
    /// radii). Integrating it over the whole spheroid reproduces
    /// [`Ellipsoid::total_area_m2`].
    ///
    /// Proof: `riemann` `Hausdorff/Spheroid.area_element_eq`
    /// (`‖rφ × rλ‖² = E·G`, the Gram determinant of the first fundamental form).
    #[must_use]
    pub fn area_element(self, lat_deg: f64) -> f64 {
        self.meridian_radius(lat_deg) * self.prime_vertical_radius(lat_deg)
            * lat_deg.to_radians().cos()
    }
}

/// The outward **unit surface normal** ("geodetic up") at `(lat_deg, lon_deg)`:
/// `(cosφ cosλ, cosφ sinλ, sinφ)`. The angle it makes with the equatorial plane
/// *is* the geodetic latitude — the defining property of geodetic latitude — and
/// it is the `up` axis of the local ENU frame ([`Ellipsoid::enu_frame`]).
///
/// Notably this is **independent of the datum** (`a`, `f`): every reference
/// ellipsoid shares the same geodetic-up direction at a given `(φ, λ)`, because
/// the normal direction is `(cosφ, sinφ)` in the meridian plane regardless of
/// the axes — the content of `riemann` `Hausdorff/Geodetic.ecef_normal_is_geodetic`
/// (the normal is `(N/a²)·(cosφ, sinφ)`) and `Hausdorff/Spheroid.surface_normal`.
#[must_use]
pub fn geodetic_normal(lat_deg: f64, lon_deg: f64) -> [f64; 3] {
    let (slat, clat) = lat_deg.to_radians().sin_cos();
    let (slon, clon) = lon_deg.to_radians().sin_cos();
    [clat * clon, clat * slon, slat]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    /// On a sphere (`f = 0`) the three latitudes coincide and both radii are `a`.
    #[test]
    fn sphere_collapses_the_latitude_system() {
        let s = Ellipsoid::SPHERE;
        for &lat in &[0.0, 12.3, 45.0, -60.0, 89.9] {
            assert!(close(s.reduced_latitude(lat), lat, 1e-12));
            assert!(close(s.geocentric_latitude(lat), lat, 1e-12));
            assert!(close(s.prime_vertical_radius(lat), s.a, 1e-6));
            assert!(close(s.meridian_radius(lat), s.a, 1e-6));
        }
    }

    /// Each conversion inverts cleanly (WGS84).
    #[test]
    fn conversions_round_trip() {
        let e = Ellipsoid::WGS84;
        for &lat in &[0.0, 12.3, 45.0, -60.0, 80.0] {
            assert!(close(e.geodetic_from_reduced(e.reduced_latitude(lat)), lat, 1e-9));
            assert!(close(e.geodetic_from_geocentric(e.geocentric_latitude(lat)), lat, 1e-9));
        }
    }

    /// For `0 < φ < 90` the spheroid squashes the latitudes toward the equator,
    /// `|ψ| < |β| < |φ|`, with the standard WGS84 values at 45°.
    #[test]
    fn latitudes_ordered_and_known() {
        let e = Ellipsoid::WGS84;
        let beta = e.reduced_latitude(45.0);
        let psi = e.geocentric_latitude(45.0);
        assert!(psi < beta && beta < 45.0);
        assert!(close(beta, 44.9038, 1e-2));
        assert!(close(psi, 44.8076, 1e-2));
    }

    /// The geodetic normal is a unit vector whose elevation (z-component) is
    /// `sin φ` — i.e. its angle above the equatorial plane is the geodetic
    /// latitude, the defining property.
    #[test]
    fn normal_is_unit_with_latitude_elevation() {
        for &(lat, lon) in &[(0.0, 0.0), (45.0, 30.0), (-60.0, 170.0), (89.0, -120.0)] {
            let n = geodetic_normal(lat, lon);
            let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(close(norm, 1.0, 1e-12));
            assert!(close(n[2], lat.to_radians().sin(), 1e-12));
        }
    }

    /// `N ≥ a` (equality at the equator); `M` brackets `a` (smaller at the
    /// equator, larger at the pole).
    #[test]
    fn radii_bracket_the_semi_major_axis() {
        let e = Ellipsoid::WGS84;
        assert!(close(e.prime_vertical_radius(0.0), e.a, 1e-6));
        assert!(e.prime_vertical_radius(45.0) > e.a);
        assert!(e.meridian_radius(0.0) < e.a);
        assert!(e.meridian_radius(90.0) > e.a);
    }

    /// The generic [`SurfaceOfRevolution::total_area`] over the ellipsoid (the
    /// trait's numeric integrator) reproduces the closed-form `total_area_m2`, and
    /// the sphere instance reproduces `4πR²`.
    #[test]
    fn surface_of_revolution_total_area() {
        use crate::space::surface::SurfaceOfRevolution;
        let e = Ellipsoid::WGS84;
        let rel = (e.total_area(200_000) - e.total_area_m2()).abs() / e.total_area_m2();
        assert!(rel < 1e-7, "WGS84 SoR area rel_err = {rel}");
        // F = 0 (orthogonal net); the trait's generic area element (√(EG)) agrees
        // with the inherent closed form M·N·cosφ.
        let (_, f, _) = e.first_fundamental_form(45.0);
        assert_eq!(f, 0.0);
        let generic = SurfaceOfRevolution::area_element(&e, 45.0);
        assert!((generic - e.area_element(45.0)).abs() / generic < 1e-12);
        // Sphere: total area 4πR².
        let s = Ellipsoid::SPHERE;
        let exact = 4.0 * core::f64::consts::PI * s.a * s.a;
        assert!((s.total_area(200_000) - exact).abs() / exact < 1e-7);
    }

    /// The local area element `M·N·cosφ`, integrated over the whole spheroid,
    /// reproduces the closed-form `total_area_m2` — two independent area
    /// computations (the first-fundamental-form integrand vs. the authalic
    /// radius) agree.
    #[test]
    fn area_element_integrates_to_total_area() {
        use core::f64::consts::{FRAC_PI_2, PI};
        let e = Ellipsoid::WGS84;
        let n = 100_000;
        let dphi = PI / n as f64; // φ from -π/2 to π/2 (radians)
        let mut sum = 0.0;
        for i in 0..n {
            let phi = -FRAC_PI_2 + (i as f64 + 0.5) * dphi; // midpoint
            sum += e.area_element(phi.to_degrees()) * dphi;
        }
        let area = 2.0 * PI * sum; // ∫ dλ = 2π
        let rel_err = (area - e.total_area_m2()).abs() / e.total_area_m2();
        assert!(rel_err < 1e-7, "rel_err = {rel_err}");
    }
}
