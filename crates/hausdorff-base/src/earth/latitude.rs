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
//! The **map projections** add three more *auxiliary latitudes* — transcendental
//! functions of `φ` (logs / a meridian-arc integral, outside the algebraic proof
//! boundary, like the Karney series) that each linearize one projection property:
//!
//! - **conformal** `χ` — the sphere latitude on which a conformal map (Mercator,
//!   Transverse Mercator, the conformal conics) is angle-preserving;
//! - **authalic** `ξ` — the equal-area latitude: a graticule cell on the authalic
//!   sphere has the spheroid's true area (Albers, ellipsoidal polygon area);
//! - **rectifying** `μ` — the equidistant latitude, *proportional to meridian-arc
//!   distance* (`μ = (π/2)·s(φ)/s_quarter`), the basis of the Transverse Mercator
//!   northing and equidistant projections.
//!
//! All six collapse to `φ` on a sphere (`f = 0`). The geodetic / geocentric /
//! reduced relations were previously re-derived inline — the `N` in [`Ellipsoid::geodetic_to_ecef`],
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

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

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
        self.meridian_radius(lat_deg)
            * self.prime_vertical_radius(lat_deg)
            * lat_deg.to_radians().cos()
    }

    /// **Conformal latitude** `χ` of the geodetic latitude `lat_deg`: the latitude
    /// on the conformal sphere whose Mercator / Transverse-Mercator / conformal-conic
    /// image is angle-preserving. `χ = π/2 − 2·arctan(t)` for the isometric building
    /// block `t` (Snyder 3-1 / 15-9); `|χ| ≤ |φ|` and `χ = φ` on a sphere.
    ///
    /// Transcendental (a log term in `e`), so outside the algebraic proof boundary;
    /// its `e → 0` limit `χ = φ` and the conformality condition are the provable
    /// parts (cf. `riemann` `Hausdorff/{Mercator, Stereographic}`).
    #[must_use]
    pub fn conformal_latitude(self, lat_deg: f64) -> f64 {
        if self.f == 0.0 {
            return lat_deg;
        }
        let e = self.e2().sqrt();
        (FRAC_PI_2 - 2.0 * conformal_t(lat_deg.to_radians(), e).atan()).to_degrees()
    }

    /// Geodetic latitude `φ` from the conformal latitude `χ` — the inverse of
    /// [`conformal_latitude`](Self::conformal_latitude), by Snyder's fixed-point
    /// iteration (a few rounds for any terrestrial latitude).
    #[must_use]
    pub fn geodetic_from_conformal(self, chi_deg: f64) -> f64 {
        if self.f == 0.0 {
            return chi_deg;
        }
        let e = self.e2().sqrt();
        let t = (FRAC_PI_4 - chi_deg.to_radians() / 2.0).tan();
        lat_from_conformal_t(t, e).to_degrees()
    }

    /// **Authalic latitude** `ξ` of the geodetic latitude `lat_deg`: the latitude
    /// on the equal-area (authalic) sphere of radius `Rq` for which a graticule
    /// cell has the spheroid's true area — `sin ξ = q(φ) / q(π/2)` (Snyder 3-11).
    /// The basis of equal-area maps (Albers) and exact ellipsoidal polygon area.
    /// `|ξ| ≤ |φ|`; `ξ = φ` on a sphere.
    ///
    /// Transcendental (`q` carries a log term); the area-preservation property and
    /// the `e → 0` limit `ξ = φ` are the provable parts.
    #[must_use]
    pub fn authalic_latitude(self, lat_deg: f64) -> f64 {
        if self.f == 0.0 {
            return lat_deg;
        }
        let e = self.e2().sqrt();
        (authalic_q(lat_deg.to_radians(), e) / authalic_q(FRAC_PI_2, e))
            .asin()
            .to_degrees()
    }

    /// Geodetic latitude `φ` from the authalic latitude `ξ` — the inverse of
    /// [`authalic_latitude`](Self::authalic_latitude), by Snyder's iteration (3-18).
    #[must_use]
    pub fn geodetic_from_authalic(self, xi_deg: f64) -> f64 {
        if self.f == 0.0 {
            return xi_deg;
        }
        let e = self.e2().sqrt();
        let q = xi_deg.to_radians().sin() * authalic_q(FRAC_PI_2, e);
        lat_from_authalic_q(q, e).to_degrees()
    }

    /// **Meridian arc length** from the equator to `lat_deg`, metres — the
    /// north–south distance along the ellipsoid surface (Snyder 3-21, a series in
    /// `e²`). Signed (negative south of the equator). The quarter-meridian
    /// `meridian_arc_m(90°)` is the pole distance.
    #[must_use]
    pub fn meridian_arc_m(self, lat_deg: f64) -> f64 {
        meridian_arc(lat_deg.to_radians(), self.a, self.e2())
    }

    /// **Rectifying latitude** `μ` of the geodetic latitude `lat_deg`: the
    /// equidistant latitude *proportional to meridian-arc distance*,
    /// `μ = (π/2)·s(φ)/s(90°)` where `s` is [`meridian_arc_m`](Self::meridian_arc_m).
    /// Equal meridian-arc steps are equal `μ` steps — the basis of the Transverse
    /// Mercator northing and equidistant projections. `μ = φ` on a sphere.
    #[must_use]
    pub fn rectifying_latitude(self, lat_deg: f64) -> f64 {
        let quarter = self.meridian_arc_m(90.0);
        if quarter == 0.0 {
            return lat_deg;
        }
        90.0 * self.meridian_arc_m(lat_deg) / quarter
    }

    /// Geodetic latitude `φ` from the rectifying latitude `μ` — the inverse of
    /// [`rectifying_latitude`](Self::rectifying_latitude), by Snyder's *footpoint
    /// latitude* series (3-26) in the third flattening `e₁`.
    #[must_use]
    pub fn geodetic_from_rectifying(self, mu_deg: f64) -> f64 {
        let s = (1.0 - self.e2()).sqrt();
        let e1 = (1.0 - s) / (1.0 + s);
        let mu = mu_deg.to_radians();
        (mu + (1.5 * e1 - 27.0 * e1.powi(3) / 32.0) * (2.0 * mu).sin()
            + (21.0 * e1 * e1 / 16.0 - 55.0 * e1.powi(4) / 32.0) * (4.0 * mu).sin()
            + (151.0 * e1.powi(3) / 96.0) * (6.0 * mu).sin()
            + (1097.0 * e1.powi(4) / 512.0) * (8.0 * mu).sin())
        .to_degrees()
    }
}

// ---- Auxiliary-latitude kernels (Snyder PP1395), shared with the map
// projections in `hausdorff-earth`. Radians in/out; `e` is the first
// eccentricity. The `Ellipsoid` methods above wrap these in the degrees +
// normalized-latitude conventions. ----

/// The conformal isometric building block `t = tan(π/4 − φ/2) / [(1 − e sinφ)/(1 +
/// e sinφ)]^{e/2}` (Snyder 15-9) — `t = tan(π/4 − χ/2)`, so `χ = π/2 − 2·arctan t`.
#[must_use]
pub fn conformal_t(phi: f64, e: f64) -> f64 {
    hausdorff_kernels::projection::conformal_t(phi, e)
}

/// Invert the conformal `t` back to geodetic latitude `φ` (radians) by Snyder's
/// fixed-point iteration (7-9). Its `e → 0` limit is the exact closed form
/// `φ = π/2 − 2·arctan t` (proof-backed by `riemann`
/// `Hausdorff/Stereographic.conformalLat_inverse`).
#[must_use]
pub fn lat_from_conformal_t(t: f64, e: f64) -> f64 {
    hausdorff_kernels::projection::lat_from_conformal_t(t, e)
}

/// The authalic-area function `q` (Snyder 3-12) — the integral `∫cosφ/(1−e²sin²φ)²`
/// in closed form; the equal-area analogue of the conformal `t`. `sin ξ = q/q(π/2)`.
#[must_use]
pub fn authalic_q(phi: f64, e: f64) -> f64 {
    hausdorff_kernels::projection::authalic_q(phi, e)
}

/// Invert the authalic `q` back to geodetic latitude `φ` (radians) by Snyder's
/// iteration (3-16).
#[must_use]
pub fn lat_from_authalic_q(q: f64, e: f64) -> f64 {
    // The pole clamp stays a scalar branch in this `f64` wrapper; the fixed-point
    // iteration itself lives in the backend-neutral kernel.
    let q_pole = authalic_q(FRAC_PI_2, e);
    if q.abs() >= q_pole - 1e-12 {
        return FRAC_PI_2.copysign(q);
    }
    hausdorff_kernels::projection::lat_from_authalic_q(q, e)
}

/// Meridian arc length from the equator to latitude `phi` (radians) on the
/// ellipsoid `(a, e²)`, metres — Snyder 3-21, the series behind
/// [`Ellipsoid::meridian_arc_m`].
fn meridian_arc(phi: f64, a: f64, e2: f64) -> f64 {
    a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2.powi(3) / 256.0) * phi
        - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2.powi(3) / 1024.0) * (2.0 * phi).sin()
        + (15.0 * e2 * e2 / 256.0 + 45.0 * e2.powi(3) / 1024.0) * (4.0 * phi).sin()
        - (35.0 * e2.powi(3) / 3072.0) * (6.0 * phi).sin())
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
            assert!(close(
                e.geodetic_from_reduced(e.reduced_latitude(lat)),
                lat,
                1e-9
            ));
            assert!(close(
                e.geodetic_from_geocentric(e.geocentric_latitude(lat)),
                lat,
                1e-9
            ));
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

    /// On a sphere (`f = 0`) the conformal / authalic / rectifying latitudes also
    /// collapse to `φ` (the `e → 0` limit).
    #[test]
    fn sphere_collapses_auxiliary_latitudes() {
        let s = Ellipsoid::SPHERE;
        for &lat in &[0.0, 12.3, 45.0, -60.0, 89.9] {
            assert!(close(s.conformal_latitude(lat), lat, 1e-12));
            assert!(close(s.authalic_latitude(lat), lat, 1e-12));
            assert!(close(s.rectifying_latitude(lat), lat, 1e-12));
        }
    }

    /// Each auxiliary latitude inverts cleanly (WGS84), both hemispheres. The
    /// conformal/authalic inverses iterate to convergence (≤ 1e-9°); the
    /// rectifying inverse is Snyder's footpoint *series* (the TM's), exact only to
    /// the series order (~sub-mm, ≤ 1e-7°).
    #[test]
    fn auxiliary_latitudes_round_trip() {
        let e = Ellipsoid::WGS84;
        for &lat in &[0.0, 12.3, 45.0, -60.0, 80.0, -88.0] {
            assert!(close(
                e.geodetic_from_conformal(e.conformal_latitude(lat)),
                lat,
                1e-9
            ));
            assert!(close(
                e.geodetic_from_authalic(e.authalic_latitude(lat)),
                lat,
                1e-9
            ));
            assert!(close(
                e.geodetic_from_rectifying(e.rectifying_latitude(lat)),
                lat,
                1e-7
            ));
        }
    }

    /// WGS84 auxiliary latitudes vs an **independent numerical integration** of the
    /// area / meridian-arc elements (2e6-step midpoint — not the Snyder series the
    /// methods use). At 45° the conformal squashes most, the authalic least, so
    /// `χ < μ < ξ < φ`.
    #[test]
    fn auxiliary_latitudes_match_independent_integration() {
        let e = Ellipsoid::WGS84;
        for &(lat, chi, xi, mu) in &[
            (30.0, 29.833_682, 29.888_997, 29.875_148),
            (45.0, 44.807_684, 44.871_703, 44.855_682),
            (60.0, 59.833_216, 59.888_786, 59.874_886),
        ] {
            assert!(close(e.conformal_latitude(lat), chi, 1e-4), "χ({lat})");
            assert!(close(e.authalic_latitude(lat), xi, 1e-4), "ξ({lat})");
            assert!(close(e.rectifying_latitude(lat), mu, 1e-4), "μ({lat})");
        }
        // The squash ordering at 45°: χ < μ < ξ < φ.
        assert!(e.conformal_latitude(45.0) < e.rectifying_latitude(45.0));
        assert!(e.rectifying_latitude(45.0) < e.authalic_latitude(45.0));
        assert!(e.authalic_latitude(45.0) < 45.0);
        // Quarter-meridian (equator → pole), WGS84 ≈ 10 001 965.73 m.
        assert!(close(e.meridian_arc_m(90.0), 10_001_965.73, 0.1));
        assert!(close(e.meridian_arc_m(-90.0), -10_001_965.73, 0.1));
    }
}
