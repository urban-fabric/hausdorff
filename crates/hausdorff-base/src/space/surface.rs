//! `surface` — the [`SurfaceOfRevolution`] trait: one generic curved-space
//! geometry (the local metric of a surface spun from a meridian profile about an
//! axis), shared by the sphere, the ellipsoid, the torus, the paraboloid, ….
//!
//! A surface of revolution `r(t, λ) = (ρ(t) cos λ, ρ(t) sin λ, z(t))` has a
//! **diagonal first fundamental form** — meridians meet parallels at right angles
//! — fully described by two functions of the meridian parameter `t`: the parallel
//! radius `√G = ρ` and the meridian arc speed `√E = √(ρ′² + z′²)`. Everything else
//! (the metric, the area element, surface area) follows generically and is
//! written **once** here as default methods. The trait's laws are machine-checked
//! in `riemann` `Hausdorff/Spheroid` over an abstract profile `(ρ, z, ρ′, z′)`:
//! `meridian_perp_parallel` (F = 0), `parallel_normSq` (G = ρ²), `meridian_normSq`
//! (E = ρ′²+z′²), `area_element_eq` (`‖rφ × rλ‖² = EG`).
//!
//! Instances: [`Torus`] here, and `Ellipsoid` (in [`crate::earth`]) — the sphere
//! being `Ellipsoid::SPHERE`.

/// A surface of revolution, described by its two first-fundamental-form radii as
/// functions of the meridian parameter `t` (degrees — geodetic latitude for the
/// ellipsoid, the tube angle for a torus).
pub trait SurfaceOfRevolution {
    /// `√E` — the **meridian arc speed**: metres of meridian arc per radian of
    /// the parameter `t`. (`riemann` `Spheroid.meridian_normSq`: `E = ρ′²+z′²`.)
    fn meridian_speed(&self, t_deg: f64) -> f64;

    /// `√G = ρ` — the **parallel radius**: the distance from the polar axis, so
    /// the circle of latitude has circumference `2πρ`. (`riemann`
    /// `Spheroid.parallel_normSq`: `G = ρ²`.)
    fn parallel_radius(&self, t_deg: f64) -> f64;

    /// The `t` range (degrees) over which the closed surface is swept (e.g.
    /// `(-90, 90)` for the ellipsoid, `(0, 360)` for a torus).
    fn meridian_domain_deg(&self) -> (f64, f64);

    /// The **first fundamental form** `(E, F, G)`. `F = 0` — the meridian/parallel
    /// net is orthogonal (`riemann` `Spheroid.meridian_perp_parallel`); `E =
    /// meridian_speed²`, `G = parallel_radius²`.
    fn first_fundamental_form(&self, t_deg: f64) -> (f64, f64, f64) {
        let (m, r) = (self.meridian_speed(t_deg), self.parallel_radius(t_deg));
        (m * m, 0.0, r * r)
    }

    /// The **surface-area element** `√(EG − F²) = √(EG) = meridian_speed ·
    /// parallel_radius` (the diagonal form collapses the Gram determinant to
    /// `E·G`; `riemann` `Spheroid.area_element_eq`).
    fn area_element(&self, t_deg: f64) -> f64 {
        self.meridian_speed(t_deg) * self.parallel_radius(t_deg)
    }

    /// The circumference `2π · parallel_radius` of the parallel at `t`.
    fn parallel_circumference(&self, t_deg: f64) -> f64 {
        std::f64::consts::TAU * self.parallel_radius(t_deg)
    }

    /// The area of the zone (band) between parameters `t0` and `t1`, by numeric
    /// midpoint integration `2π ∫ area_element dt` over `steps` panels. One
    /// generic integrator over the proof-backed area element serves every surface.
    fn zonal_area(&self, t0_deg: f64, t1_deg: f64, steps: usize) -> f64 {
        let steps = steps.max(1);
        let dt = (t1_deg - t0_deg) / steps as f64;
        let dt_rad = dt.to_radians();
        let mut sum = 0.0;
        for i in 0..steps {
            let t = t0_deg + (i as f64 + 0.5) * dt;
            sum += self.area_element(t) * dt_rad;
        }
        (std::f64::consts::TAU * sum).abs()
    }

    /// The **total surface area** — [`zonal_area`](Self::zonal_area) over the full
    /// [`meridian_domain_deg`](Self::meridian_domain_deg).
    fn total_area(&self, steps: usize) -> f64 {
        let (a, b) = self.meridian_domain_deg();
        self.zonal_area(a, b, steps)
    }
}

/// A **torus** of revolution: a tube of radius `minor` whose centre circles the
/// polar axis at radius `major`. The meridian profile is the tube circle
/// `ρ = major + minor·cos t`, `z = minor·sin t` (`t` the tube angle), a genuine
/// non-ellipsoidal [`SurfaceOfRevolution`] — the abstract first-fundamental-form
/// theorems apply directly (`riemann` `Spheroid.{torus_meridian_normSq,
/// torus_area_element}`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Torus {
    /// Major radius `R` (axis → tube centre).
    pub major: f64,
    /// Minor radius `r` (tube radius).
    pub minor: f64,
}

impl SurfaceOfRevolution for Torus {
    /// `√E = r` (constant), `riemann` `torus_meridian_normSq`.
    fn meridian_speed(&self, _t_deg: f64) -> f64 {
        self.minor
    }

    /// `ρ = R + r·cos t`.
    fn parallel_radius(&self, t_deg: f64) -> f64 {
        self.major + self.minor * t_deg.to_radians().cos()
    }

    /// The tube closes over a full `0..360°` turn.
    fn meridian_domain_deg(&self) -> (f64, f64) {
        (0.0, 360.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn torus_first_fundamental_form_and_area() {
        let t = Torus {
            major: 3.0,
            minor: 1.0,
        };
        // F = 0 (orthogonal); E = r² = 1; G = ρ² = (R + r cos t)².
        let (e, f, g) = t.first_fundamental_form(0.0);
        assert_eq!(f, 0.0);
        assert!((e - 1.0).abs() < 1e-12);
        assert!((g - (3.0 + 1.0_f64).powi(2)).abs() < 1e-12); // ρ = 4 at t=0
        // Total area of a torus is the closed form 4π²Rr.
        let exact = 4.0 * std::f64::consts::PI.powi(2) * t.major * t.minor;
        let area = t.total_area(20_000);
        assert!(
            (area - exact).abs() / exact < 1e-7,
            "area {area} vs {exact}"
        );
        // The parallel circumference at the outer equator (t=0): 2π·4.
        assert!((t.parallel_circumference(0.0) - std::f64::consts::TAU * 4.0).abs() < 1e-12);
    }
}
