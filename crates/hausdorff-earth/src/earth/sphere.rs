//! Planet-flexible metric layer: turn the generic, unit-less spherical geometry
//! (angles from `hausdorff-base`'s `space::s1`, areas in steradians from `hausdorff-sphere`'s `space::s2`) into
//! physical lengths and areas on a sphere of a chosen radius.
//!
//! `s1`/`s2` live on the *unit* circle/sphere and carry no radius; [`crate::earth`]
//! is Earth-specific. [`Sphere`] is the general bridge — pick a radius (Earth,
//! Mars, the Moon, or any custom value) and it converts angles → metres and
//! steradians → m². This is the only place a radius (and hence a physical unit)
//! enters; the geometry itself stays unit-less and planet-agnostic.

use crate::earth::{AreaSquareMeters, Bearing, Ellipsoid, GeoPoint, LengthMeters};
use hausdorff_base::earth::latitude::authalic_q;
use hausdorff_base::space::s1::Angle;
use hausdorff_sphere::space::s2::{S2Point, polygon_signed_area_vos};

/// A sphere of a given radius — the bridge from unit-less spherical geometry to
/// physical distances and areas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere {
    radius_m: f64,
}

impl Sphere {
    /// Earth — GRS80 / WGS84 mean radius R₁, matching [`crate::earth::EARTH_RADIUS_M`].
    pub const EARTH: Sphere = Sphere {
        radius_m: crate::earth::EARTH_RADIUS_M,
    };
    /// Mars — IAU mean radius.
    pub const MARS: Sphere = Sphere {
        radius_m: 3_389_500.0,
    };
    /// The Moon — mean radius.
    pub const MOON: Sphere = Sphere {
        radius_m: 1_737_400.0,
    };

    /// A sphere of the given radius.
    #[must_use]
    pub fn new(radius: LengthMeters) -> Self {
        Self {
            radius_m: radius.value(),
        }
    }

    /// A sphere of `radius_m` metres.
    #[must_use]
    pub const fn from_meters(radius_m: f64) -> Self {
        Self { radius_m }
    }

    /// The radius.
    #[must_use]
    pub fn radius(self) -> LengthMeters {
        LengthMeters::new(self.radius_m)
    }

    /// Arc length of a central `angle`: `radius × angle`. The general form of a
    /// great-circle distance — `Sphere::EARTH.length(...)` is the Earth case.
    ///
    /// The radian measure of `angle` is dimensionless, so it scales the radius
    /// through `LengthMeters`' scalar multiplication.
    #[must_use]
    pub fn length(self, angle: Angle) -> LengthMeters {
        self.radius() * angle.radians()
    }

    /// Surface area subtended by a solid angle in `steradians`: `radius² ×
    /// steradians` (e.g. a full sphere is `4π` steradians).
    ///
    /// `radius × radius` is the `LengthMeters × LengthMeters =
    /// AreaSquareMeters` algebraic connection — the area subtended by one
    /// steradian — then scaled by the dimensionless solid angle. Both carriers
    /// are Hausdorff's own, implementing `euclid-units`' dimension contracts.
    #[must_use]
    pub fn area(self, steradians: f64) -> AreaSquareMeters {
        let area_per_steradian = self.radius() * self.radius();
        area_per_steradian * steradians
    }

    /// Surface area subtended by a solid angle in `steradians` as raw m² — the
    /// `f64` kernel behind [`area`](Self::area), for callers that want the bare
    /// scalar rather than the typed [`AreaSquareMeters`].
    #[must_use]
    pub fn area_m2(self, steradians: f64) -> f64 {
        self.area(steradians).value()
    }

    /// **Native spherical area (m²) of a WGS84 polygon `ring`** — convert each
    /// `(lat, lon)` to its unit-sphere point and sum the van Oosterom fan
    /// excesses, then scale by `radius²`. No map projection: the geodesic
    /// (great-circle) polygon is measured directly on the sphere.
    ///
    /// Returns the **unsigned** area; a ring of fewer than three points is `0`.
    /// This reads geodetic latitude as spherical latitude (the spherical model);
    /// for the true datum area on the spheroid use
    /// [`Ellipsoid::ring_perimeter_area_m`](crate::earth::Ellipsoid::ring_perimeter_area_m).
    /// Proof-backed by `riemann` Hausdorff/SphereArea (`polygon_signed_area_vos`,
    /// the fan sum of van Oosterom triangle excesses).
    #[must_use]
    pub fn polygon_area_m2(self, ring: &[GeoPoint]) -> f64 {
        let pts: Vec<S2Point> = ring.iter().map(|&p| S2Point::from_geo_point(p)).collect();
        self.area_m2(polygon_signed_area_vos(&pts).abs())
    }

    /// Great-circle (haversine) distance between two lat/lon points on this
    /// sphere. `Sphere::EARTH.distance(...)` matches
    /// [`crate::earth::haversine_meters`]; other radii give the same geodesic
    /// on that planet. The points are read as lat/lon degrees on *this* sphere.
    #[must_use]
    pub fn distance(self, a: GeoPoint, b: GeoPoint) -> LengthMeters {
        LengthMeters::new(self.distance_m(a, b))
    }

    /// [`distance`](Self::distance) as raw metres — the `f64` kernel.
    #[must_use]
    pub fn distance_m(self, a: GeoPoint, b: GeoPoint) -> f64 {
        crate::earth::haversine_on_m(a, b, self.radius_m)
    }

    /// The point reached by travelling `distance` from `origin` along the
    /// initial compass `bearing` on this sphere — the forward (direct)
    /// geodesic, inverse-pair to [`distance`](Self::distance).
    #[must_use]
    pub fn destination(
        self,
        origin: GeoPoint,
        bearing: Bearing,
        distance: LengthMeters,
    ) -> GeoPoint {
        self.destination_m(origin, bearing, distance.value())
    }

    /// [`destination`](Self::destination) with the distance as raw metres.
    #[must_use]
    pub fn destination_m(self, origin: GeoPoint, bearing: Bearing, distance_m: f64) -> GeoPoint {
        crate::earth::destination_on_m(origin, bearing, distance_m, self.radius_m)
    }
}

/// The **authalic (equal-area) sphere** of an [`Ellipsoid`] — the sphere of radius
/// `Rq` with the *same total surface area* as the ellipsoid. The authalic latitude
/// ξ ([`Ellipsoid::authalic_latitude`]) is the area-preserving map onto it, so a
/// polygon warped `φ → ξ` and measured here has its **true ellipsoidal area**.
///
/// This is the proof-clean route to datum polygon area: it goes through the
/// **proven** spherical-excess kernel (`riemann` Hausdorff/SphereArea,
/// [`Sphere::polygon_area_m2`]'s van Oosterom fan sum) rather than the
/// numerically-validated-but-unproven Karney series
/// ([`Ellipsoid::ring_perimeter_area_m`]), which it reproduces to ≈ 1e-4 %.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AuthalicSphere {
    ellipsoid: Ellipsoid,
    radius_m: f64,
}

impl AuthalicSphere {
    /// The authalic sphere of `ellipsoid`: radius `Rq = a·√(q_p/2)`, where
    /// `q_p = q(π/2)` is the authalic function at the pole, so `4π Rq²` is exactly
    /// the ellipsoid's total surface area. (For a sphere, `Rq = a`.)
    #[must_use]
    pub fn of(ellipsoid: Ellipsoid) -> Self {
        let radius_m = if ellipsoid.f == 0.0 {
            ellipsoid.a
        } else {
            let e = ellipsoid.e2().sqrt();
            ellipsoid.a * (authalic_q(std::f64::consts::FRAC_PI_2, e) / 2.0).sqrt()
        };
        Self {
            ellipsoid,
            radius_m,
        }
    }

    /// The authalic radius `Rq`.
    #[must_use]
    pub fn radius(self) -> LengthMeters {
        LengthMeters::new(self.radius_m)
    }

    /// **True ellipsoidal area (m²) of a WGS84 polygon `ring`**, via the authalic
    /// sphere. Each vertex's geodetic latitude is warped to its authalic latitude
    /// ξ (longitude unchanged), the van Oosterom spherical excesses are summed,
    /// and the result is scaled by `Rq²`. Because `φ → ξ` is area-preserving, this
    /// reproduces the datum polygon area through the **proven** spherical-excess
    /// path — agreeing with the Karney geodesic area
    /// ([`Ellipsoid::ring_perimeter_area_m`]) to ≈ 1e-4 % while staying inside the
    /// closed-form proof boundary.
    ///
    /// Returns the **unsigned** area; a ring of fewer than three vertices is `0`.
    /// The ring closes itself (the edge from the last vertex back to the first).
    #[must_use]
    pub fn polygon_area_m2(self, ring: &[GeoPoint]) -> f64 {
        let pts: Vec<S2Point> = ring
            .iter()
            .map(|&p| {
                let xi = self.ellipsoid.authalic_latitude(p.lat);
                S2Point::from_geo_point(GeoPoint::from_lat_lon(xi, p.lon))
            })
            .collect();
        Sphere::from_meters(self.radius_m).area_m2(polygon_signed_area_vos(&pts).abs())
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use super::*;
    use crate::earth::haversine_m;

    #[test]
    fn earth_distance_matches_haversine() {
        let a = GeoPoint::from_lat_lon(48.85, 2.35);
        let b = GeoPoint::from_lat_lon(51.5, -0.12);
        assert_eq!(Sphere::EARTH.distance_m(a, b), haversine_m(a, b));
    }

    #[test]
    fn quarter_meridian_scales_with_radius() {
        let equator = GeoPoint::from_lat_lon(0.0, 0.0);
        let pole = GeoPoint::from_lat_lon(90.0, 0.0);
        for sphere in [Sphere::EARTH, Sphere::MARS, Sphere::MOON] {
            let quarter = sphere.distance_m(equator, pole);
            let expected = FRAC_PI_2 * sphere.radius().value();
            assert!((quarter - expected).abs() < 1e-6, "{quarter} vs {expected}");
        }
    }

    #[test]
    fn destination_inverts_distance_on_mars() {
        let origin = GeoPoint::from_lat_lon(-14.57, 175.47); // near Olympus Mons' antipode
        let bearing = Bearing::from_degrees(42.0);
        let dest = Sphere::MARS.destination_m(origin, bearing, 250_000.0);
        let back = Sphere::MARS.distance_m(origin, dest);
        assert!((back - 250_000.0).abs() < 1e-6, "round trip {back}");
    }

    #[test]
    fn spherical_polygon_area_from_wgs84() {
        use hausdorff_base::earth::Ellipsoid;
        let g = GeoPoint::from_lat_lon;
        // A ~1°×1° box near the equator (read natively on the sphere).
        let bx = [g(0.0, 0.0), g(0.0, 1.0), g(1.0, 1.0), g(1.0, 0.0)];
        let sph = Sphere::EARTH.polygon_area_m2(&bx);
        // Right ballpark (~1.2e10 m²) and positive.
        assert!(sph > 1.0e10 && sph < 1.5e10, "spherical area {sph}");
        // Within ~2% of the ellipsoidal (Karney) area — the same geodesic
        // polygon; the gap is the sphere-vs-spheroid flattening.
        let (_, ell) = Ellipsoid::WGS84.ring_perimeter_area_unsigned_m(&bx);
        assert!(
            (sph - ell).abs() / ell < 0.02,
            "spherical {sph} vs ellipsoidal {ell}"
        );
        // A degenerate ring is zero.
        assert_eq!(Sphere::EARTH.polygon_area_m2(&bx[..2]), 0.0);
    }

    #[test]
    fn authalic_radius_and_total_area() {
        // WGS84 authalic radius ≈ 6 371 007.18 m; 4π Rq² is the total surface area.
        let au = AuthalicSphere::of(Ellipsoid::WGS84);
        let rq = au.radius().value();
        assert!((rq - 6_371_007.181).abs() < 1e-2, "Rq = {rq}");
        let total = Ellipsoid::WGS84.total_area_m2();
        let rel = (4.0 * std::f64::consts::PI * rq * rq - total).abs() / total;
        assert!(rel < 1e-7, "4πRq² vs total area rel = {rel}");
        // A sphere is its own authalic sphere (Rq = a).
        assert_eq!(
            AuthalicSphere::of(Ellipsoid::SPHERE).radius().value(),
            Ellipsoid::SPHERE.a
        );
    }

    #[test]
    fn authalic_area_matches_karney_and_beats_the_plain_sphere() {
        let g = GeoPoint::from_lat_lon;
        let au = AuthalicSphere::of(Ellipsoid::WGS84);
        // A graticule cell, a mid-latitude box, and a continent-scale quad.
        for ring in [
            vec![g(45.0, 10.0), g(45.0, 11.0), g(46.0, 11.0), g(46.0, 10.0)],
            vec![g(30.0, 0.0), g(30.0, 10.0), g(40.0, 10.0), g(40.0, 0.0)],
            vec![g(0.0, 0.0), g(0.0, 40.0), g(50.0, 40.0), g(50.0, 0.0)],
        ] {
            let (_, karney) = Ellipsoid::WGS84.ring_perimeter_area_unsigned_m(&ring);
            let authalic = au.polygon_area_m2(&ring);
            let plain = Sphere::EARTH.polygon_area_m2(&ring);
            // Authalic reproduces the geodesic (Karney) area to < 0.01 % …
            assert!(
                (authalic - karney).abs() / karney < 1e-4,
                "authalic {authalic} vs karney {karney}"
            );
            // … and is strictly closer to it than the un-warped spherical area.
            assert!(
                (authalic - karney).abs() < (plain - karney).abs(),
                "authalic {authalic} should beat plain {plain} (karney {karney})"
            );
        }
        // A degenerate ring is zero.
        assert_eq!(au.polygon_area_m2(&[g(0.0, 0.0), g(1.0, 1.0)]), 0.0);
    }

    #[test]
    fn authalic_sphere_of_a_sphere_is_the_plain_sphere() {
        // For a true sphere ξ = φ and Rq = a, so the authalic area equals the
        // native spherical area at that radius (the construction degenerates).
        let g = GeoPoint::from_lat_lon;
        let ring = [g(10.0, 0.0), g(10.0, 20.0), g(30.0, 20.0), g(30.0, 0.0)];
        let e = Ellipsoid::SPHERE;
        let authalic = AuthalicSphere::of(e).polygon_area_m2(&ring);
        let sphere = Sphere::from_meters(e.a).polygon_area_m2(&ring);
        let rel = (authalic - sphere).abs() / sphere;
        assert!(rel < 1e-12, "authalic {authalic} vs sphere {sphere}");
    }
}
