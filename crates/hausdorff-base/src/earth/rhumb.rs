//! Rhumb-line (loxodrome) geodesy over [`EARTH_RADIUS_M`]: the constant-bearing
//! course that appears straight on a Mercator chart. The raw `f64`-metre
//! primitives ([`rhumb_distance_m`], [`rhumb_destination_m`]) are the kernels;
//! [`rhumb_distance_meters`] and [`rhumb_destination`] are the `euclid`-typed
//! forms callers should prefer so earth-distance quantities carry units rather
//! than bare `f64`.
//!
//! The math follows the Mercator/isometric-latitude formulation used by
//! georust `geo` (itself adapted from Turf.js and Chris Veness's spherical
//! geodesy tools), so results conform to `geo`'s `Rhumb` metric: longitude
//! deltas over 180° take the shorter rhumb across the antimeridian, and an
//! east-west course (Δψ → 0) falls back to the `cos φ₁` stretch factor.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use euclid::LengthMeters;

use super::{Bearing, EARTH_RADIUS_M, GeoPoint};

/// |Δψ| at or below this is an east-west course: the Mercator stretch factor
/// `Δφ/Δψ` becomes 0/0, so its limit `cos φ₁` is used instead. Matches `geo`.
const EW_COURSE_THRESHOLD: f64 = 1e-11;

/// The shared inverse-problem quantities between two points: `φ₁` and the
/// latitude (`Δφ`), longitude (`Δλ`), and isometric-latitude (`Δψ`) deltas.
struct RhumbInverse {
    phi1: f64,
    delta_lambda: f64,
    delta_phi: f64,
    delta_psi: f64,
}

impl RhumbInverse {
    fn new(from: GeoPoint, to: GeoPoint) -> Self {
        let phi1 = from.lat.to_radians();
        let phi2 = to.lat.to_radians();
        let mut delta_lambda = (to.lon - from.lon).to_radians();
        if delta_lambda > PI {
            delta_lambda -= 2.0 * PI;
        }
        if delta_lambda < -PI {
            delta_lambda += 2.0 * PI;
        }

        let delta_psi = ((phi2 / 2.0 + FRAC_PI_4).tan() / (phi1 / 2.0 + FRAC_PI_4).tan()).ln();
        let delta_phi = phi2 - phi1;

        Self {
            phi1,
            delta_lambda,
            delta_phi,
            delta_psi,
        }
    }

    /// The constant course angle in radians, in `(-π, π]`.
    fn theta(&self) -> f64 {
        self.delta_lambda.atan2(self.delta_psi)
    }

    /// The angular rhumb length in radians (distance over the sphere radius).
    fn delta(&self) -> f64 {
        let q = if self.delta_psi.abs() > EW_COURSE_THRESHOLD {
            self.delta_phi / self.delta_psi
        } else {
            self.phi1.cos()
        };
        (self.delta_phi * self.delta_phi + q * q * self.delta_lambda * self.delta_lambda).sqrt()
    }
}

/// The direct problem in radians: travel the angular distance `delta` from
/// `(λ₁, φ₁)` on the constant course `theta`. A latitude carried past a pole
/// is reflected back, matching `geo`.
fn destination_radians(delta: f64, lambda1: f64, phi1: f64, theta: f64) -> GeoPoint {
    let delta_phi = delta * theta.cos();
    let mut phi2 = phi1 + delta_phi;

    if phi2.abs() > FRAC_PI_2 {
        phi2 = if phi2 > 0.0 { PI - phi2 } else { -PI - phi2 };
    }

    let delta_psi = ((phi2 / 2.0 + FRAC_PI_4).tan() / (phi1 / 2.0 + FRAC_PI_4).tan()).ln();
    let q = if delta_psi.abs() > EW_COURSE_THRESHOLD {
        delta_phi / delta_psi
    } else {
        phi1.cos()
    };

    let lambda2 = lambda1 + (delta * theta.sin()) / q;
    let lon_deg = ((lambda2.to_degrees() + 540.0) % 360.0) - 180.0;
    GeoPoint::from_lat_lon(phi2.to_degrees(), lon_deg)
}

/// Rhumb-line distance between `a` and `b` in **metres**, on a sphere of
/// radius [`EARTH_RADIUS_M`]. The loxodrome inverse — at least the
/// great-circle [`super::haversine_m`] distance, equal only when the course
/// is a meridian or the equator.
#[must_use]
pub fn rhumb_distance_m(a: GeoPoint, b: GeoPoint) -> f64 {
    RhumbInverse::new(a, b).delta() * EARTH_RADIUS_M
}

/// The constant compass bearing of the rhumb line from `from` to `to` —
/// unlike a great circle's initial bearing, it holds for the whole course.
#[must_use]
pub fn rhumb_bearing(from: GeoPoint, to: GeoPoint) -> Bearing {
    Bearing::from_degrees(RhumbInverse::new(from, to).theta().to_degrees())
}

/// The point reached by travelling `distance_m` metres from `origin` on the
/// constant compass `bearing`, on a sphere of radius [`EARTH_RADIUS_M`]. The
/// forward (direct) loxodrome — inverse-pair to [`rhumb_distance_m`] /
/// [`rhumb_bearing`].
///
/// A latitude carried past a pole is reflected back and longitude is
/// normalized into `[-180, 180)`, so results stay valid across the poles and
/// antimeridian.
#[must_use]
pub fn rhumb_destination_m(origin: GeoPoint, bearing: Bearing, distance_m: f64) -> GeoPoint {
    destination_radians(
        distance_m / EARTH_RADIUS_M,
        origin.lon.to_radians(),
        origin.lat.to_radians(),
        bearing.radians(),
    )
}

/// The point `ratio` of the way (by rhumb distance) from `from` to `to` along
/// their rhumb line: `0` is `from`, `1` is `to`, values outside `[0, 1]`
/// extrapolate along the same constant course.
#[must_use]
pub fn rhumb_point_at_ratio(from: GeoPoint, to: GeoPoint, ratio: f64) -> GeoPoint {
    let inverse = RhumbInverse::new(from, to);
    destination_radians(
        ratio * inverse.delta(),
        from.lon.to_radians(),
        inverse.phi1,
        inverse.theta(),
    )
}

/// Rhumb-line distance between `a` and `b` as a typed length — a
/// [`LengthMeters`] wrapper over [`rhumb_distance_m`]. Prefer this over the
/// raw `f64` kernel so callers can't drop or confuse the unit.
#[must_use]
pub fn rhumb_distance_meters(a: GeoPoint, b: GeoPoint) -> LengthMeters {
    LengthMeters::new(rhumb_distance_m(a, b))
}

/// The point reached by travelling `distance` from `origin` on the constant
/// compass `bearing` — a [`LengthMeters`]-typed wrapper over
/// [`rhumb_destination_m`]. Inverse-pair to [`rhumb_distance_meters`] /
/// [`rhumb_bearing`].
#[must_use]
pub fn rhumb_destination(origin: GeoPoint, bearing: Bearing, distance: LengthMeters) -> GeoPoint {
    rhumb_destination_m(origin, bearing, distance.as_meters())
}
