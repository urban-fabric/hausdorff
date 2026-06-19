//! Spherical geodesics over [`EARTH_RADIUS_M`]. The raw `f64`-metre primitives
//! ([`haversine_m`], [`destination_m`]) are the kernels; [`haversine_meters`]
//! and [`destination`] are the `euclid`-typed forms callers should prefer so
//! earth-distance quantities carry units rather than bare `f64`.

use euclid::LengthMeters;
pub use hausdorff_kernels::earth::EARTH_RADIUS_M;

use super::{Bearing, GeoPoint};

/// Earth diameter (`2·`[`EARTH_RADIUS_M`]), folded once so [`haversine_m`]'s hot
/// path is a single multiply rather than `2.0 * R * …`.
const EARTH_DIAMETER_M: f64 = 2.0 * EARTH_RADIUS_M;

/// The haversine kernel parameterized by sphere *diameter* — `inline(always)`
/// so the Earth wrapper keeps its folded [`EARTH_DIAMETER_M`] constant and the
/// hot path stays a single multiply.
///
/// Proof-backed by `riemann` `Hausdorff/Spherical` (`haversine_correct`: the
/// haversine `h` is `sin²(Δσ/2)`, so `2·asin√h` is the central angle / great-circle
/// distance; `chord_sq_eq_four_haversine`) and `Hausdorff/HaversineBounds`
/// (`h ∈ [0, 1]`, distance `∈ [0, πR]` — so the `h.min(1.0)` clamp is a no-op on
/// the true value).
#[inline(always)]
fn haversine_diameter_m(a: GeoPoint, b: GeoPoint, diameter_m: f64) -> f64 {
    let lat1 = a.lat.to_radians();
    let lat2 = b.lat.to_radians();
    let sin_dlat = ((lat2 - lat1) * 0.5).sin();
    let sin_dlon = ((b.lon - a.lon).to_radians() * 0.5).sin();
    let h = sin_dlat * sin_dlat + lat1.cos() * lat2.cos() * sin_dlon * sin_dlon;
    diameter_m * h.min(1.0).sqrt().asin()
}

/// Great-circle distance between `a` and `b` in **metres**, on a sphere of
/// radius [`EARTH_RADIUS_M`].
#[must_use]
pub fn haversine_m(a: GeoPoint, b: GeoPoint) -> f64 {
    haversine_diameter_m(a, b, EARTH_DIAMETER_M)
}

/// Great-circle distance between `a` and `b` in metres on a sphere of
/// `radius_m` metres — the planet-flexible kernel behind
/// [`crate::earth::sphere::Sphere::distance`].
#[must_use]
pub fn haversine_on_m(a: GeoPoint, b: GeoPoint, radius_m: f64) -> f64 {
    haversine_diameter_m(a, b, 2.0 * radius_m)
}

/// The point reached by travelling `distance_m` metres from `origin` along the
/// initial compass `bearing`, on a sphere of radius [`EARTH_RADIUS_M`]. The
/// forward (direct) geodesic — inverse-pair to [`haversine_m`].
///
/// Latitude is clamped to the valid `asin` domain and longitude is normalized
/// into `[-180, 180]`, so results stay valid across the poles and antimeridian.
#[must_use]
pub fn destination_m(origin: GeoPoint, bearing: Bearing, distance_m: f64) -> GeoPoint {
    destination_on_m(origin, bearing, distance_m, EARTH_RADIUS_M)
}

/// [`destination_m`] on a sphere of `radius_m` metres — the planet-flexible
/// kernel behind [`crate::earth::sphere::Sphere::destination`].
#[must_use]
pub fn destination_on_m(
    origin: GeoPoint,
    bearing: Bearing,
    distance_m: f64,
    radius_m: f64,
) -> GeoPoint {
    let lat1 = origin.lat.to_radians();
    let lon1 = origin.lon.to_radians();
    let brg = bearing.radians();
    let dr = distance_m / radius_m;

    let sin_lat2 = (lat1.sin() * dr.cos() + lat1.cos() * dr.sin() * brg.cos()).clamp(-1.0, 1.0);
    let lat2 = sin_lat2.asin();
    let lon2 = lon1 + (brg.sin() * dr.sin() * lat1.cos()).atan2(dr.cos() - lat1.sin() * sin_lat2);

    let lon_deg = (lon2.to_degrees() + 540.0).rem_euclid(360.0) - 180.0;
    GeoPoint::from_lat_lon(lat2.to_degrees(), lon_deg)
}

/// Great-circle distance between `a` and `b` as a typed length — a
/// [`LengthMeters`] wrapper over [`haversine_m`]. Prefer this over the raw
/// `f64` kernel so callers can't drop or confuse the unit.
#[must_use]
pub fn haversine_meters(a: GeoPoint, b: GeoPoint) -> LengthMeters {
    LengthMeters::new(haversine_m(a, b))
}

/// The point reached by travelling `distance` from `origin` along the initial
/// compass `bearing` — a [`LengthMeters`]-typed wrapper over [`destination_m`].
/// Inverse-pair to [`haversine_meters`].
#[must_use]
pub fn destination(origin: GeoPoint, bearing: Bearing, distance: LengthMeters) -> GeoPoint {
    destination_m(origin, bearing, distance.as_meters())
}
