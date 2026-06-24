//! Spherical geodesics over [`EARTH_RADIUS_M`]. The raw `f64`-metre primitives
//! ([`haversine_m`], [`destination_m`]) are thin `GeoPoint`-typed wrappers over
//! the backend-neutral [`hausdorff_kernels::earth`] kernels — written once over
//! any [`euclid::numeric::Real`], so the same formula runs on `f64`, `f32`,
//! and tensors. [`haversine_meters`] and [`destination`] are the `euclid`-typed
//! forms callers should prefer so earth-distance quantities carry units rather
//! than bare `f64`.

pub use hausdorff_kernels::earth::EARTH_RADIUS_M;

use super::{Bearing, GeoPoint, LengthMeters};

/// Earth diameter (`2·`[`EARTH_RADIUS_M`]), folded once so [`haversine_m`]'s hot
/// path is a single multiply rather than `2.0 * R * …`.
const EARTH_DIAMETER_M: f64 = 2.0 * EARTH_RADIUS_M;

/// The haversine kernel parameterized by sphere *diameter*, delegating to the
/// generic [`hausdorff_kernels::earth::haversine_diameter_m`]. `inline(always)`
/// so the Earth wrapper keeps its folded [`EARTH_DIAMETER_M`] constant and the
/// hot path stays a single multiply.
///
/// Proof-backed by `riemann` `Hausdorff/Spherical` (`haversine_correct`: the
/// haversine `h` is `sin²(Δσ/2)`, so `2·asin√h` is the central angle / great-circle
/// distance; `chord_sq_eq_four_haversine`) and `Hausdorff/HaversineBounds`
/// (`h ∈ [0, 1]`, distance `∈ [0, πR]` — so the `[0, 1]` clamp is a no-op on the
/// true value).
#[inline(always)]
fn haversine_diameter_m(a: GeoPoint, b: GeoPoint, diameter_m: f64) -> f64 {
    hausdorff_kernels::earth::haversine_diameter_m(a.lat, a.lon, b.lat, b.lon, diameter_m)
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
    let (lat, lon) = hausdorff_kernels::earth::destination_on_m(
        origin.lat,
        origin.lon,
        bearing.degrees(),
        distance_m,
        radius_m,
    );
    GeoPoint::from_lat_lon(lat, lon)
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
    destination_m(origin, bearing, distance.value())
}
