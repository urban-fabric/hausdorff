//! Spherical geodesics over [`EARTH_RADIUS_M`]. The raw `f64`-metre primitives
//! ([`haversine_m`], [`destination_m`]) are the kernels; [`haversine_meters`]
//! and [`destination`] are the `euclid`-typed forms callers should prefer so
//! earth-distance quantities carry units rather than bare `f64`.

use euclid::LengthMeters;

use super::{Bearing, GeoPoint};

/// Mean Earth radius in metres (WGS84 mean radius R₁). The single source of
/// truth for the workspace's spherical geodesics — galois and kolmogorov-geo
/// re-export this rather than redefining it.
pub const EARTH_RADIUS_M: f64 = 6_371_008.8;

/// Earth diameter (`2·`[`EARTH_RADIUS_M`]), folded once so [`haversine_m`]'s hot
/// path is a single multiply rather than `2.0 * R * …`.
const EARTH_DIAMETER_M: f64 = 2.0 * EARTH_RADIUS_M;

/// Great-circle distance between `a` and `b` in **metres**, on a sphere of
/// radius [`EARTH_RADIUS_M`].
#[must_use]
pub fn haversine_m(a: GeoPoint, b: GeoPoint) -> f64 {
    let lat1 = a.lat.to_radians();
    let lat2 = b.lat.to_radians();
    let sin_dlat = ((lat2 - lat1) * 0.5).sin();
    let sin_dlon = ((b.lon - a.lon).to_radians() * 0.5).sin();
    let h = sin_dlat * sin_dlat + lat1.cos() * lat2.cos() * sin_dlon * sin_dlon;
    // 2·R·asin(√h) — equivalent to the atan2 form for h ∈ [0, 1], with one fewer
    // sqrt and a cheaper single-argument inverse trig. (`min(1.0)` guards the
    // asin domain against rounding past 1 for near-antipodal points.)
    EARTH_DIAMETER_M * h.min(1.0).sqrt().asin()
}

/// The point reached by travelling `distance_m` metres from `origin` along the
/// initial compass `bearing`, on a sphere of radius [`EARTH_RADIUS_M`]. The
/// forward (direct) geodesic — inverse-pair to [`haversine_m`].
///
/// Latitude is clamped to the valid `asin` domain and longitude is normalized
/// into `[-180, 180]`, so results stay valid across the poles and antimeridian.
#[must_use]
pub fn destination_m(origin: GeoPoint, bearing: Bearing, distance_m: f64) -> GeoPoint {
    let lat1 = origin.lat.to_radians();
    let lon1 = origin.lon.to_radians();
    let brg = bearing.radians();
    let dr = distance_m / EARTH_RADIUS_M;

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
