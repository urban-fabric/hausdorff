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

/// Great-circle distance between `a` and `b` in **metres**, on a sphere of
/// radius [`EARTH_RADIUS_M`].
#[must_use]
pub fn haversine_m(a: GeoPoint, b: GeoPoint) -> f64 {
    let lat1 = a.lat.to_radians();
    let lat2 = b.lat.to_radians();
    let dlat = lat2 - lat1;
    let dlon = (b.lon - a.lon).to_radians();
    let h = (dlat * 0.5).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon * 0.5).sin().powi(2);
    2.0 * EARTH_RADIUS_M * h.sqrt().atan2((1.0 - h).sqrt())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn destination_is_inverse_of_haversine() {
        let berlin = GeoPoint::from_lat_lon(52.52, 13.405);
        for deg in [0.0, 37.0, 90.0, 180.0, 271.5, 359.0] {
            let dest = destination_m(berlin, Bearing::from_degrees(deg), 1_500.0);
            assert!(
                (haversine_m(berlin, dest) - 1_500.0).abs() < 0.5,
                "deg {deg}"
            );
        }
    }
}
