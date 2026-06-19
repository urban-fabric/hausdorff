//! Great-circle Earth formulas over any [`LazyReal`] scalar.
//!
//! This module is the backend-portable sibling of the concrete [`GeoPoint`]
//! geodesic API. It takes latitude/longitude and bearings as degree-valued
//! scalars directly, and returns scalars instead of constructing Hausdorff data
//! records. That keeps the storage/serde `GeoPoint` contract concrete while the
//! formula kernels can be traced by tensor, columnar, and scalar backends.

use euclid::numeric::LazyReal;

use super::EARTH_RADIUS_M;

/// Great-circle distance in metres on a sphere with the supplied `diameter_m`.
///
/// Coordinate arguments are latitude/longitude in decimal degrees.
#[must_use]
pub fn haversine_diameter_m<T: LazyReal>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
    diameter_m: T,
) -> T {
    let half = from_lat_deg.scalar_like(0.5);
    let zero = from_lat_deg.scalar_like(0.0);
    let one = from_lat_deg.scalar_like(1.0);

    let lat1 = from_lat_deg.to_radians();
    let lat2 = to_lat_deg.to_radians();
    let sin_dlat = ((lat2.clone() - lat1.clone()) * half.clone()).sin();
    let sin_dlon = ((to_lon_deg - from_lon_deg).to_radians() * half).sin();
    let h = sin_dlat.clone() * sin_dlat + lat1.cos() * lat2.cos() * sin_dlon.clone() * sin_dlon;

    diameter_m * h.clamp(zero, one).sqrt().asin()
}

/// Great-circle distance in metres on a sphere of `radius_m`.
///
/// Coordinate arguments are latitude/longitude in decimal degrees.
#[must_use]
pub fn haversine_on_m<T: LazyReal>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
    radius_m: T,
) -> T {
    let two = from_lat_deg.scalar_like(2.0);
    haversine_diameter_m(
        from_lat_deg,
        from_lon_deg,
        to_lat_deg,
        to_lon_deg,
        radius_m * two,
    )
}

/// Great-circle distance in metres on Earth.
///
/// Coordinate arguments are latitude/longitude in decimal degrees.
#[must_use]
pub fn haversine_m<T: LazyReal>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
) -> T {
    let diameter_m = from_lat_deg.scalar_like(2.0 * EARTH_RADIUS_M);
    haversine_diameter_m(
        from_lat_deg,
        from_lon_deg,
        to_lat_deg,
        to_lon_deg,
        diameter_m,
    )
}

/// Initial great-circle compass bearing in degrees, normalized into `[0, 360)`.
///
/// Coordinate arguments are latitude/longitude in decimal degrees.
#[must_use]
pub fn initial_bearing_degrees<T: LazyReal>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
) -> T {
    let full_turn = from_lat_deg.scalar_like(360.0);

    let lat1 = from_lat_deg.to_radians();
    let lat2 = to_lat_deg.to_radians();
    let dlon = (to_lon_deg - from_lon_deg).to_radians();
    let y = dlon.clone().sin() * lat2.clone().cos();
    let x = lat1.clone().cos() * lat2.clone().sin() - lat1.sin() * lat2.cos() * dlon.cos();

    y.atan2(x).to_degrees().rem_euclid(full_turn)
}

/// Destination latitude/longitude in decimal degrees after travelling
/// `distance_m` on a sphere of `radius_m` from an origin along `bearing_deg`.
#[must_use]
pub fn destination_on_m<T: LazyReal>(
    origin_lat_deg: T,
    origin_lon_deg: T,
    bearing_deg: T,
    distance_m: T,
    radius_m: T,
) -> (T, T) {
    let minus_one = origin_lat_deg.scalar_like(-1.0);
    let one = origin_lat_deg.scalar_like(1.0);
    let lon_wrap_offset = origin_lat_deg.scalar_like(540.0);
    let full_turn = origin_lat_deg.scalar_like(360.0);
    let half_turn = origin_lat_deg.scalar_like(180.0);

    let lat1 = origin_lat_deg.to_radians();
    let lon1 = origin_lon_deg.to_radians();
    let brg = bearing_deg.to_radians();
    let dr = distance_m / radius_m;

    let sin_lat2 = (lat1.clone().sin() * dr.clone().cos()
        + lat1.clone().cos() * dr.clone().sin() * brg.clone().cos())
    .clamp(minus_one, one);
    let lat2 = sin_lat2.clone().asin();
    let lon2 = lon1
        + (brg.clone().sin() * dr.clone().sin() * lat1.clone().cos())
            .atan2(dr.cos() - lat1.sin() * sin_lat2);

    let lon_deg = (lon2.to_degrees() + lon_wrap_offset).rem_euclid(full_turn) - half_turn;
    (lat2.to_degrees(), lon_deg)
}

/// Destination latitude/longitude in decimal degrees after travelling
/// `distance_m` on Earth from an origin along `bearing_deg`.
#[must_use]
pub fn destination_m<T: LazyReal>(
    origin_lat_deg: T,
    origin_lon_deg: T,
    bearing_deg: T,
    distance_m: T,
) -> (T, T) {
    let radius_m = origin_lat_deg.scalar_like(EARTH_RADIUS_M);
    destination_on_m(
        origin_lat_deg,
        origin_lon_deg,
        bearing_deg,
        distance_m,
        radius_m,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::earth::{self, Bearing, GeoPoint};

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual {actual} != expected {expected} within {tolerance}"
        );
    }

    #[test]
    fn haversine_lazy_matches_concrete() {
        let from = GeoPoint::from_lat_lon(37.7749, -122.4194);
        let to = GeoPoint::from_lat_lon(40.7128, -74.0060);

        let actual = haversine_m(from.lat, from.lon, to.lat, to.lon);
        let expected = earth::haversine_m(from, to);

        assert_close(actual, expected, 1.0e-9);
    }

    #[test]
    fn haversine_lazy_runs_for_f32() {
        let distance = haversine_m(37.7749_f32, -122.4194_f32, 40.7128_f32, -74.0060_f32);

        assert!(distance.is_finite());
    }

    #[test]
    fn initial_bearing_lazy_matches_concrete() {
        let from = GeoPoint::from_lat_lon(37.7749, -122.4194);
        let to = GeoPoint::from_lat_lon(40.7128, -74.0060);

        let actual = initial_bearing_degrees(from.lat, from.lon, to.lat, to.lon);
        let expected = Bearing::initial(from, to).degrees();

        assert_close(actual, expected, 1.0e-12);
    }

    #[test]
    fn destination_lazy_matches_concrete() {
        let origin = GeoPoint::from_lat_lon(37.7749, -122.4194);
        let bearing = Bearing::from_degrees(73.0);
        let distance_m = 42_000.0;

        let (actual_lat, actual_lon) =
            destination_m(origin.lat, origin.lon, bearing.degrees(), distance_m);
        let expected = earth::destination_m(origin, bearing, distance_m);

        assert_close(actual_lat, expected.lat, 1.0e-12);
        assert_close(actual_lon, expected.lon, 1.0e-12);
    }
}
