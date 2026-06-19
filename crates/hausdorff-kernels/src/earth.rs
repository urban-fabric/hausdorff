//! Earth and spherical geodesic kernels over any [`LazyReal`] scalar.

use euclid::numeric::LazyReal;

/// Mean Earth radius in metres (WGS84 mean radius R1).
pub const EARTH_RADIUS_M: f64 = 6_371_008.8;

/// A latitude/longitude input whose coordinates are decimal degrees.
///
/// Backend adapters can implement this for expression columns, tensor handles,
/// or small wrapper structs without depending on Hausdorff's concrete `GeoPoint`
/// storage contract.
pub trait LatLon<T: LazyReal> {
    /// Latitude in decimal degrees.
    fn lat_deg(&self) -> T;

    /// Longitude in decimal degrees.
    fn lon_deg(&self) -> T;
}

/// A bearing input whose value is decimal degrees clockwise from north.
pub trait BearingDegrees<T: LazyReal> {
    /// Bearing in decimal degrees.
    fn degrees(&self) -> T;
}

impl<T: LazyReal> LatLon<T> for (T, T) {
    fn lat_deg(&self) -> T {
        self.0.clone()
    }

    fn lon_deg(&self) -> T {
        self.1.clone()
    }
}

impl<T: LazyReal> BearingDegrees<T> for T {
    fn degrees(&self) -> T {
        self.clone()
    }
}

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

/// Great-circle distance in metres between two latitude/longitude inputs.
#[must_use]
pub fn haversine_between_m<T, A, B>(from: &A, to: &B) -> T
where
    T: LazyReal,
    A: LatLon<T>,
    B: LatLon<T>,
{
    haversine_m(from.lat_deg(), from.lon_deg(), to.lat_deg(), to.lon_deg())
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

/// Initial great-circle compass bearing between two latitude/longitude inputs.
#[must_use]
pub fn initial_bearing_between_degrees<T, A, B>(from: &A, to: &B) -> T
where
    T: LazyReal,
    A: LatLon<T>,
    B: LatLon<T>,
{
    initial_bearing_degrees(from.lat_deg(), from.lon_deg(), to.lat_deg(), to.lon_deg())
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

/// Destination latitude/longitude after travelling from a latitude/longitude
/// input along a bearing input.
#[must_use]
pub fn destination_from_m<T, P, B>(origin: &P, bearing: &B, distance_m: T) -> (T, T)
where
    T: LazyReal,
    P: LatLon<T>,
    B: BearingDegrees<T>,
{
    destination_m(
        origin.lat_deg(),
        origin.lon_deg(),
        bearing.degrees(),
        distance_m,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual {actual} != expected {expected} within {tolerance}"
        );
    }

    #[test]
    fn haversine_scalar_matches_reference_distance() {
        let distance_m = haversine_m(37.7749_f64, -122.4194, 40.7128, -74.0060);

        assert_close(distance_m, 4_129_091.868_394_031, 1.0e-6);
    }

    #[test]
    fn haversine_runs_for_f32() {
        let distance_m = haversine_m(37.7749_f32, -122.4194, 40.7128, -74.0060);

        assert!(distance_m.is_finite());
    }

    #[test]
    fn lat_lon_trait_matches_scalar_haversine() {
        let from = (37.7749_f64, -122.4194_f64);
        let to = (40.7128_f64, -74.0060_f64);

        let actual: f64 = haversine_between_m(&from, &to);
        let expected = haversine_m(from.0, from.1, to.0, to.1);

        assert_close(actual, expected, 1.0e-12);
    }

    #[test]
    fn bearing_trait_matches_scalar_destination() {
        let origin = (37.7749_f64, -122.4194_f64);
        let bearing = 73.0_f64;
        let distance_m = 42_000.0_f64;

        let actual = destination_from_m(&origin, &bearing, distance_m);
        let expected = destination_m(origin.0, origin.1, bearing, distance_m);

        assert_close(actual.0, expected.0, 1.0e-12);
        assert_close(actual.1, expected.1, 1.0e-12);
    }
}
