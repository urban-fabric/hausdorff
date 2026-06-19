//! Backend-neutral Earth kernels plus adapters for Hausdorff's eager records.

pub use hausdorff_kernels::earth::{
    BearingDegrees, LatLon, destination_from_m, destination_m, destination_on_m,
    haversine_between_m, haversine_diameter_m, haversine_m, haversine_on_m,
    initial_bearing_between_degrees, initial_bearing_degrees,
};

use super::{Bearing, GeoPoint};

impl LatLon<f64> for GeoPoint {
    fn lat_deg(&self) -> f64 {
        self.lat
    }

    fn lon_deg(&self) -> f64 {
        self.lon
    }
}

impl BearingDegrees<f64> for Bearing {
    fn degrees(&self) -> f64 {
        (*self).degrees()
    }
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
    fn haversine_trait_adapter_matches_concrete() {
        let from = GeoPoint::from_lat_lon(37.7749, -122.4194);
        let to = GeoPoint::from_lat_lon(40.7128, -74.0060);

        let actual: f64 = haversine_between_m(&from, &to);
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

    #[test]
    fn destination_trait_adapter_matches_concrete() {
        let origin = GeoPoint::from_lat_lon(37.7749, -122.4194);
        let bearing = Bearing::from_degrees(73.0);
        let distance_m = 42_000.0;

        let (actual_lat, actual_lon): (f64, f64) =
            destination_from_m(&origin, &bearing, distance_m);
        let expected = earth::destination_m(origin, bearing, distance_m);

        assert_close(actual_lat, expected.lat, 1.0e-12);
        assert_close(actual_lon, expected.lon, 1.0e-12);
    }
}
