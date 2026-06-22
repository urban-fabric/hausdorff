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
fn typed_sphere_matches_mean_earth_convenience() {
    let from = (37.7749_f64, -122.4194_f64);
    let to = (40.7128_f64, -74.0060_f64);
    let sphere = Sphere::<f64>::mean_earth();

    let actual = haversine_between_on_sphere_m(&from, &to, &sphere);
    let expected = haversine_between_m(&from, &to);

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

#[test]
fn typed_sphere_destination_matches_mean_earth_convenience() {
    let origin = (37.7749_f64, -122.4194_f64);
    let bearing = 73.0_f64;
    let distance_m = 42_000.0_f64;
    let sphere = Sphere::<f64>::mean_earth();

    let actual = destination_from_on_sphere_m(&origin, &bearing, distance_m, &sphere);
    let expected = destination_from_m(&origin, &bearing, distance_m);

    assert_close(actual.0, expected.0, 1.0e-12);
    assert_close(actual.1, expected.1, 1.0e-12);
}
