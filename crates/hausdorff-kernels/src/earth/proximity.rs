//! Latitude/longitude validity and proximity kernels over scalar columns.

use euclid::numeric::{Boolean as BooleanMask, Compare, Real, RealConstants};

use super::{
    GeospatialConstants, LatLon, Sphere, haversine_diameter_m, haversine_m, haversine_on_m,
};

/// The Boolean mask produced by comparisons on scalar `T`.
pub type Boolean<T> = <T as Compare>::Mask;

#[inline]
fn mask_and<M: BooleanMask>(a: M, b: M) -> M {
    a.and(&b)
}

/// Mask that is true where `lat_deg` is finite and in `[-90, 90]`.
#[inline]
#[must_use]
pub fn valid_latitude_degrees<T: Real + RealConstants>(lat_deg: T) -> Boolean<T> {
    let max = T::quarter_turn_degrees();
    let min = -max.clone();
    let finite = lat_deg.clone().is_finite();
    let in_range = mask_and(min.le(&lat_deg), lat_deg.le(&max));

    mask_and(finite, in_range)
}

/// Mask that is true where `lon_deg` is finite and in `[-180, 180]`.
#[inline]
#[must_use]
pub fn valid_longitude_degrees<T: Real + RealConstants>(lon_deg: T) -> Boolean<T> {
    let max = T::half_turn_degrees();
    let min = -max.clone();
    let finite = lon_deg.clone().is_finite();
    let in_range = mask_and(min.le(&lon_deg), lon_deg.le(&max));

    mask_and(finite, in_range)
}

/// Mask that is true where `(lat_deg, lon_deg)` is a finite, in-range WGS84
/// latitude/longitude coordinate.
#[inline]
#[must_use]
pub fn valid_latlon_degrees<T: Real + RealConstants>(lat_deg: T, lon_deg: T) -> Boolean<T> {
    mask_and(
        valid_latitude_degrees(lat_deg),
        valid_longitude_degrees(lon_deg),
    )
}

/// Normalize longitude degrees into `[-180, 180)`.
#[inline]
#[must_use]
pub fn normalize_longitude_degrees<T: Real + RealConstants>(lon_deg: T) -> T {
    (lon_deg + T::three_half_turn_degrees()).rem_euclid(T::full_turn_degrees())
        - T::half_turn_degrees()
}

/// Normalize compass bearing degrees into `[0, 360)`.
#[inline]
#[must_use]
pub fn normalize_bearing_degrees<T: Real + RealConstants>(bearing_deg: T) -> T {
    bearing_deg.rem_euclid(T::full_turn_degrees())
}

/// Mask that is true where two coordinates are within `radius_m` on a sphere of
/// `sphere_radius_m`.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn haversine_within_radius_on_m<T: Real>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
    radius_m: T,
    sphere_radius_m: T,
) -> Boolean<T> {
    let distance_m = haversine_on_m(
        from_lat_deg,
        from_lon_deg,
        to_lat_deg,
        to_lon_deg,
        sphere_radius_m,
    );

    distance_m.le(&radius_m)
}

/// Mask that is true where two coordinates are within `radius_m` on a typed
/// sphere, reusing pre-lifted spherical constants.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn haversine_within_radius_on_sphere_m<T: Real>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
    radius_m: T,
    sphere: &Sphere<T>,
) -> Boolean<T> {
    let distance_m = haversine_diameter_m(
        from_lat_deg,
        from_lon_deg,
        to_lat_deg,
        to_lon_deg,
        sphere.diameter_m.clone(),
    );

    distance_m.le(&radius_m)
}

/// Mask that is true where two coordinates are within `radius_m` on Earth.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn haversine_within_radius_m<T: GeospatialConstants>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
    radius_m: T,
) -> Boolean<T> {
    let distance_m = haversine_m(from_lat_deg, from_lon_deg, to_lat_deg, to_lon_deg);

    distance_m.le(&radius_m)
}

/// Mask that is true where two latitude/longitude inputs are within `radius_m`
/// on Earth.
#[inline]
#[must_use]
pub fn haversine_between_within_radius_m<T, A, B>(from: &A, to: &B, radius_m: T) -> Boolean<T>
where
    T: GeospatialConstants,
    A: LatLon<T>,
    B: LatLon<T>,
{
    haversine_within_radius_m(
        from.lat_deg(),
        from.lon_deg(),
        to.lat_deg(),
        to.lon_deg(),
        radius_m,
    )
}

/// Mask that is true where both coordinates are valid WGS84 lon/lat points and
/// the Earth haversine distance is within `radius_m`.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn valid_haversine_within_radius_m<T: GeospatialConstants + RealConstants>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
    radius_m: T,
) -> Boolean<T> {
    let valid_from = valid_latlon_degrees(from_lat_deg.clone(), from_lon_deg.clone());
    let valid_to = valid_latlon_degrees(to_lat_deg.clone(), to_lon_deg.clone());
    let within =
        haversine_within_radius_m(from_lat_deg, from_lon_deg, to_lat_deg, to_lon_deg, radius_m);

    mask_and(mask_and(valid_from, valid_to), within)
}

/// Mask that is true where two latitude/longitude inputs are valid WGS84
/// lon/lat points and within `radius_m` on Earth.
#[inline]
#[must_use]
pub fn valid_haversine_between_within_radius_m<T, A, B>(from: &A, to: &B, radius_m: T) -> Boolean<T>
where
    T: GeospatialConstants + RealConstants,
    A: LatLon<T>,
    B: LatLon<T>,
{
    valid_haversine_within_radius_m(
        from.lat_deg(),
        from.lon_deg(),
        to.lat_deg(),
        to.lon_deg(),
        radius_m,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn latlon_validity_is_mask_valued() {
        assert!(valid_latlon_degrees::<f64>(90.0, -180.0));
        assert!(valid_latlon_degrees::<f64>(-90.0, 180.0));
        assert!(!valid_latlon_degrees::<f64>(90.0001, 0.0));
        assert!(!valid_latlon_degrees::<f64>(0.0, 180.0001));
        assert!(!valid_latlon_degrees::<f64>(f64::NAN, 0.0));
        assert!(!valid_latlon_degrees::<f64>(0.0, f64::INFINITY));
    }

    #[test]
    fn normalizers_use_backend_constants() {
        assert_eq!(normalize_longitude_degrees(181.0_f64), -179.0);
        assert_eq!(normalize_longitude_degrees(-181.0_f64), 179.0);
        assert_eq!(normalize_bearing_degrees(-10.0_f64), 350.0);
        assert_eq!(normalize_bearing_degrees(370.0_f64), 10.0);
    }

    #[test]
    fn within_radius_matches_haversine_distance() {
        let sf = (37.7749_f64, -122.4194_f64);
        let oakland = (37.8044_f64, -122.2712_f64);
        let nyc = (40.7128_f64, -74.0060_f64);

        assert!(haversine_within_radius_m(
            sf.0, sf.1, oakland.0, oakland.1, 20_000.0
        ));
        assert!(!haversine_within_radius_m(
            sf.0, sf.1, nyc.0, nyc.1, 20_000.0
        ));
        assert!(haversine_between_within_radius_m(&sf, &oakland, 20_000.0));
        assert!(!haversine_between_within_radius_m(&sf, &nyc, 20_000.0));

        let sphere = Sphere::<f64>::mean_earth();
        assert_eq!(
            haversine_within_radius_on_sphere_m(
                sf.0, sf.1, oakland.0, oakland.1, 20_000.0, &sphere
            ),
            haversine_within_radius_m(sf.0, sf.1, oakland.0, oakland.1, 20_000.0)
        );
    }

    #[test]
    fn valid_within_radius_rejects_invalid_coordinates() {
        assert!(!valid_haversine_within_radius_m::<f64>(
            91.0, 0.0, 0.0, 0.0, 20_000.0
        ));
        assert!(!valid_haversine_within_radius_m::<f64>(
            0.0,
            0.0,
            0.0,
            f64::NAN,
            20_000.0
        ));
        assert!(!valid_haversine_between_within_radius_m(
            &(91.0_f64, 0.0),
            &(0.0_f64, 0.0),
            20_000.0
        ));
    }

    #[test]
    fn proximity_kernels_run_for_f32() {
        assert!(valid_latlon_degrees::<f32>(45.0, -122.0));
        assert!(haversine_within_radius_m::<f32>(
            37.7749, -122.4194, 37.8044, -122.2712, 20_000.0
        ));
    }
}
