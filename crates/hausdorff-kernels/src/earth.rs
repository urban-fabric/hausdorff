//! Earth and spherical geodesic kernels over any [`Real`] scalar.

use euclid_numeric::Real;

mod bbox;
mod proximity;
pub use bbox::{
    bbox_center_latlon, bbox_contains_bbox, bbox_contains_latlon, bbox_is_empty, bboxes_intersect,
};
pub use proximity::{
    haversine_between_within_radius_m, haversine_within_radius_m, haversine_within_radius_on_m,
    haversine_within_radius_on_sphere_m, normalize_bearing_degrees, normalize_longitude_degrees,
    valid_haversine_between_within_radius_m, valid_haversine_within_radius_m,
    valid_latitude_degrees, valid_latlon_degrees, valid_longitude_degrees,
};

/// Mean Earth radius in metres (WGS84 mean radius R1).
pub const EARTH_RADIUS_M: f64 = 6_371_008.8;

/// Geospatial constants that are not part of ordinary real analysis.
///
/// These are required methods so dataframe/tensor/expression backends can
/// represent domain constants natively instead of inheriting a default
/// `T::scalar(...)` lift in the generic kernels.
pub trait GeospatialConstants: Real {
    /// Mean Earth radius in metres (WGS84 mean radius R1).
    fn mean_earth_radius_m() -> Self;
}

impl GeospatialConstants for f64 {
    #[inline]
    fn mean_earth_radius_m() -> Self {
        EARTH_RADIUS_M
    }
}

impl GeospatialConstants for f32 {
    #[inline]
    fn mean_earth_radius_m() -> Self {
        EARTH_RADIUS_M as f32
    }
}

/// Backend-typed spherical Earth constants.
///
/// Build this once for a backend and reuse it across point batches instead of
/// lifting `EARTH_RADIUS_M` inside every distance/destination call.
#[derive(Clone, Debug, PartialEq)]
pub struct Sphere<T> {
    /// Sphere radius in metres.
    pub radius_m: T,
    /// Sphere diameter in metres.
    pub diameter_m: T,
}

impl<T: Real> Sphere<T> {
    /// Build spherical constants from a backend-typed radius.
    #[must_use]
    pub fn from_radius_m(radius_m: T) -> Self {
        Self {
            radius_m: radius_m.clone(),
            diameter_m: radius_m * T::two(),
        }
    }
}

impl<T: GeospatialConstants> Sphere<T> {
    /// Mean Earth sphere (WGS84 mean radius R1), using `T`'s native constant.
    #[must_use]
    pub fn mean_earth() -> Self {
        Self::from_radius_m(T::mean_earth_radius_m())
    }
}

/// A latitude/longitude input whose coordinates are decimal degrees.
///
/// Backend adapters can implement this for expression columns, tensor handles,
/// or small wrapper structs without depending on Hausdorff's concrete `GeoPoint`
/// storage contract.
pub trait LatLon<T: Real> {
    /// Latitude in decimal degrees.
    fn lat_deg(&self) -> T;

    /// Longitude in decimal degrees.
    fn lon_deg(&self) -> T;
}

/// A bearing input whose value is decimal degrees clockwise from north.
pub trait BearingDegrees<T: Real> {
    /// Bearing in decimal degrees.
    fn degrees(&self) -> T;
}

impl<T: Real> LatLon<T> for (T, T) {
    fn lat_deg(&self) -> T {
        self.0.clone()
    }

    fn lon_deg(&self) -> T {
        self.1.clone()
    }
}

impl<T: Real> BearingDegrees<T> for T {
    fn degrees(&self) -> T {
        self.clone()
    }
}

/// Great-circle distance in metres on a sphere with the supplied `diameter_m`.
///
/// Coordinate arguments are latitude/longitude in decimal degrees.
#[must_use]
pub fn haversine_diameter_m<T: Real>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
    diameter_m: T,
) -> T {
    let half = T::half();
    let zero = T::zero();
    let one = T::one();

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
pub fn haversine_on_m<T: Real>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
    radius_m: T,
) -> T {
    let two = T::two();
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
pub fn haversine_m<T: GeospatialConstants>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
) -> T {
    let diameter_m = T::mean_earth_radius_m() * T::two();
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
    T: GeospatialConstants,
    A: LatLon<T>,
    B: LatLon<T>,
{
    haversine_m(from.lat_deg(), from.lon_deg(), to.lat_deg(), to.lon_deg())
}

/// Great-circle distance between two latitude/longitude inputs on a typed
/// sphere, reusing pre-lifted spherical constants.
#[must_use]
pub fn haversine_between_on_sphere_m<T, A, B>(from: &A, to: &B, sphere: &Sphere<T>) -> T
where
    T: Real,
    A: LatLon<T>,
    B: LatLon<T>,
{
    haversine_diameter_m(
        from.lat_deg(),
        from.lon_deg(),
        to.lat_deg(),
        to.lon_deg(),
        sphere.diameter_m.clone(),
    )
}

/// Initial great-circle compass bearing in degrees, normalized into `[0, 360)`.
///
/// Coordinate arguments are latitude/longitude in decimal degrees.
#[must_use]
pub fn initial_bearing_degrees<T: Real>(
    from_lat_deg: T,
    from_lon_deg: T,
    to_lat_deg: T,
    to_lon_deg: T,
) -> T {
    let full_turn = T::full_turn_degrees();

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
    T: Real,
    A: LatLon<T>,
    B: LatLon<T>,
{
    initial_bearing_degrees(from.lat_deg(), from.lon_deg(), to.lat_deg(), to.lon_deg())
}

/// Destination latitude/longitude in decimal degrees after travelling
/// `distance_m` on a sphere of `radius_m` from an origin along `bearing_deg`.
#[must_use]
pub fn destination_on_m<T: Real>(
    origin_lat_deg: T,
    origin_lon_deg: T,
    bearing_deg: T,
    distance_m: T,
    radius_m: T,
) -> (T, T) {
    let minus_one = T::neg_one();
    let one = T::one();
    let lon_wrap_offset = T::three_half_turn_degrees();
    let full_turn = T::full_turn_degrees();
    let half_turn = T::half_turn_degrees();

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
pub fn destination_m<T: GeospatialConstants>(
    origin_lat_deg: T,
    origin_lon_deg: T,
    bearing_deg: T,
    distance_m: T,
) -> (T, T) {
    let radius_m = T::mean_earth_radius_m();
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
    T: GeospatialConstants,
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

/// Destination latitude/longitude after travelling from a latitude/longitude
/// input along a bearing input on a typed sphere.
#[must_use]
pub fn destination_from_on_sphere_m<T, P, B>(
    origin: &P,
    bearing: &B,
    distance_m: T,
    sphere: &Sphere<T>,
) -> (T, T)
where
    T: Real,
    P: LatLon<T>,
    B: BearingDegrees<T>,
{
    destination_on_m(
        origin.lat_deg(),
        origin.lon_deg(),
        bearing.degrees(),
        distance_m,
        sphere.radius_m.clone(),
    )
}

#[cfg(test)]
mod tests;
