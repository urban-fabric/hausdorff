//! The S2 `Point` — a direction on the unit sphere — and conversions to and
//! from the workspace's WGS84 [`GeoPoint`]. Mirrors golang/geo's `s2.Point`
//! together with its `LatLng` conversions.

use hausdorff_base::earth::GeoPoint;

use super::vector::R3Vector;

/// A point on the unit sphere, stored as a 3D direction vector. Shares
/// `hausdorff::earth::GeoPoint` at the lat/lon boundary so S2 and the rest of the
/// workspace speak the same coordinate type.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct S2Point(pub(super) R3Vector);

impl S2Point {
    pub(super) fn from_vector(v: R3Vector) -> Self {
        Self(v)
    }

    pub(super) fn vector(self) -> R3Vector {
        self.0
    }

    /// The point for a latitude/longitude given in **radians**.
    fn from_lat_lng_radians(lat: f64, lng: f64) -> Self {
        let cos_lat = lat.cos();
        Self(R3Vector::new(
            lng.cos() * cos_lat,
            lng.sin() * cos_lat,
            lat.sin(),
        ))
    }

    /// This point's (latitude, longitude) in **radians**.
    fn to_lat_lng_radians(self) -> (f64, f64) {
        let v = self.0;
        let lat = v.z.atan2((v.x * v.x + v.y * v.y).sqrt());
        let lng = v.y.atan2(v.x);
        (lat, lng)
    }

    /// The unit-sphere point for a WGS84 coordinate (decimal degrees).
    #[must_use]
    pub fn from_geo_point(p: GeoPoint) -> Self {
        Self::from_lat_lng_radians(p.lat.to_radians(), p.lon.to_radians())
    }

    /// This point as a WGS84 coordinate (decimal degrees).
    #[must_use]
    pub fn to_geo_point(self) -> GeoPoint {
        let (lat, lng) = self.to_lat_lng_radians();
        GeoPoint::from_lat_lon(lat.to_degrees(), lng.to_degrees())
    }
}
