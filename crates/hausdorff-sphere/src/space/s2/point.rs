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

    // --- S² as a statistical manifold ---------------------------------------
    //
    // The operations a probability distribution *over the sphere* is built from —
    // its support is the unit sphere, its metric is the geodesic angle, and a
    // tangent frame carries tangent-space constructions. Surfaced here so a
    // distribution (von Mises–Fisher and kin, in `bayes`) can be written against
    // S² as a space. This is the concrete S² instance of the `Sphere`/`Manifold`
    // abstraction the higher-dimensional spaces will later share.

    /// The inner product `p · q` — the **cosine of the geodesic (central) angle**
    /// between the two directions, and the core quantity of spherical
    /// distributions (the von Mises–Fisher density is `∝ exp(κ·μ·x)`).
    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        self.0.dot(other.0)
    }

    /// The point in the **direction of** `v` — its projection onto the sphere
    /// (`v.normalize()`). The boundary back onto S² after ℝ³ vector arithmetic
    /// (summing or averaging directions, stepping along a tangent). `v` must be
    /// non-zero.
    #[must_use]
    pub fn from_direction(v: R3Vector) -> Self {
        Self(v.normalize())
    }

    /// This point as its ambient **unit vector** in ℝ³ — the embedding
    /// coordinates, for vector arithmetic in the ambient/tangent space.
    #[must_use]
    pub fn as_vector(self) -> R3Vector {
        self.0
    }

    /// An orthonormal basis `(e₁, e₂)` of the **tangent plane** `T_pS²` at this
    /// point (both orthogonal to `self`); `self, e₁, e₂` form a right-handed
    /// frame. The frame a tangent-space construction is expressed in — e.g.
    /// sampling a uniform direction around `self`, or a wrapped/tangent Gaussian.
    #[must_use]
    pub fn tangent_basis(self) -> [R3Vector; 2] {
        let e1 = self.0.ortho();
        let e2 = self.0.cross(e1);
        [e1, e2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_is_cosine_of_geodesic_angle() {
        let p = S2Point::from_geo_point(GeoPoint::from_lat_lon(0.0, 0.0));
        let q = S2Point::from_geo_point(GeoPoint::from_lat_lon(0.0, 90.0));
        assert!((p.dot(p) - 1.0).abs() < 1e-12); // same point
        assert!(p.dot(q).abs() < 1e-12); // 90° apart → 0
        let anti = S2Point::from_geo_point(GeoPoint::from_lat_lon(0.0, 180.0));
        assert!((p.dot(anti) + 1.0).abs() < 1e-12); // antipodal → −1
    }

    #[test]
    fn from_direction_normalizes_and_round_trips() {
        let p = S2Point::from_geo_point(GeoPoint::from_lat_lon(37.0, -122.0));
        // A scaled copy of the underlying vector lands on the same point.
        let scaled = p.as_vector().mul(5.0);
        let q = S2Point::from_direction(scaled);
        assert!(q.as_vector().is_unit());
        assert!((p.dot(q) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn tangent_basis_is_orthonormal_and_perpendicular() {
        let p = S2Point::from_geo_point(GeoPoint::from_lat_lon(20.0, 50.0));
        let [e1, e2] = p.tangent_basis();
        let v = p.as_vector();
        assert!(e1.is_unit() && e2.is_unit());
        assert!(e1.dot(v).abs() < 1e-12); // e1 ⊥ p
        assert!(e2.dot(v).abs() < 1e-12); // e2 ⊥ p
        assert!(e1.dot(e2).abs() < 1e-12); // e1 ⊥ e2
    }
}
