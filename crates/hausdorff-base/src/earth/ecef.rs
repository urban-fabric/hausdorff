//! Geocentric (**ECEF**) coordinates and local **ENU** frames on a reference
//! [`Ellipsoid`] — exact closed-form conversions between geodetic
//! (lat/lon/height) and Earth-Centred-Earth-Fixed Cartesian metres, and the
//! local east/north/up tangent frame anchored at an origin.
//!
//! A local ENU frame is often the cleanest "flat metres" workspace for
//! city-scale analysis: unlike UTM it has no zone seams and no scale factor at
//! the origin (distortion grows only quadratically with distance, ~d²/2R "up"
//! drop from curvature).

use crate::earth::GeoPoint;
use crate::earth::ellipsoid::Ellipsoid;

/// An Earth-Centred-Earth-Fixed Cartesian point, metres: `x` toward
/// (lat 0, lon 0), `y` toward (lat 0, lon 90°E), `z` toward the north pole.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EcefPoint {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// A displacement in a local [`EnuFrame`], metres east/north/up of the origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Enu {
    pub east: f64,
    pub north: f64,
    pub up: f64,
}

impl Ellipsoid {
    /// Geodetic → ECEF: the Cartesian position of `p` at `height_m` metres
    /// above this ellipsoid. Exact (closed form).
    #[must_use]
    pub fn geodetic_to_ecef(self, p: GeoPoint, height_m: f64) -> EcefPoint {
        let (slat, clat) = p.lat.to_radians().sin_cos();
        let (slon, clon) = p.lon.to_radians().sin_cos();
        // Prime-vertical radius of curvature (proof-backed; see `latitude`).
        let n = self.prime_vertical_radius(p.lat);
        EcefPoint {
            x: (n + height_m) * clat * clon,
            y: (n + height_m) * clat * slon,
            z: (n * (1.0 - self.e2()) + height_m) * slat,
        }
    }

    /// ECEF → geodetic: the latitude/longitude and ellipsoidal height of a
    /// Cartesian point. Bowring's method, iterated to convergence (2–3
    /// rounds; sub-nanometre for any terrestrial point).
    #[must_use]
    pub fn ecef_to_geodetic(self, p: EcefPoint) -> (GeoPoint, f64) {
        let e2 = self.e2();
        let ep2 = self.ep2();
        let b = self.b();
        let rho = p.x.hypot(p.y);
        let lon = p.y.atan2(p.x).to_degrees();

        if rho == 0.0 {
            // On the polar axis: latitude is ±90, height measured from the pole.
            let lat = if p.z >= 0.0 { 90.0 } else { -90.0 };
            return (GeoPoint::from_lat_lon(lat, lon), p.z.abs() - b);
        }

        // Bowring: iterate the parametric latitude β.
        let mut beta = (p.z * self.a).atan2(rho * b);
        let mut lat = 0.0;
        for _ in 0..4 {
            let (sb, cb) = beta.sin_cos();
            lat = (p.z + ep2 * b * sb * sb * sb).atan2(rho - e2 * self.a * cb * cb * cb);
            let new_beta = ((1.0 - self.f) * lat.tan()).atan();
            if (new_beta - beta).abs() < 1e-15 {
                break;
            }
            beta = new_beta;
        }
        let (slat, clat) = lat.sin_cos();
        // Robust at all latitudes (no division by cos φ).
        let height = rho * clat + p.z * slat - self.a * (1.0 - e2 * slat * slat).sqrt();
        (GeoPoint::from_lat_lon(lat.to_degrees(), lon), height)
    }

    /// True iff the ECEF point `p` lies on or inside this reference ellipsoid
    /// (ellipsoidal height ≤ 0) — the division-free test
    /// `b²(x² + y²) + a²z² ≤ a²b²` with semi-major `a` and semi-minor `b`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/Ellipsoid.inEllipsoid` (the solid
    /// spheroid is `inEllipsoid a a b`); the region is convex
    /// (`inEllipsoid_convex`). The 3-D mirror of `space::r2::Ellipse`.
    #[must_use]
    pub fn ecef_contains(self, p: EcefPoint) -> bool {
        let (a2, b2) = (self.a * self.a, self.b() * self.b());
        b2 * (p.x * p.x + p.y * p.y) + a2 * p.z * p.z <= a2 * b2
    }

    /// The local east/north/up tangent frame anchored at `origin` (at
    /// `height_m` above this ellipsoid).
    #[must_use]
    pub fn enu_frame(self, origin: GeoPoint, height_m: f64) -> EnuFrame {
        let (slat, clat) = origin.lat.to_radians().sin_cos();
        let (slon, clon) = origin.lon.to_radians().sin_cos();
        EnuFrame {
            ellipsoid: self,
            origin: self.geodetic_to_ecef(origin, height_m),
            east: [-slon, clon, 0.0],
            north: [-slat * clon, -slat * slon, clat],
            // The `up` axis is the geodetic surface normal (proof-backed).
            up: crate::earth::latitude::geodetic_normal(origin.lat, origin.lon),
        }
    }
}

/// A local east/north/up Cartesian frame anchored at a geodetic origin
/// ([`Ellipsoid::enu_frame`]): geodetic and ECEF positions become metre
/// displacements from the origin, exact in both directions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnuFrame {
    ellipsoid: Ellipsoid,
    origin: EcefPoint,
    east: [f64; 3],
    north: [f64; 3],
    up: [f64; 3],
}

impl EnuFrame {
    /// The ENU displacement of an ECEF point from this frame's origin.
    #[must_use]
    pub fn ecef_to_enu(&self, p: EcefPoint) -> Enu {
        let d = [
            self.origin.x - p.x,
            self.origin.y - p.y,
            self.origin.z - p.z,
        ];
        let d = [-d[0], -d[1], -d[2]];
        let dot = |r: [f64; 3]| r[0] * d[0] + r[1] * d[1] + r[2] * d[2];
        Enu {
            east: dot(self.east),
            north: dot(self.north),
            up: dot(self.up),
        }
    }

    /// The ECEF point at an ENU displacement from this frame's origin.
    #[must_use]
    pub fn enu_to_ecef(&self, enu: Enu) -> EcefPoint {
        let at =
            |i: usize| self.east[i] * enu.east + self.north[i] * enu.north + self.up[i] * enu.up;
        EcefPoint {
            x: self.origin.x + at(0),
            y: self.origin.y + at(1),
            z: self.origin.z + at(2),
        }
    }

    /// The ENU displacement of a geodetic point (at `height_m` above the
    /// frame's ellipsoid).
    #[must_use]
    pub fn to_enu(&self, p: GeoPoint, height_m: f64) -> Enu {
        self.ecef_to_enu(self.ellipsoid.geodetic_to_ecef(p, height_m))
    }

    /// The geodetic point (and ellipsoidal height) at an ENU displacement.
    #[must_use]
    pub fn from_enu(&self, enu: Enu) -> (GeoPoint, f64) {
        self.ellipsoid.ecef_to_geodetic(self.enu_to_ecef(enu))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(lat: f64, lon: f64) -> GeoPoint {
        GeoPoint::from_lat_lon(lat, lon)
    }

    #[test]
    fn ecef_anchors() {
        let e = Ellipsoid::WGS84;
        let p = e.geodetic_to_ecef(g(0.0, 0.0), 0.0);
        assert!((p.x - e.a).abs() < 1e-9 && p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
        let p = e.geodetic_to_ecef(g(0.0, 90.0), 0.0);
        assert!(p.x.abs() < 1e-9 && (p.y - e.a).abs() < 1e-9);
        let p = e.geodetic_to_ecef(g(90.0, 0.0), 0.0);
        assert!(p.x.abs() < 1e-9 && (p.z - e.b()).abs() < 1e-9);
    }

    #[test]
    fn ecef_contains_classifies_height_sign() {
        let e = Ellipsoid::WGS84;
        for &lat in &[-89.0, -23.5, 0.0, 45.0, 89.0] {
            for &lon in &[-120.0, 0.0, 121.5] {
                // Below the surface is inside; above is outside. (The surface
                // itself sits on the boundary to within float rounding, so it is
                // not asserted — `≤` can flip by an ulp there.)
                assert!(e.ecef_contains(e.geodetic_to_ecef(g(lat, lon), -1000.0)));
                assert!(!e.ecef_contains(e.geodetic_to_ecef(g(lat, lon), 1000.0)));
            }
        }
        // The geocentre is well inside.
        assert!(e.ecef_contains(EcefPoint { x: 0.0, y: 0.0, z: 0.0 }));
    }

    #[test]
    fn ecef_round_trips_everywhere() {
        let e = Ellipsoid::WGS84;
        for &lat in &[-89.9, -45.0, -0.1, 0.0, 23.5, 48.8566, 89.9, 90.0] {
            for &lon in &[-179.9, -120.0, 0.0, 2.3522, 121.5, 180.0] {
                for &h in &[-100.0, 0.0, 8848.0, 35_786_000.0] {
                    let (back, h_back) = e.ecef_to_geodetic(e.geodetic_to_ecef(g(lat, lon), h));
                    assert!((back.lat - lat).abs() < 1e-9, "lat {lat} → {}", back.lat);
                    if lat.abs() < 90.0 {
                        let dlon = (back.lon - lon).rem_euclid(360.0);
                        assert!(dlon.min(360.0 - dlon) < 1e-9, "lon {lon} → {}", back.lon);
                    }
                    assert!((h_back - h).abs() < 1e-6, "h {h} → {h_back}");
                }
            }
        }
    }

    #[test]
    fn enu_axes_point_the_right_way() {
        let e = Ellipsoid::WGS84;
        let frame = e.enu_frame(g(48.8566, 2.3522), 0.0);
        // A point slightly east: positive east, tiny north/up.
        let d = frame.to_enu(g(48.8566, 2.3622), 0.0);
        assert!(d.east > 500.0 && d.north.abs() < 5.0 && d.up < 0.0, "{d:?}");
        // Slightly north: positive north.
        let d = frame.to_enu(g(48.8666, 2.3522), 0.0);
        assert!(d.north > 1000.0 && d.east.abs() < 5.0, "{d:?}");
        // Straight up: positive up only.
        let d = frame.to_enu(g(48.8566, 2.3522), 100.0);
        assert!((d.up - 100.0).abs() < 1e-9 && d.east.abs() < 1e-9 && d.north.abs() < 1e-9);
    }

    #[test]
    fn enu_round_trip_and_short_baseline_distance() {
        let e = Ellipsoid::WGS84;
        let origin = g(48.8566, 2.3522);
        let frame = e.enu_frame(origin, 35.0);
        let target = g(48.8666, 2.3700);
        let enu = frame.to_enu(target, 35.0);
        let (back, h) = frame.from_enu(enu);
        assert!((back.lat - target.lat).abs() < 1e-12 && (back.lon - target.lon).abs() < 1e-12);
        assert!((h - 35.0).abs() < 1e-8);

        // For a ~1.7 km baseline *on the surface*, the ENU chord length
        // matches the geodesic distance to well under a millimetre (chord-vs-
        // arc error is cubic in distance). At altitude the chord is longer by
        // h/R — that's physics, not error, so compare at h = 0.
        let surface = e.enu_frame(origin, 0.0).to_enu(target, 0.0);
        let chord =
            (surface.east * surface.east + surface.north * surface.north + surface.up * surface.up)
                .sqrt();
        let geodesic = e.distance_m(origin, target);
        assert!((chord - geodesic).abs() < 1e-3, "{chord} vs {geodesic}");
    }
}
