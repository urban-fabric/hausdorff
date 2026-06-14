//! The S2 `Point` — a direction on the unit sphere — and conversions to and
//! from the workspace's WGS84 [`GeoPoint`]. Mirrors golang/geo's `s2.Point`
//! together with its `LatLng` conversions.

use hausdorff_base::earth::GeoPoint;
use hausdorff_base::space::s1::Angle;
use hausdorff_base::space::sphere::Sphere;

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

    /// The **great-circle (central-angle) distance** to `other`, in **radians**
    /// — the geodesic metric on the unit sphere, `acos(self · other)` computed
    /// robustly as `atan2(|p×q|, p·q)`. In `[0, π]`, zero iff the points
    /// coincide, `π` for antipodes; multiply by the radius for arc length.
    ///
    /// This is the distance that makes `S2Point` a genuine metric space (the
    /// [`MetricSpace`](hausdorff_planar::algebra::MetricSpace) instance via
    /// [`SphereMetric`](super::metric::SphereMetric)). Proof-backed by `riemann`
    /// `Hausdorff/SphereMetric`: `sphereDist` is `angle`, with the metric axioms
    /// `sphereDist_comm` (symmetry), `sphereDist_triangle` (triangle
    /// inequality), `sphereDist_nonneg`, `sphereDist_self` (self-distance zero),
    /// `sphereDist_le_pi`, and `cos_sphereDist_unit` (`cos dist = p·q`).
    #[must_use]
    pub fn distance(self, other: Self) -> f64 {
        self.0.angle(other.0).radians()
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

    // --- Great-circle navigation: the local ENU frame & direct/inverse problems

    /// The local **East–North–Up orthonormal frame** `(east, north, up)` at this
    /// point, the horizontal compass frame great-circle navigation is expressed
    /// in. With `p = self`:
    ///
    /// * `up = p` (the outward normal),
    /// * `east = normalize(ẑ × p) = normalize((−p.y, p.x, 0))`,
    /// * `north = p × east`.
    ///
    /// At a pole the East direction is undefined (`ẑ × p` vanishes); there this
    /// returns the degenerate frame with `east = north = 0` (and `up = p`), so
    /// callers can detect the singularity by a zero `east`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/Bearing`: this `east` is Lean's
    /// [`east`] `(−sin l, cos l, 0)` and this `north = p × east` is Lean's
    /// [`north`] `(−sin φ cos l, −sin φ sin l, cos φ)`, the frame proven
    /// orthonormal and horizontal by `east_unit`, `north_unit`, `east_perp_north`,
    /// `east_perp_pos`, and `north_perp_pos`.
    #[must_use]
    pub fn enu_frame(self) -> [R3Vector; 3] {
        let p = self.0;
        // ẑ × p = (−p.y, p.x, 0); zero exactly at the poles.
        let east = R3Vector::new(-p.y, p.x, 0.0).normalize();
        let north = p.cross(east);
        [east, north, p]
    }

    /// The **initial bearing** (forward azimuth) of the great-circle path from
    /// `self` to `other`, as an [`Angle`] measured **clockwise from North** in the
    /// local ENU frame at `self`. The direct counterpart is [`destination`]; this
    /// pair is the inverse/direct problem on unit directions (the `S2Point` layer,
    /// complementing [`Bearing::initial`] which works on `GeoPoint` at an Earth
    /// radius).
    ///
    /// The geodesic's initial heading is the tangent toward `other`,
    /// `t = other.vector().project_onto_tangent(self.vector())` (the proven
    /// `tangent p q = q − (p·q) p`); the bearing is `atan2(t·east, t·north)`.
    /// When `self` is a pole, or `other` coincides with / is antipodal to `self`
    /// (the frame or `t` is ~zero), the bearing is undefined and `Angle::ZERO`
    /// is returned.
    ///
    /// Proof-backed by `riemann`
    /// `Hausdorff/Bearing.initial_bearing_tangent` — `t·east` and `t·north` are
    /// *exactly* the two `atan2` arguments of the great-circle initial bearing —
    /// together with `dot_east`/`dot_north` (the frame coordinates) and
    /// `tangent_perp_pos` (`t` is genuinely tangent).
    ///
    /// [`destination`]: Self::destination
    /// [`Bearing::initial`]: hausdorff_base::earth::Bearing::initial
    #[must_use]
    pub fn bearing_to(self, other: Self) -> Angle {
        let [east, north, _up] = self.enu_frame();
        let t = other.0.project_onto_tangent(self.0);
        // Degenerate: pole (east == 0) or coincident/antipodal target (t == 0).
        if east.norm2() == 0.0 || t.norm2() == 0.0 {
            return Angle::ZERO;
        }
        Angle::from_radians(t.dot(east).atan2(t.dot(north)))
    }

    /// The **direct problem**: the point reached by travelling angular `distance`
    /// δ from `self` along the great circle with initial `bearing` θ (clockwise
    /// from North). The inverse counterpart is [`bearing_to`]; together they round
    /// trip, `p.destination(p.bearing_to(q), q-central-angle) ≈ q`.
    ///
    /// In the local ENU frame at `self`, the unit step direction is
    /// `t_dir = sin θ · east + cos θ · north`, and the result is
    /// `(cos δ · p + sin δ · t_dir).normalize()`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/Destination`: `dest = cos δ · p +
    /// sin δ · tdir` with `tdir = sin θ · east + cos θ · north` lands on the
    /// sphere (`dest_unit`) at exactly central angle δ from the origin
    /// (`dest_dot_origin`: `dest · p = cos δ`), the bearing tangent being a unit
    /// tangent vector (`tdir_unit`, `tdir_perp_p`).
    ///
    /// [`bearing_to`]: Self::bearing_to
    #[must_use]
    pub fn destination(self, bearing: Angle, distance: Angle) -> Self {
        let [east, north, _up] = self.enu_frame();
        // tdir = sin θ · east + cos θ · north  (unit, ⟂ self).
        let t_dir = east.mul(bearing.sin()).add(north.mul(bearing.cos()));
        // dest = cos δ · p + sin δ · tdir, renormalized to absorb f64 residue.
        let dest = self.0.mul(distance.cos()).add(t_dir.mul(distance.sin()));
        Self::from_vector(dest.normalize())
    }
}

/// The 2-sphere as the fixed `n = 3` instance of the general [`Sphere`] manifold, so a
/// distribution written over `Sphere` runs on geographic `S2Point` data unchanged.
impl Sphere for S2Point {
    fn ambient_dim(&self) -> usize {
        3
    }
    fn dot(&self, other: &Self) -> f64 {
        S2Point::dot(*self, *other)
    }
    fn write_ambient(&self, out: &mut [f64]) {
        out[0] = self.0.x;
        out[1] = self.0.y;
        out[2] = self.0.z;
    }
    fn from_ambient(coords: &[f64]) -> Self {
        Self::from_direction(R3Vector::new(coords[0], coords[1], coords[2]))
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
    fn distance_is_the_central_angle_metric() {
        use std::f64::consts::PI;
        let p = S2Point::from_geo_point(GeoPoint::from_lat_lon(0.0, 0.0));
        let q = S2Point::from_geo_point(GeoPoint::from_lat_lon(0.0, 90.0));
        let anti = S2Point::from_geo_point(GeoPoint::from_lat_lon(0.0, 180.0));
        // Self-distance is zero (sphereDist_self).
        assert!(p.distance(p).abs() < 1e-12);
        // 90° apart → π/2; antipodal → π (sphereDist_le_pi).
        assert!((p.distance(q) - PI / 2.0).abs() < 1e-12);
        assert!((p.distance(anti) - PI).abs() < 1e-12);
        // Symmetry (sphereDist_comm).
        assert!((p.distance(q) - q.distance(p)).abs() < 1e-12);
        // cos(dist) == p · q (cos_sphereDist_unit).
        assert!((p.distance(q).cos() - p.dot(q)).abs() < 1e-12);
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

    fn pt(x: f64, y: f64, z: f64) -> S2Point {
        S2Point::from_direction(R3Vector::new(x, y, z))
    }

    #[test]
    fn enu_frame_is_orthonormal_horizontal_and_matches_lean() {
        // At (1,0,0) (φ=0, l=0): east = (0,1,0), north = (0,0,1), up = (1,0,0).
        let p = pt(1.0, 0.0, 0.0);
        let [east, north, up] = p.enu_frame();
        assert!((east.x - 0.0).abs() < 1e-12 && (east.y - 1.0).abs() < 1e-12);
        assert!((north.z - 1.0).abs() < 1e-12);
        assert!((up.x - 1.0).abs() < 1e-12);
        // A generic point: orthonormal, horizontal, right-handed (north = p×east),
        // and east matches Lean's (−sin l, cos l, 0).
        let q = S2Point::from_geo_point(GeoPoint::from_lat_lon(20.0, 50.0));
        let [e, n, u] = q.enu_frame();
        assert!(e.is_unit() && n.is_unit() && u.is_unit());
        assert!(e.dot(n).abs() < 1e-12 && e.dot(u).abs() < 1e-12 && n.dot(u).abs() < 1e-12);
        let l = 50.0_f64.to_radians();
        assert!((e.x + l.sin()).abs() < 1e-12 && (e.y - l.cos()).abs() < 1e-12);
        // north == p × east (Lean `north = p × east`).
        let np = q.as_vector().cross(e);
        assert!(
            (n.x - np.x).abs() < 1e-12 && (n.y - np.y).abs() < 1e-12 && (n.z - np.z).abs() < 1e-12
        );
        // A pole has a degenerate (zero East/North) frame.
        let pole = pt(0.0, 0.0, 1.0);
        let [pe, pn, _pu] = pole.enu_frame();
        assert_eq!(pe.norm2(), 0.0);
        assert_eq!(pn.norm2(), 0.0);
    }

    #[test]
    fn bearing_to_cardinal_directions() {
        use std::f64::consts::FRAC_PI_2;
        let p = pt(1.0, 0.0, 0.0); // on the equator at lon 0
        // Due east: (0,1,0) → +90°.
        assert!((p.bearing_to(pt(0.0, 1.0, 0.0)).radians() - FRAC_PI_2).abs() < 1e-12);
        // Due north: (0,0,1) → 0°.
        assert!(p.bearing_to(pt(0.0, 0.0, 1.0)).radians().abs() < 1e-12);
        // Degenerate: coincident and pole origin return Angle::ZERO.
        assert_eq!(p.bearing_to(p), Angle::ZERO);
        let pole = pt(0.0, 0.0, 1.0);
        assert_eq!(pole.bearing_to(pt(1.0, 0.0, 0.0)), Angle::ZERO);
    }

    #[test]
    fn destination_steps_along_bearing() {
        use std::f64::consts::FRAC_PI_2;
        let p = pt(1.0, 0.0, 0.0);
        let delta = Angle::from_radians(0.5);
        // Due east → (cos δ, sin δ, 0).
        let east = p.destination(Angle::from_radians(FRAC_PI_2), delta);
        assert!((east.as_vector().x - delta.cos()).abs() < 1e-12);
        assert!((east.as_vector().y - delta.sin()).abs() < 1e-12);
        assert!(east.as_vector().z.abs() < 1e-12);
        // Due north → (cos δ, 0, sin δ).
        let north = p.destination(Angle::ZERO, delta);
        assert!((north.as_vector().x - delta.cos()).abs() < 1e-12);
        assert!(north.as_vector().y.abs() < 1e-12);
        assert!((north.as_vector().z - delta.sin()).abs() < 1e-12);
        // Destination sits at exactly central angle δ from the origin (dest_dot_origin).
        assert!((p.distance(east) - 0.5).abs() < 1e-12);
        assert!(east.as_vector().is_unit() && north.as_vector().is_unit()); // dest_unit
    }

    #[test]
    fn direct_undoes_inverse() {
        // p.destination(p.bearing_to(q), dist(p,q)) ≈ q for generic pairs.
        let cases = [
            (
                S2Point::from_geo_point(GeoPoint::from_lat_lon(37.0, -122.0)),
                S2Point::from_geo_point(GeoPoint::from_lat_lon(51.5, -0.13)),
            ),
            (
                S2Point::from_geo_point(GeoPoint::from_lat_lon(-33.9, 151.2)),
                S2Point::from_geo_point(GeoPoint::from_lat_lon(35.7, 139.7)),
            ),
        ];
        for (p, q) in cases {
            let brg = p.bearing_to(q);
            let dist = Angle::from_radians(p.distance(q));
            let back = p.destination(brg, dist);
            assert!(
                (back.dot(q) - 1.0).abs() < 1e-12,
                "round trip mismatch: dot = {}",
                back.dot(q)
            );
        }
    }
}
