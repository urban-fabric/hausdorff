//! A first-class **great circle** on the unit sphere — the oriented plane
//! through the origin, named by its unit normal (its pole).
//!
//! A [`GreatCircle`] unifies the great-circle operations that
//! [`arc`](super::arc) exposes as free functions over an `a → b` pair: a
//! point's signed [`cross_track`](GreatCircle::cross_track) distance, its
//! [`project`](GreatCircle::project)ion onto the circle, the in-plane
//! [`tangent`](GreatCircle::tangent_at) direction, and the circle's Clairaut
//! [`northernmost_latitude`](GreatCircle::northernmost_latitude). It is a
//! convenience wrapper, *not* a reimplementation: the geometry is the same
//! pole/cross-track/tangent math used in `arc.rs`, with the pole fixed once.

use hausdorff_base::space::s1::Angle;

use super::S2Point;
use super::vector::R3Vector;

/// An oriented **great circle**: the intersection of the unit sphere with the
/// plane through the origin whose **unit normal** is `normal` (the circle's
/// *pole*). The orientation is the normal's sign — travel around the circle is
/// counter-clockwise as seen from the pole side.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct GreatCircle {
    /// The unit pole — normal to the circle's plane.
    normal: R3Vector,
}

impl GreatCircle {
    /// The great circle through two points, oriented so travel `a → b` is
    /// counter-clockwise around the pole. The pole is `normalize(a × b)`;
    /// returns `None` when `a` and `b` are parallel or antipodal (the cross
    /// product is zero, so the connecting great circle is not unique).
    ///
    /// The same pole computation [`arc`](super::arc)'s functions use.
    ///
    /// Proof-backed by `riemann` `Hausdorff/Clairaut.arc_on_plane` (both `a` and
    /// `b` lie on the plane with normal `a × b`, so that plane carries the arc).
    #[must_use]
    pub fn from_points(a: S2Point, b: S2Point) -> Option<Self> {
        let n = a.vector().cross(b.vector());
        if n.norm2() > 0.0 {
            Some(Self {
                normal: n.normalize(),
            })
        } else {
            None
        }
    }

    /// The great circle with the given pole; `n` is normalized. (If `n` is the
    /// zero vector the normal stays zero and the circle is degenerate.)
    #[must_use]
    pub fn from_normal(n: R3Vector) -> Self {
        Self {
            normal: n.normalize(),
        }
    }

    /// The circle's **pole** — the unit normal to its plane.
    #[must_use]
    pub fn pole(&self) -> R3Vector {
        self.normal
    }

    /// The circle's unit **normal** (its [`pole`](Self::pole)).
    #[must_use]
    pub fn normal(&self) -> R3Vector {
        self.normal
    }

    /// Whether `p` lies **on** the great circle, within `tol`: its signed
    /// distance to the plane, `|p·normal|`, is `≤ tol`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/CrossTrack.foot_perp_normal` (a point
    /// is on the circle exactly when its component along the pole vanishes).
    #[must_use]
    pub fn contains(&self, p: S2Point, tol: f64) -> bool {
        p.vector().dot(self.normal).abs() <= tol
    }

    /// The **signed cross-track distance** of `p` to the circle — the angle
    /// `asin(clamp(p·normal, −1, 1))`, positive on the pole side. Its magnitude
    /// is the angular distance from `p` to the nearest point of the circle; the
    /// sign matches [`arc::cross_track`](super::arc::cross_track) for an arc whose
    /// pole is this `normal`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/CrossTrack`: `crossTrack_pythagoras`
    /// (the cross-/along-track decomposition) and `crossTrackSin_sq_le_one`
    /// (`sin²` of the cross-track is `≤ 1`, so the `asin` argument is valid —
    /// the explicit clamp only absorbs `f64` rounding).
    #[must_use]
    pub fn cross_track(&self, p: S2Point) -> Angle {
        Angle::from_radians(p.vector().dot(self.normal).clamp(-1.0, 1.0).asin())
    }

    /// The **nearest point on the circle** to `p` — `p` projected onto the
    /// circle's plane and renormalized to the sphere. Returns `None` when `p` is
    /// the pole (or antipole): then the projection is the zero vector and every
    /// point of the circle is equidistant.
    ///
    /// Proof-backed by `riemann` `Hausdorff/CrossTrack.foot_perp_normal` (the
    /// foot of the perpendicular onto the circle's plane is orthogonal to the
    /// pole, i.e. lies on the circle).
    #[must_use]
    pub fn project(&self, p: S2Point) -> Option<S2Point> {
        let foot = p.vector().project_onto_tangent(self.normal);
        // A non-degenerate foot has length ≈ cos(cross-track) ∈ (0, 1]; near the
        // pole it collapses to (numerical) zero. The tolerance absorbs the f64
        // residue of projecting an only-approximately-unit pole onto itself.
        if foot.norm2() > 1e-24 {
            Some(S2Point::from_vector(foot.normalize()))
        } else {
            None
        }
    }

    /// The in-plane **tangent direction** at `p` — `p`'s ambient vector with its
    /// pole component removed (`project_onto_tangent(normal)`), the direction of
    /// travel along the circle when `p` is on it. The zero vector when `p` is the
    /// pole. Not renormalized (callers that need a unit tangent normalize it).
    ///
    /// Proof-backed by `riemann` `Hausdorff/CrossTrack.foot_perp_normal` (the
    /// tangent is the foot of the perpendicular onto the plane, ⟂ the pole).
    #[must_use]
    pub fn tangent_at(&self, p: S2Point) -> R3Vector {
        p.vector().project_onto_tangent(self.normal)
    }

    /// The **northernmost latitude** the circle reaches, in radians — its
    /// Clairaut maximum `asin(√(1 − normal_z²))`. (The circle is symmetric about
    /// the equator, so its southernmost latitude is the negation.)
    ///
    /// Proof-backed by `riemann` `Hausdorff/Clairaut`: `latitude_bound` and
    /// `latitude_le_sqrt` (every point of the circle has `|sin(lat)| ≤ √(1 − n_z²)`,
    /// attained at the Clairaut vertex).
    #[must_use]
    pub fn northernmost_latitude(&self) -> f64 {
        (1.0 - self.normal.z * self.normal.z).max(0.0).sqrt().asin()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::earth::GeoPoint;

    fn sp(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    const FRAC_PI_2: f64 = std::f64::consts::FRAC_PI_2;

    /// The equator: the plane z = 0, pole = +ẑ.
    fn equator() -> GreatCircle {
        GreatCircle::from_normal(R3Vector::new(0.0, 0.0, 1.0))
    }

    #[test]
    fn equator_contains_equatorial_points_not_poles() {
        let eq = equator();
        // Points on the equator are contained.
        assert!(eq.contains(sp(0.0, 0.0), 1e-12));
        assert!(eq.contains(sp(0.0, 90.0), 1e-12));
        assert!(eq.contains(sp(0.0, -135.0), 1e-12));
        // The poles are not.
        assert!(!eq.contains(sp(90.0, 0.0), 1e-12));
        assert!(!eq.contains(sp(-90.0, 0.0), 1e-12));
        // An off-equator point is not.
        assert!(!eq.contains(sp(10.0, 20.0), 1e-12));
    }

    #[test]
    fn cross_track_of_pole_is_quarter_turn() {
        let eq = equator();
        // The north pole is +π/2 from the equator (pole side); south is −π/2.
        assert!((eq.cross_track(sp(90.0, 0.0)).radians() - FRAC_PI_2).abs() < 1e-12);
        assert!((eq.cross_track(sp(-90.0, 0.0)).radians() + FRAC_PI_2).abs() < 1e-12);
        // A point at 30°N is +30° (left of eastward travel along the equator).
        assert!((eq.cross_track(sp(30.0, 45.0)).degrees() - 30.0).abs() < 1e-12);
        // On the circle → 0.
        assert!(eq.cross_track(sp(0.0, 12.0)).radians().abs() < 1e-12);
    }

    #[test]
    fn cross_track_matches_arc_free_function() {
        // For an eastward equator arc the arc-pole is +ẑ, so the type and the
        // free function must agree.
        let a = sp(0.0, 0.0);
        let b = sp(0.0, 30.0);
        let gc = GreatCircle::from_points(a, b).unwrap();
        let p = sp(20.0, 10.0);
        let via_type = gc.cross_track(p).radians();
        let via_free = super::super::arc::cross_track(p, a, b).radians();
        assert!(
            (via_type - via_free).abs() < 1e-12,
            "{via_type} vs {via_free}"
        );
    }

    #[test]
    fn project_lands_on_circle_and_is_idempotent() {
        let gc = GreatCircle::from_normal(R3Vector::new(0.2, -0.5, 0.84).normalize());
        let p = sp(37.0, -122.0);
        let proj = gc.project(p).unwrap();
        // The projection lies on the circle.
        assert!(gc.contains(proj, 1e-12));
        // Projecting an already-on-circle point returns (essentially) itself.
        let proj2 = gc.project(proj).unwrap();
        assert!((proj.dot(proj2) - 1.0).abs() < 1e-12);
        // The pole has no projection.
        let pole = S2Point::from_vector(gc.pole());
        assert!(gc.project(pole).is_none());
        // The tangent at the pole is the zero vector.
        assert!(gc.tangent_at(pole).norm2() < 1e-24);
        // The tangent at an on-circle point is non-zero and ⟂ the pole.
        let t = gc.tangent_at(proj);
        assert!(t.norm2() > 1e-6);
        assert!(t.dot(gc.pole()).abs() < 1e-12);
    }

    #[test]
    fn northernmost_latitude_is_complement_of_pole_tilt() {
        // A pole at "latitude" α above the equatorial plane — n = (cos α, 0, sin α)
        // — gives a circle whose Clairaut max latitude is π/2 − α (the pole and
        // its circle are a quarter-turn apart). α = π/2 is the north pole → the
        // equator (max lat 0); α = 0 is an equatorial pole → a meridian (max lat π/2).
        for alpha_deg in [0.0, 10.0, 45.0, 90.0] {
            let alpha = alpha_deg * std::f64::consts::PI / 180.0;
            let n = R3Vector::new(alpha.cos(), 0.0, alpha.sin());
            let gc = GreatCircle::from_normal(n);
            let expected = FRAC_PI_2 - alpha;
            assert!(
                (gc.northernmost_latitude() - expected).abs() < 1e-12,
                "α={alpha_deg}: {} vs {expected}",
                gc.northernmost_latitude()
            );
        }
        // The equator (pole +ẑ) tops out at the equator (latitude 0).
        assert!(equator().northernmost_latitude().abs() < 1e-12);
        // A meridian circle (pole in the equatorial plane) reaches the pole.
        let meridian = GreatCircle::from_normal(R3Vector::new(0.0, 1.0, 0.0));
        assert!((meridian.northernmost_latitude() - FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn from_points_degenerate_is_none() {
        let a = sp(10.0, 20.0);
        // Same point → no unique circle.
        assert!(GreatCircle::from_points(a, a).is_none());
        // Antipodal → no unique circle.
        let anti = S2Point::from_vector(a.as_vector().mul(-1.0));
        assert!(GreatCircle::from_points(a, anti).is_none());
    }
}
