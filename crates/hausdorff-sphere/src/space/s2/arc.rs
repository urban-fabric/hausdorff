//! Great-circle **arc** kernels on the unit sphere: interpolation along an
//! arc, cross-track / along-track decomposition of a point against an arc,
//! the closest point on an arc segment, and the latitude extremes of an arc
//! (Clairaut's relation). All results are angles ([`Angle`], radians) or unit
//! points — radius-free, like the rest of `s2`. The lat/lon (and metres)
//! layer over these lives in [`hausdorff_base::earth::arc`].

use std::f64::consts::PI;

use hausdorff_base::space::s1::Angle;

use super::vector::R3Vector;
use super::{Rect, S2Point};

/// The point a fraction `t` of the way along the great-circle arc `a → b`
/// (spherical linear interpolation): `t = 0` is `a`, `t = 1` is `b`; values
/// outside `[0, 1]` extrapolate along the same great circle.
///
/// For (numerically) antipodal endpoints the connecting great circle is not
/// unique; an arbitrary one through both points is used.
#[must_use]
pub fn interpolate(a: S2Point, b: S2Point, t: f64) -> S2Point {
    let va = a.vector();
    let theta = va.angle(b.vector()).radians();
    // Rotate `a` by t·θ about the arc's pole (Rodrigues, with unit axis ⊥ a):
    // p = a·cos(tθ) + (n̂ × a)·sin(tθ).
    let axis = arc_pole(a, b);
    let (s, c) = (t * theta).sin_cos();
    S2Point::from_vector(va.mul(c).add(axis.cross(va).mul(s)).normalize())
}

/// `S2Point` interpolates along the **great circle** (SLERP) — the spherical
/// instance of the shared
/// [`Interpolate`](hausdorff_planar::algebra::Interpolate) trait. The endpoint
/// laws `lerp(a, b, 0) = a`, `lerp(a, b, 1) = b` are the riemann `Hausdorff/Slerp`
/// theorems (`slerp_zero` / `slerp_angle`), so generic `Interpolate` code
/// (`midpoint`, `sample`) runs on the sphere with no special case.
impl hausdorff_planar::algebra::Interpolate for S2Point {
    type Scalar = f64;

    #[inline]
    fn lerp(self, q: Self, t: f64) -> Self {
        interpolate(self, q, t)
    }
}

/// The midpoint of the arc `a → b`.
#[must_use]
pub fn midpoint(a: S2Point, b: S2Point) -> S2Point {
    interpolate(a, b, 0.5)
}

/// The pole (unit normal) of the great circle through `a` and `b`, oriented so
/// travel `a → b` is counter-clockwise around it — i.e. the pole is on the
/// **left** of the direction of travel. Falls back to an arbitrary
/// perpendicular axis when `a` and `b` are parallel or antiparallel.
fn arc_pole(a: S2Point, b: S2Point) -> R3Vector {
    let n = a.vector().cross(b.vector());
    if n.norm2() > 0.0 {
        n.normalize()
    } else {
        a.vector().ortho()
    }
}

/// Signed cross-track angle of `p` against the great circle through `a → b`:
/// positive when `p` is **left** of the direction of travel, negative right.
/// The magnitude is the angular distance from `p` to the circle (not clamped
/// to the segment — see [`distance_to_arc`]).
#[must_use]
pub fn cross_track(p: S2Point, a: S2Point, b: S2Point) -> Angle {
    let n = arc_pole(a, b);
    Angle::from_radians(p.vector().dot(n).clamp(-1.0, 1.0).asin())
}

/// Along-track angle of `p` against the arc `a → b`: the signed angular
/// distance from `a` to the foot of `p`'s perpendicular on the great circle,
/// positive toward `b`. Together with [`cross_track`] this decomposes `p`
/// relative to the arc.
#[must_use]
pub fn along_track(p: S2Point, a: S2Point, b: S2Point) -> Angle {
    let n = arc_pole(a, b);
    // Project p onto the circle's plane; measure the signed angle from a.
    let vp = p.vector();
    let foot = vp.project_onto_tangent(n);
    if foot.norm2() == 0.0 {
        // p is the circle's pole: every foot is equidistant; call it 0.
        return Angle::from_radians(0.0);
    }
    let foot = foot.normalize();
    let va = a.vector();
    let sign = if n.dot(va.cross(foot)) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Angle::from_radians(sign * va.angle(foot).radians())
}

/// The closest point to `p` on the arc **segment** `a → b` (the perpendicular
/// foot when it falls within the segment, else the nearer endpoint).
#[must_use]
pub fn closest_point(p: S2Point, a: S2Point, b: S2Point) -> S2Point {
    let arc_len = a.vector().angle(b.vector()).radians();
    let along = along_track(p, a, b).radians();
    if along <= 0.0 {
        a
    } else if along >= arc_len {
        b
    } else {
        interpolate(a, b, along / arc_len)
    }
}

/// Angular distance from `p` to the arc **segment** `a → b` — the distance to
/// [`closest_point`], so endpoints clamp (unlike the unclamped
/// [`cross_track`] magnitude).
#[must_use]
pub fn distance_to_arc(p: S2Point, a: S2Point, b: S2Point) -> Angle {
    p.vector().angle(closest_point(p, a, b).vector())
}

/// A conservative lat/lng rectangle containing the minor arc segment `a -> b`.
///
/// Latitude uses [`latitude_extremes`] (Clairaut's relation) so a high-riding
/// great-circle arc does not escape the endpoint latitude band. Longitude is
/// monotone along a non-pole-touching minor arc; if the arc touches a pole, the
/// rectangle widens to every longitude.
///
/// Proof-backed by `riemann` `Hausdorff/Clairaut.latitude_bound` for the
/// latitude band and by `Hausdorff/SphereLongitude`/`SphereRect` for the
/// wrap-aware longitude interval.
#[must_use]
pub fn rect_bound(a: S2Point, b: S2Point) -> Rect {
    const PAD_DEG: f64 = 1e-9;
    const POLE_EPS_DEG: f64 = 1e-12;

    if a.vector().angle(b.vector()).radians() >= PI - 1e-15 {
        return Rect::full();
    }

    let ag = a.to_geo_point();
    let bg = b.to_geo_point();
    let (min_lat, max_lat) = latitude_extremes(a, b);
    let min_lat = (min_lat.degrees() - PAD_DEG).max(-90.0);
    let max_lat = (max_lat.degrees() + PAD_DEG).min(90.0);
    if max_lat >= 90.0 - POLE_EPS_DEG || min_lat <= -90.0 + POLE_EPS_DEG {
        return Rect::from_degrees(min_lat, -180.0, max_lat, 180.0);
    }

    let delta = (bg.lon - ag.lon + 540.0).rem_euclid(360.0) - 180.0;
    let (lng_lo, lng_hi) = if delta >= 0.0 {
        (ag.lon - PAD_DEG, ag.lon + delta + PAD_DEG)
    } else {
        (ag.lon + delta - PAD_DEG, ag.lon + PAD_DEG)
    };
    Rect::from_degrees(
        min_lat,
        normalize_degrees(lng_lo),
        max_lat,
        normalize_degrees(lng_hi),
    )
}

fn normalize_degrees(degrees: f64) -> f64 {
    let x = (degrees + 180.0).rem_euclid(360.0) - 180.0;
    if x <= -180.0 { 180.0 } else { x }
}

/// The minimum and maximum **latitude** attained on the arc segment `a → b`,
/// in radians. By Clairaut's relation the great circle's extreme latitude is
/// `±asin(√(1 − n_z²))` (with `n` the circle's pole); the extreme counts only
/// if its meridian-tangent point falls within the segment, otherwise the
/// endpoint latitudes bound the arc.
#[must_use]
pub fn latitude_extremes(a: S2Point, b: S2Point) -> (Angle, Angle) {
    let (va, vb) = (a.vector(), b.vector());
    let lat = |v: R3Vector| v.z.atan2((v.x * v.x + v.y * v.y).sqrt());
    let mut min_lat = lat(va).min(lat(vb));
    let mut max_lat = lat(va).max(lat(vb));

    let n = arc_pole(a, b);
    // The northernmost point of the full circle: the unit vector in the
    // circle's plane with the largest z — normalize(ẑ − n·n_z), i.e. ẑ
    // projected onto the circle's tangent plane.
    let vertex = R3Vector::new(0.0, 0.0, 1.0).project_onto_tangent(n);
    if vertex.norm2() > 0.0 {
        let north = vertex.normalize();
        if on_segment(north, a, b, n) {
            max_lat = max_lat.max(lat(north));
        }
        let south = north.mul(-1.0);
        if on_segment(south, a, b, n) {
            min_lat = min_lat.min(lat(south));
        }
    }
    (Angle::from_radians(min_lat), Angle::from_radians(max_lat))
}

/// True iff `v` (a point on the great circle with pole `n`) lies on the minor
/// arc from `a` to `b`: it must be reached rotating from `a` toward `b` and
/// not past `b`.
fn on_segment(v: R3Vector, a: S2Point, b: S2Point, n: R3Vector) -> bool {
    n.dot(a.vector().cross(v)) >= 0.0 && n.dot(v.cross(b.vector())) >= 0.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::earth::GeoPoint;

    fn sp(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    #[test]
    fn interpolate_endpoints_and_midpoint() {
        let a = sp(0.0, 0.0);
        let b = sp(0.0, 90.0);
        let mid = midpoint(a, b).to_geo_point();
        assert!(mid.lat.abs() < 1e-12 && (mid.lon - 45.0).abs() < 1e-12);
        let q = interpolate(a, b, 1.0).to_geo_point();
        assert!((q.lon - 90.0).abs() < 1e-12);
    }

    #[test]
    fn interpolate_trait_impl_walks_the_great_circle() {
        use hausdorff_planar::algebra::{Interpolate, sample};
        let a = sp(0.0, 0.0);
        let b = sp(0.0, 90.0);
        // Endpoint laws + great-circle midpoint through the shared trait.
        assert!(a.lerp(b, 0.0).to_geo_point().lon.abs() < 1e-12);
        assert!((a.lerp(b, 1.0).to_geo_point().lon - 90.0).abs() < 1e-12);
        let mid = a.midpoint(b).to_geo_point();
        assert!(mid.lat.abs() < 1e-12 && (mid.lon - 45.0).abs() < 1e-12);
        // The generic `sample` walks the great circle (lon 0, 30, 60, 90).
        let pts = sample(a, b, 3);
        assert_eq!(pts.len(), 4);
        assert!((pts[1].to_geo_point().lon - 30.0).abs() < 1e-12);
        assert!((pts[2].to_geo_point().lon - 60.0).abs() < 1e-12);
    }

    #[test]
    fn cross_track_sign_and_magnitude() {
        // Arc west→east along the equator; a point at 1°N is 1° to the LEFT.
        let a = sp(0.0, 0.0);
        let b = sp(0.0, 10.0);
        let xt = cross_track(sp(1.0, 5.0), a, b);
        assert!((xt.degrees() - 1.0).abs() < 1e-12, "{}", xt.degrees());
        let xt_south = cross_track(sp(-2.0, 5.0), a, b);
        assert!((xt_south.degrees() + 2.0).abs() < 1e-12);
    }

    #[test]
    fn along_track_and_closest_point_clamp() {
        let a = sp(0.0, 0.0);
        let b = sp(0.0, 10.0);
        // Abeam at lon 4: foot at (0, 4).
        let at = along_track(sp(3.0, 4.0), a, b);
        assert!((at.degrees() - 4.0).abs() < 1e-2, "{}", at.degrees());
        // Behind a: negative along-track, clamps to a.
        assert!(along_track(sp(0.0, -5.0), a, b).degrees() < 0.0);
        let c = closest_point(sp(1.0, -5.0), a, b).to_geo_point();
        assert!(c.lat.abs() < 1e-12 && c.lon.abs() < 1e-12);
        // Beyond b clamps to b.
        let c = closest_point(sp(1.0, 15.0), a, b).to_geo_point();
        assert!((c.lon - 10.0).abs() < 1e-12);
    }

    #[test]
    fn latitude_extremes_finds_clairaut_vertex() {
        // A 90°-of-longitude arc at lat 45°: its midpoint rides higher than
        // either endpoint (the Clairaut vertex).
        let a = sp(45.0, 0.0);
        let b = sp(45.0, 90.0);
        let (min, max) = latitude_extremes(a, b);
        assert!((min.degrees() - 45.0).abs() < 1e-9, "{}", min.degrees());
        assert!(max.degrees() > 45.0 + 5.0, "{}", max.degrees());
        // The vertex latitude must match the sampled maximum.
        let sampled = (0..=1000)
            .map(|i| interpolate(a, b, i as f64 / 1000.0).to_geo_point().lat)
            .fold(f64::MIN, f64::max);
        assert!(
            (max.degrees() - sampled).abs() < 1e-4,
            "{} vs {sampled}",
            max.degrees()
        );
    }

    #[test]
    fn rect_bound_contains_high_latitude_arc_samples() {
        let a = sp(45.0, 0.0);
        let b = sp(45.0, 90.0);
        let bound = rect_bound(a, b);
        assert!(
            bound.hi().lat > 50.0,
            "bound should include the Clairaut vertex"
        );
        for i in 0..=100 {
            let p = interpolate(a, b, f64::from(i) / 100.0).to_geo_point();
            assert!(bound.contains_lat_lng(p), "{p:?} escaped {bound:?}");
        }
    }

    #[test]
    fn rect_bound_wraps_antimeridian() {
        let a = sp(10.0, 170.0);
        let b = sp(10.0, -170.0);
        let bound = rect_bound(a, b);
        assert!(bound.contains_lat_lng(a.to_geo_point()));
        assert!(bound.contains_lat_lng(b.to_geo_point()));
        for i in 0..=100 {
            let p = interpolate(a, b, f64::from(i) / 100.0).to_geo_point();
            assert!(bound.contains_lat_lng(p), "{p:?} escaped {bound:?}");
        }
    }
}
