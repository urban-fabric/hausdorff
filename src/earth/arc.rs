//! Great-circle **arcs** between lat/lon points — the Earth-facing layer over
//! the unit-sphere kernels in [`crate::space::s2::arc`] and the edge
//! predicates in [`crate::space::s2::edges`]. Everything here is exact on the
//! sphere (closed forms, no iteration); for the ellipsoidal counterparts see
//! [`Ellipsoid`](crate::earth::Ellipsoid).
//!
//! Distances come back in metres via [`Sphere`](crate::earth::sphere::Sphere)
//! methods (radius-aware: `Sphere::EARTH` or any other planet); the
//! point-valued operations (interpolation, intersection, closest point) are
//! radius-free and live as free functions.

use crate::earth::sphere::Sphere;
use crate::earth::{BBox, GeoPoint};
use crate::space::s2::arc as s2arc;
use crate::space::s2::{Crossing, Loop, S2Point, crossing_sign, intersection as s2_intersection};

fn sp(p: GeoPoint) -> S2Point {
    S2Point::from_geo_point(p)
}

/// The point a fraction `t` along the great-circle arc `a → b` (`t = 0` is
/// `a`, `t = 1` is `b`; values outside `[0, 1]` extrapolate). For antipodal
/// endpoints the circle is not unique and an arbitrary one is used.
#[must_use]
pub fn interpolate(a: GeoPoint, b: GeoPoint, t: f64) -> GeoPoint {
    s2arc::interpolate(sp(a), sp(b), t).to_geo_point()
}

/// The midpoint of the great-circle arc `a → b`.
#[must_use]
pub fn midpoint(a: GeoPoint, b: GeoPoint) -> GeoPoint {
    s2arc::midpoint(sp(a), sp(b)).to_geo_point()
}

/// The closest point to `p` on the arc **segment** `a → b` (the perpendicular
/// foot, clamped to the nearer endpoint when it falls outside the segment).
#[must_use]
pub fn closest_point_on_arc(p: GeoPoint, a: GeoPoint, b: GeoPoint) -> GeoPoint {
    s2arc::closest_point(sp(p), sp(a), sp(b)).to_geo_point()
}

/// Where the arc segments `a0 → a1` and `b0 → b1` cross, if they cross at a
/// point interior to both. Touching at a shared endpoint (or any other
/// boundary-degenerate contact) returns `None` — test those cases explicitly
/// with [`crate::space::s2::crossing_sign`] if they matter.
#[must_use]
pub fn intersection(a0: GeoPoint, a1: GeoPoint, b0: GeoPoint, b1: GeoPoint) -> Option<GeoPoint> {
    let (a0, a1, b0, b1) = (sp(a0), sp(a1), sp(b0), sp(b1));
    match crossing_sign(a0, a1, b0, b1) {
        Crossing::Cross => Some(s2_intersection(a0, a1, b0, b1).to_geo_point()),
        Crossing::Maybe | Crossing::DoNotCross => None,
    }
}

/// True iff `p` lies inside the spherical polygon whose boundary is `ring`
/// (lat/lon vertices, closed or unclosed; consecutive vertices joined by
/// great-circle arcs). The interior is taken to be the **smaller** side of the
/// ring regardless of winding — a planar-`r2` polygon ring straddling the
/// antimeridian or a pole tests correctly here where a lon/lat point-in-
/// polygon would not. Rings with fewer than three distinct vertices contain
/// nothing.
#[must_use]
pub fn ring_contains(ring: &[GeoPoint], p: GeoPoint) -> bool {
    let mut pts: Vec<S2Point> = ring.iter().map(|&v| sp(v)).collect();
    if pts.len() >= 2 && ring.first() == ring.last() {
        pts.pop();
    }
    if pts.len() < 3 {
        return false;
    }
    let l = Loop::new(pts);
    // Normalize to the smaller side: a loop enclosing more than half the
    // sphere is the same boundary wound the other way.
    let l = if l.area() > 2.0 * std::f64::consts::PI {
        let mut v = l.vertices().to_vec();
        v.reverse();
        Loop::new(v)
    } else {
        l
    };
    l.contains_point(sp(p))
}

/// The lat/lon bounding box of the arc **segment** `a → b`. The latitude
/// bound includes the Clairaut vertex when the arc rides poleward of both
/// endpoints — the case a from-the-vertices bbox under-covers. The box is
/// conservative where a min/max box cannot represent the truth: an arc
/// crossing the antimeridian (or touching a pole) widens to the full
/// longitude range.
#[must_use]
pub fn arc_bbox(a: GeoPoint, b: GeoPoint) -> BBox {
    // Padding absorbs the round-off of points *computed* on the arc (e.g.
    // interpolation), so the box stays a true bound — the RectBounder trick.
    // 1e-9° is ~0.1 mm on Earth.
    const PAD_DEG: f64 = 1e-9;

    let (min_lat, max_lat) = s2arc::latitude_extremes(sp(a), sp(b));
    let min_lat = (min_lat.degrees() - PAD_DEG).max(-90.0);
    let max_lat = (max_lat.degrees() + PAD_DEG).min(90.0);

    // Longitude is monotonic along a (non-meridian) minor arc, so the span is
    // a.lon → a.lon + Δ with Δ the wrapped difference.
    let delta = (b.lon - a.lon + 540.0).rem_euclid(360.0) - 180.0;
    let (lon_w, lon_e) = if delta >= 0.0 {
        (a.lon - PAD_DEG, a.lon + delta + PAD_DEG)
    } else {
        (a.lon + delta - PAD_DEG, a.lon + PAD_DEG)
    };
    if max_lat >= 90.0 - 1e-12 || min_lat <= -90.0 + 1e-12 || lon_w < -180.0 || lon_e > 180.0 {
        // Pole-touching or antimeridian-crossing: a min/max box can't say it
        // tighter.
        return BBox::new(min_lat, -180.0, max_lat, 180.0);
    }
    BBox::new(min_lat, lon_w, max_lat, lon_e)
}

/// Arc-segment distances in metres on this sphere — the metric layer over
/// [`crate::space::s2::arc`]. On `Sphere::EARTH` these are the spherical
/// (great-circle) forms; they are closed-form and exact on the sphere.
impl Sphere {
    /// Signed cross-track distance from `p` to the **great circle** through
    /// `a → b`: positive when `p` is left of the direction of travel,
    /// negative right. Not clamped to the segment — see
    /// [`distance_to_arc_m`](Self::distance_to_arc_m).
    #[must_use]
    pub fn cross_track_m(self, p: GeoPoint, a: GeoPoint, b: GeoPoint) -> f64 {
        self.length(s2arc::cross_track(sp(p), sp(a), sp(b)))
            .as_meters()
    }

    /// Signed along-track distance: how far along the arc direction (from
    /// `a`, positive toward `b`) the perpendicular foot of `p` falls.
    #[must_use]
    pub fn along_track_m(self, p: GeoPoint, a: GeoPoint, b: GeoPoint) -> f64 {
        self.length(s2arc::along_track(sp(p), sp(a), sp(b)))
            .as_meters()
    }

    /// Distance from `p` to the arc **segment** `a → b` (perpendicular foot,
    /// clamped to the nearer endpoint) — the spherical point-to-segment
    /// distance.
    #[must_use]
    pub fn distance_to_arc_m(self, p: GeoPoint, a: GeoPoint, b: GeoPoint) -> f64 {
        self.length(s2arc::distance_to_arc(sp(p), sp(a), sp(b)))
            .as_meters()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::earth::haversine_m;

    fn g(lat: f64, lon: f64) -> GeoPoint {
        GeoPoint::from_lat_lon(lat, lon)
    }

    #[test]
    fn equator_meridian_arcs_cross_at_origin() {
        let x = intersection(g(0.0, -10.0), g(0.0, 10.0), g(-10.0, 0.0), g(10.0, 0.0))
            .expect("arcs cross");
        assert!(x.lat.abs() < 1e-12 && x.lon.abs() < 1e-12, "{x:?}");
        // Parallel-ish arcs don't cross.
        assert!(intersection(g(0.0, -10.0), g(0.0, 10.0), g(5.0, -10.0), g(5.0, 10.0)).is_none());
    }

    #[test]
    fn ring_contains_is_winding_insensitive_and_pole_safe() {
        let ccw = [
            g(-10.0, -10.0),
            g(-10.0, 10.0),
            g(10.0, 10.0),
            g(10.0, -10.0),
        ];
        let cw: Vec<GeoPoint> = ccw.iter().rev().copied().collect();
        for ring in [&ccw[..], &cw[..]] {
            assert!(ring_contains(ring, g(0.0, 0.0)));
            assert!(!ring_contains(ring, g(20.0, 0.0)));
        }
        // A ring around the north pole: lon/lat point-in-polygon fails here,
        // the spherical test must not.
        let polar: Vec<GeoPoint> = (0..12).map(|i| g(80.0, -180.0 + 30.0 * i as f64)).collect();
        assert!(ring_contains(&polar, g(89.0, 45.0)));
        assert!(!ring_contains(&polar, g(70.0, 45.0)));
    }

    #[test]
    fn cross_track_distance_on_earth() {
        // 1° north of an equatorial arc ≈ one degree of meridian arc.
        let d = Sphere::EARTH.cross_track_m(g(1.0, 5.0), g(0.0, 0.0), g(0.0, 10.0));
        let one_degree = haversine_m(g(0.0, 5.0), g(1.0, 5.0));
        assert!((d - one_degree).abs() < 1e-6, "{d} vs {one_degree}");
        // Right of track → negative.
        assert!(Sphere::EARTH.cross_track_m(g(-1.0, 5.0), g(0.0, 0.0), g(0.0, 10.0)) < 0.0);
    }

    #[test]
    fn distance_to_arc_clamps_to_endpoints() {
        let (a, b) = (g(0.0, 0.0), g(0.0, 10.0));
        // Abeam: perpendicular distance.
        let abeam = Sphere::EARTH.distance_to_arc_m(g(2.0, 5.0), a, b);
        assert!((abeam - Sphere::EARTH.cross_track_m(g(2.0, 5.0), a, b)).abs() < 1e-6);
        // Past the end: distance to the endpoint, not the circle.
        let beyond = Sphere::EARTH.distance_to_arc_m(g(0.0, 15.0), a, b);
        assert!((beyond - haversine_m(g(0.0, 15.0), b)).abs() < 1e-6);
    }

    #[test]
    fn arc_bbox_includes_clairaut_bulge() {
        let (a, b) = (g(45.0, 0.0), g(45.0, 90.0));
        let bbox = arc_bbox(a, b);
        // Sample the arc; every point must be inside.
        for i in 0..=100 {
            let p = interpolate(a, b, f64::from(i) / 100.0);
            assert!(bbox.contains(p), "{p:?} outside {bbox:?}");
        }
        // And the vertices-only box would NOT contain the midpoint.
        let naive = BBox::from_coords(&[a, b]);
        assert!(!naive.contains(midpoint(a, b)));
    }

    #[test]
    fn arc_bbox_widens_across_the_antimeridian() {
        let bbox = arc_bbox(g(10.0, 175.0), g(12.0, -175.0));
        assert!(bbox.contains(g(11.0, 179.5)) && bbox.contains(g(11.0, -179.5)));
    }
}
