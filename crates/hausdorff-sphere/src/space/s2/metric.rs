//! The **central-angle metric** on the unit sphere — the spherical instance of
//! the workspace's shared [`MetricSpace`](hausdorff_planar::algebra::MetricSpace)
//! abstraction, mirroring how the planar `Euclidean` / `Haversine` metrics plug
//! into the same generic set/curve measures.
//!
//! [`SphereMetric`] is the radius-free great-circle metric: the distance between
//! two [`S2Point`]s is the central angle between their direction vectors, in
//! **radians**. Because it implements
//! [`Distance<S2Point, S2Point>`](hausdorff_planar::space::r2::Distance) it
//! becomes a `MetricSpace<S2Point>` for free via the blanket impl, so the
//! generic [`directed_hausdorff`](hausdorff_planar::algebra::directed_hausdorff)
//! / [`hausdorff_distance`](hausdorff_planar::algebra::hausdorff_distance) /
//! [`frechet_distance`](hausdorff_planar::algebra::frechet_distance) algorithms
//! run on spherical point sets with no special case.
//!
//! Proof-backed by `riemann` `Hausdorff/SphereMetric`, which proves the
//! great-circle distance is a *genuine* metric: `sphereDist` is Mathlib's
//! `angle`, with `sphereDist_comm` (symmetry), `sphereDist_triangle` (the
//! triangle inequality — the geometric heart), `sphereDist_nonneg`,
//! `sphereDist_self` (self-distance zero), `sphereDist_eq_zero_iff` (identity of
//! indiscernibles — `dist = 0` *iff* the points coincide, the separation axiom,
//! with `sphereDist_pos_of_ne` its positive-definite form), `sphereDist_le_pi`
//! (with `sphereDist_eq_pi_iff` pinning the maximum: `dist = π` *iff* the points
//! are antipodal — so both extremes are characterized, `0`↔equal and `π`↔antipodal),
//! and `cos_sphereDist_unit` (`cos dist = p·q`).

use hausdorff_planar::space::r2::Distance;

use super::S2Point;

/// The unit-sphere **great-circle metric**: distance is the central angle (in
/// **radians**, `[0, π]`) between two directions — radius-free, like the rest of
/// `s2`. Multiply by the sphere radius for arc length in metres.
///
/// The spherical sibling of [`Euclidean`](hausdorff_planar::space::r2::Euclidean)
/// / [`Haversine`](hausdorff_planar::space::r2::Haversine); a zero-sized marker
/// that makes `S2Point` a [`MetricSpace`](hausdorff_planar::algebra::MetricSpace).
#[derive(Clone, Copy, Debug, Default)]
pub struct SphereMetric;

impl Distance<S2Point, S2Point> for SphereMetric {
    /// The central-angle (great-circle) distance in radians — [`S2Point::distance`]
    /// (riemann `sphereDist`, `Hausdorff/SphereMetric`).
    #[inline]
    fn distance(&self, a: S2Point, b: S2Point) -> f64 {
        a.distance(b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::earth::GeoPoint;
    use std::f64::consts::PI;

    fn pt(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    // The metric distance, called through the `Distance` trait `SphereMetric`
    // implements (the same one the `MetricSpace` blanket impl wraps).
    fn d(a: S2Point, b: S2Point) -> f64 {
        Distance::distance(&SphereMetric, a, b)
    }

    #[test]
    fn sphere_metric_axioms() {
        let p = pt(0.0, 0.0);
        let q = pt(0.0, 90.0);
        let r = pt(0.0, 45.0);
        // Self-distance is zero (sphereDist_self).
        assert!(d(p, p).abs() < 1e-12);
        // Symmetry (sphereDist_comm).
        assert!((d(p, q) - d(q, p)).abs() < 1e-12);
        // A concrete triangle inequality (sphereDist_triangle): the equator
        // detour via r is not shorter than the direct hop p → q.
        assert!(d(p, q) <= d(p, r) + d(r, q) + 1e-12);
        // The direct distance is π/2; non-negative throughout (sphereDist_nonneg).
        assert!((d(p, q) - PI / 2.0).abs() < 1e-12);
    }

    #[test]
    fn generic_hausdorff_runs_on_the_sphere() {
        use hausdorff_planar::algebra::hausdorff_distance;
        // The shared `MetricSpace`-generic Hausdorff distance instantiates at the
        // sphere through `SphereMetric` (via its `Distance` impl) with no special
        // case.
        let m = SphereMetric;
        let a = [pt(0.0, 0.0), pt(0.0, 10.0)];
        let b = [pt(0.0, 0.0), pt(0.0, 10.0)];
        // Identical sets → zero Hausdorff distance.
        assert!(hausdorff_distance(&m, &a, &b).unwrap().abs() < 1e-12);
        // A 1° latitude shift → ~1° (≈ 0.01745 rad) directed bound.
        let c = [pt(1.0, 0.0), pt(1.0, 10.0)];
        let h = hausdorff_distance(&m, &a, &c).unwrap();
        assert!((h - 1.0_f64.to_radians()).abs() < 1e-3, "{h}");
    }
}
