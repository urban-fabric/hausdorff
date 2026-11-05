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

use super::coords::MAX_LEVEL;
use crate::space::s2::S2Point;

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

/// An S2 cell **length metric** (Google S2's `LengthMetric`): a one-dimensional
/// cell measure — minimum width, average edge length, … — that scales as
/// `deriv · 2^(−level)` and is returned in **radians**, like [`SphereMetric`]
/// (multiply by the sphere radius for metres). The deriv constants are those of
/// S2's **quadratic** projection, the one this module uses (`coords`).
///
/// The level helpers turn a target arc length into the S2 level whose cells meet
/// it. The motivating use: partition a point set so every pair within a search
/// radius `r` lies in one cell or its neighbours, by picking
/// [`MIN_WIDTH.get_max_level(r)`](LengthMetric::get_max_level) — because
/// [`MIN_WIDTH`]`.value(level)` lower-bounds *every* cell's width at that level,
/// `value(level) >= r` guarantees the 1-ring covers the radius.
#[derive(Clone, Copy, Debug)]
pub struct LengthMetric {
    deriv: f64,
}

impl LengthMetric {
    /// The metric value at `level`, in radians: `deriv · 2^(−level)`.
    #[must_use]
    pub fn value(self, level: u32) -> f64 {
        self.deriv * 2f64.powi(-(level as i32))
    }

    /// The smallest (coarsest) level whose value is **at most** `value` radians,
    /// clamped to `[0, MAX_LEVEL]`. Non-positive `value` clamps to the finest
    /// level. Exact at level boundaries (corrected for `log2` rounding).
    #[must_use]
    pub fn get_min_level(self, value: f64) -> u32 {
        if value <= 0.0 {
            return MAX_LEVEL;
        }
        let mut level = (self.deriv / value)
            .log2()
            .ceil()
            .clamp(0.0, f64::from(MAX_LEVEL)) as u32;
        while level > 0 && self.value(level - 1) <= value {
            level -= 1;
        }
        while level < MAX_LEVEL && self.value(level) > value {
            level += 1;
        }
        level
    }

    /// The largest (finest) level whose value is **at least** `value` radians,
    /// clamped to `[0, MAX_LEVEL]`. Non-positive `value` clamps to the finest
    /// level. Exact at level boundaries (corrected for `log2` rounding).
    #[must_use]
    pub fn get_max_level(self, value: f64) -> u32 {
        if value <= 0.0 {
            return MAX_LEVEL;
        }
        let mut level = (self.deriv / value)
            .log2()
            .floor()
            .clamp(0.0, f64::from(MAX_LEVEL)) as u32;
        while level < MAX_LEVEL && self.value(level + 1) >= value {
            level += 1;
        }
        while level > 0 && self.value(level) < value {
            level -= 1;
        }
        level
    }

    /// The level whose value is closest to `value` radians (S2's geometric-mean
    /// rule: `get_min_level(√2 · value)`).
    #[must_use]
    pub fn get_closest_level(self, value: f64) -> u32 {
        self.get_min_level(std::f64::consts::SQRT_2 * value)
    }
}

/// S2 **minimum cell width** (quadratic-projection deriv `2√2/3`): a guaranteed
/// lower bound on every cell's width at a given level. Pick a partition level
/// with [`get_max_level`](LengthMetric::get_max_level) so that a cell plus its
/// neighbours still covers a target radius.
pub const MIN_WIDTH: LengthMetric = LengthMetric {
    deriv: 2.0 * std::f64::consts::SQRT_2 / 3.0,
};

/// S2 **average cell edge length** (quadratic-projection deriv `1.459213447877158`).
pub const AVG_EDGE: LengthMetric = LengthMetric {
    deriv: 1.459_213_447_877_158,
};

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

    #[test]
    fn length_metric_halves_each_level_and_round_trips() {
        for level in 0..MAX_LEVEL {
            let ratio = MIN_WIDTH.value(level + 1) / MIN_WIDTH.value(level);
            assert!((ratio - 0.5).abs() < 1e-12, "level {level}");
        }
        // The value at a level maps back to exactly that level (both directions).
        for level in 0..=MAX_LEVEL {
            let v = MIN_WIDTH.value(level);
            assert_eq!(MIN_WIDTH.get_max_level(v), level, "max @ {level}");
            assert_eq!(MIN_WIDTH.get_min_level(v), level, "min @ {level}");
        }
    }

    #[test]
    fn get_max_level_brackets_the_radius_for_completeness() {
        // The chosen level's width is >= the radius and the next finer level's is
        // <: exactly the bracket that makes a cell-plus-1-ring cover the radius.
        for &r in &[1e-1, 1e-2, 1e-3, 1e-5_f64] {
            let level = MIN_WIDTH.get_max_level(r);
            assert!(MIN_WIDTH.value(level) >= r, "r={r} level={level}");
            if level < MAX_LEVEL {
                assert!(MIN_WIDTH.value(level + 1) < r, "r={r} level={level}");
            }
        }
        // Clamping: non-positive → finest; a radius larger than a face cell → coarsest.
        assert_eq!(MIN_WIDTH.get_max_level(0.0), MAX_LEVEL);
        assert_eq!(MIN_WIDTH.get_max_level(100.0), 0);
        // AVG_EDGE is coarser (longer) than MIN_WIDTH at every level.
        for level in 0..=MAX_LEVEL {
            assert!(
                AVG_EDGE.value(level) > MIN_WIDTH.value(level),
                "level {level}"
            );
        }
    }
}
