//! Set- and curve-similarity measures specialized to the Earth's great-circle
//! metric: the [Hausdorff distance] between two point sets and the discrete
//! [Fréchet distance] between two paths — both in **metres** on a sphere of
//! radius [`EARTH_RADIUS_M`](super::EARTH_RADIUS_M).
//!
//! These are the metric-generic algorithms in [`crate::algebra`] applied to the
//! [`Haversine`](crate::space::r2::Haversine) metric over [`GeoPoint`]; the generic
//! forms work over any [`MetricSpace`](crate::algebra::MetricSpace), so the same
//! code computes the planar/Euclidean variants too.
//!
//! Distance to an empty set/path is undefined, so these return `Option` and
//! yield `None` for an empty operand — consistent with the generic kernels.
//!
//! Ported from georust `geo` (`hausdorff_distance`, `frechet_distance`); see the
//! crate `NOTICE`.
//!
//! [Hausdorff distance]: https://en.wikipedia.org/wiki/Hausdorff_distance
//! [Fréchet distance]: https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance

use euclid::LengthMeters;

use crate::space::r2::Haversine;

use super::GeoPoint;

/// Directed Hausdorff distance from `a` to `b` in **metres**: the greatest, over
/// every point of `a`, of the great-circle distance to the nearest point of `b`.
///
/// This is asymmetric — `directed_hausdorff_m(a, b)` generally differs from
/// `directed_hausdorff_m(b, a)`. Returns `None` if either set is empty.
#[must_use]
pub fn directed_hausdorff_m(a: &[GeoPoint], b: &[GeoPoint]) -> Option<f64> {
    crate::algebra::directed_hausdorff(&Haversine, a, b)
}

/// (Undirected) Hausdorff distance between two point sets in **metres**: the
/// maximum of the two directed distances. The classic measure of how badly two
/// point sets mismatch — the namesake of this crate.
///
/// Returns `None` if either set is empty.
#[must_use]
pub fn hausdorff_distance_m(a: &[GeoPoint], b: &[GeoPoint]) -> Option<f64> {
    crate::algebra::hausdorff_distance(&Haversine, a, b)
}

/// Typed [`LengthMeters`] form of [`hausdorff_distance_m`] — prefer this so the
/// unit travels with the value.
#[must_use]
pub fn hausdorff_distance(a: &[GeoPoint], b: &[GeoPoint]) -> Option<LengthMeters> {
    hausdorff_distance_m(a, b).map(LengthMeters::new)
}

/// Discrete Fréchet distance between two polylines in **metres** — the classic
/// "dog walker on a leash" similarity of curves, after Eiter & Mannila's
/// *Computing Discrete Fréchet Distance* (1994).
///
/// Walks both paths monotonically and reports the shortest leash that lets both
/// reach the end; sensitive to ordering, unlike [`hausdorff_distance_m`].
/// Returns `None` if either path is empty.
///
/// Uses an `O(m·n)`-time, `O(min(m, n))`-space dynamic program over the
/// great-circle metric.
#[must_use]
pub fn frechet_distance_m(p: &[GeoPoint], q: &[GeoPoint]) -> Option<f64> {
    crate::algebra::frechet_distance(&Haversine, p, q)
}

/// Typed [`LengthMeters`] form of [`frechet_distance_m`].
#[must_use]
pub fn frechet_distance(p: &[GeoPoint], q: &[GeoPoint]) -> Option<LengthMeters> {
    frechet_distance_m(p, q).map(LengthMeters::new)
}
