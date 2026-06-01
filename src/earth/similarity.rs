//! Set- and curve-similarity measures: the [Hausdorff distance] between two
//! point sets and the discrete [Fréchet distance] between two paths.
//!
//! Both use the crate's spherical great-circle metric ([`haversine_m`]), so
//! results are metres on a sphere of radius [`EARTH_RADIUS_M`] — consistent with
//! the rest of `geo` rather than the planar/Euclidean form upstream uses.
//!
//! Ported from georust `geo` (`hausdorff_distance`, `frechet_distance`); see the
//! crate `NOTICE`.
//!
//! [Hausdorff distance]: https://en.wikipedia.org/wiki/Hausdorff_distance
//! [Fréchet distance]: https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance

use euclid::LengthMeters;

use super::{GeoPoint, haversine_m};

/// Directed Hausdorff distance from `a` to `b` in **metres**: the greatest, over
/// every point of `a`, of the great-circle distance to the nearest point of `b`.
///
/// This is asymmetric — `directed_hausdorff_m(a, b)` generally differs from
/// `directed_hausdorff_m(b, a)`. Returns `0.0` if either set is empty.
#[must_use]
pub fn directed_hausdorff_m(a: &[GeoPoint], b: &[GeoPoint]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    a.iter().fold(0.0_f64, |max_min, &p| {
        let nearest = b
            .iter()
            .fold(f64::INFINITY, |best, &q| best.min(haversine_m(p, q)));
        max_min.max(nearest)
    })
}

/// (Undirected) Hausdorff distance between two point sets in **metres**: the
/// maximum of the two directed distances. The classic measure of how badly two
/// point sets mismatch — the namesake of this crate.
///
/// Returns `0.0` if either set is empty.
#[must_use]
pub fn hausdorff_distance_m(a: &[GeoPoint], b: &[GeoPoint]) -> f64 {
    directed_hausdorff_m(a, b).max(directed_hausdorff_m(b, a))
}

/// Typed [`LengthMeters`] form of [`hausdorff_distance_m`] — prefer this so the
/// unit travels with the value.
#[must_use]
pub fn hausdorff_distance(a: &[GeoPoint], b: &[GeoPoint]) -> LengthMeters {
    LengthMeters::new(hausdorff_distance_m(a, b))
}

/// Discrete Fréchet distance between two polylines in **metres** — the classic
/// "dog walker on a leash" similarity of curves, after Eiter & Mannila's
/// *Computing Discrete Fréchet Distance* (1994).
///
/// Walks both paths monotonically and reports the shortest leash that lets both
/// reach the end; sensitive to ordering, unlike [`hausdorff_distance_m`].
/// Returns `0.0` if either path is empty.
///
/// Uses an `O(m·n)`-time, `O(min(m, n))`-space dynamic program over the
/// great-circle metric.
#[must_use]
pub fn frechet_distance_m(p: &[GeoPoint], q: &[GeoPoint]) -> f64 {
    let (m, n) = (p.len(), q.len());
    if m == 0 || n == 0 {
        return 0.0;
    }

    // One DP row of length `n`; `row[j]` holds the coupling measure ca(i, j).
    // We update it in place left-to-right, carrying the diagonal ca(i-1, j-1).
    let mut row = vec![0.0_f64; n];
    row[0] = haversine_m(p[0], q[0]);
    for j in 1..n {
        row[j] = row[j - 1].max(haversine_m(p[0], q[j]));
    }
    for i in 1..m {
        let mut diag = row[0]; // ca(i-1, 0)
        row[0] = row[0].max(haversine_m(p[i], q[0]));
        for j in 1..n {
            let up = row[j]; // ca(i-1, j) before overwrite
            let best_prev = row[j].min(row[j - 1]).min(diag);
            row[j] = best_prev.max(haversine_m(p[i], q[j]));
            diag = up; // becomes ca(i-1, j) == ca(i-1, (j+1)-1) next step
        }
    }
    row[n - 1]
}

/// Typed [`LengthMeters`] form of [`frechet_distance_m`].
#[must_use]
pub fn frechet_distance(p: &[GeoPoint], q: &[GeoPoint]) -> LengthMeters {
    LengthMeters::new(frechet_distance_m(p, q))
}
