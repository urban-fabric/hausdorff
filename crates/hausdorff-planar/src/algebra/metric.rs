//! Metric spaces and the set/curve measures defined over *any* of them.
//!
//! After `spire`'s `algebra.MetricSpace`: a metric is anything that assigns a
//! non-negative distance to a pair of points. The crate's concrete metrics —
//! [`Euclidean`](crate::space::r2::Euclidean) and [`Haversine`](crate::space::r2::Haversine) —
//! already implement [`Distance`], so they become [`MetricSpace`]s for free via
//! the blanket impl below.
//!
//! The namesake measures of this crate — the [Hausdorff distance] and the
//! discrete [Fréchet distance] — are *defined over an arbitrary metric space*,
//! so they live here as generic algorithms rather than being welded to one
//! point type. The [`hausdorff_base::earth`] layer re-exposes the great-circle
//! specializations as thin wrappers.
//!
//! [Hausdorff distance]: https://en.wikipedia.org/wiki/Hausdorff_distance
//! [Fréchet distance]: https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance

use crate::space::r2::Distance;

/// A metric space over points of type `P`: a non-negative, symmetric distance.
///
/// Modelled on `spire`'s `MetricSpace[V, R]` (here the scalar `R` is fixed to
/// `f64`). Every [`Distance<P, P>`](crate::space::r2::Distance) is automatically a
/// `MetricSpace<P>`, so `Euclidean`/`Haversine` need no extra impl.
pub trait MetricSpace<P> {
    /// The distance between `a` and `b` in this metric's own unit (metres for
    /// the great-circle metrics, the plane's own unit for `Euclidean`).
    fn distance(&self, a: P, b: P) -> f64;
}

impl<M, P> MetricSpace<P> for M
where
    M: Distance<P, P>,
{
    #[inline]
    fn distance(&self, a: P, b: P) -> f64 {
        Distance::distance(self, a, b)
    }
}

/// Directed Hausdorff distance from `a` to `b` under `metric`: the greatest,
/// over every point of `a`, of the distance to the nearest point of `b`.
///
/// Asymmetric — `directed_hausdorff(m, a, b)` generally differs from
/// `directed_hausdorff(m, b, a)`. Returns `None` if either set is empty.
#[must_use]
pub fn directed_hausdorff<M, P>(metric: &M, a: &[P], b: &[P]) -> Option<f64>
where
    M: MetricSpace<P>,
    P: Copy,
{
    if a.is_empty() || b.is_empty() {
        return None;
    }
    Some(a.iter().fold(0.0_f64, |max_min, &p| {
        let nearest = b
            .iter()
            .fold(f64::INFINITY, |best, &q| best.min(metric.distance(p, q)));
        max_min.max(nearest)
    }))
}

/// (Undirected) Hausdorff distance between two point sets under `metric`: the
/// maximum of the two directed distances. Returns `None` if either set is empty.
#[must_use]
pub fn hausdorff_distance<M, P>(metric: &M, a: &[P], b: &[P]) -> Option<f64>
where
    M: MetricSpace<P>,
    P: Copy,
{
    directed_hausdorff(metric, a, b)
        .and_then(|ab| directed_hausdorff(metric, b, a).map(|ba| ab.max(ba)))
}

/// Discrete Fréchet distance between two polylines under `metric` — the
/// "dog walker on a leash" curve similarity, after Eiter & Mannila (1994).
///
/// Sensitive to ordering, unlike [`hausdorff_distance`]. Returns `None` if either
/// path is empty. `O(m·n)` time, `O(n)` space. Allocates one DP row per call;
/// for a batch of comparisons reuse a [`FrechetWorkspace`] to amortize it.
#[must_use]
pub fn frechet_distance<M, P>(metric: &M, p: &[P], q: &[P]) -> Option<f64>
where
    M: MetricSpace<P>,
    P: Copy,
{
    FrechetWorkspace::new().distance(metric, p, q)
}

/// Reusable scratch space for repeated discrete-Fréchet computations.
///
/// The dynamic program needs one `f64` row of width `|q|`. Allocating it per call
/// is wasteful when comparing many curve pairs (trajectory matching, nearest-curve
/// search), so a workspace amortizes the allocation across calls — the
/// faer/RustFFT "reusable scratch" pattern. The free [`frechet_distance`] is the
/// allocate-internally convenience form. Cheap to clone and `Send`/`Sync`, so one
/// per worker thread reuses its buffer across a batch.
#[derive(Debug, Default, Clone)]
pub struct FrechetWorkspace {
    row: Vec<f64>,
}

impl FrechetWorkspace {
    /// A new, empty workspace; its buffer grows to fit the first comparison and
    /// is reused (without reallocating) by subsequent same-or-smaller ones.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Discrete Fréchet distance between `p` and `q` under `metric`, reusing this
    /// workspace's row buffer. Identical result to [`frechet_distance`].
    /// Returns `None` if either path is empty.
    pub fn distance<M, P>(&mut self, metric: &M, p: &[P], q: &[P]) -> Option<f64>
    where
        M: MetricSpace<P>,
        P: Copy,
    {
        let (m, n) = (p.len(), q.len());
        if m == 0 || n == 0 {
            return None;
        }

        let row = &mut self.row;
        row.clear();
        row.resize(n, 0.0);
        row[0] = metric.distance(p[0], q[0]);
        for j in 1..n {
            row[j] = row[j - 1].max(metric.distance(p[0], q[j]));
        }
        for &pi in &p[1..] {
            let mut diag = row[0];
            row[0] = row[0].max(metric.distance(pi, q[0]));
            for j in 1..n {
                let up = row[j];
                let best_prev = row[j].min(row[j - 1]).min(diag);
                row[j] = best_prev.max(metric.distance(pi, q[j]));
                diag = up;
            }
        }
        Some(row[n - 1])
    }
}
