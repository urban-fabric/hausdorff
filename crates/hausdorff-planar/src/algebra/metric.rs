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
//! Both namesake kernels are machine-checked in `riemann`, over the *exact*
//! nested folds computed here:
//! - **`Hausdorff/SetDistance`** proves the finite Hausdorff distance is a
//!   genuine metric on nonempty point sets — symmetric (`hausdorffDist_comm`),
//!   nonnegative, the triangle inequality (`hausdorffDist_triangle`), and
//!   vanishing *iff the sets are equal* (`hausdorffDist_eq_zero_iff`) — with
//!   `directed_le_iff` giving the certifying-checker reading of the directed
//!   distance and `directed_exists_witness` proving the max-min is *attained*
//!   (the basis for [`directed_hausdorff_witness`]).
//! - **`Hausdorff/Frechet`** proves the Eiter–Mannila recurrence
//!   ([`FrechetWorkspace::distance`]) is the discrete Fréchet distance:
//!   `frechet_isLeast` — the DP value is the *least* cost over every monotone
//!   coupling (attained by `frechet_reach`, a lower bound by `frechet_min`).
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
///
/// Proof-backed by `riemann` `Hausdorff/SetDistance.directedHausdorff`: the
/// `sup'`/`inf'` there *are* this outer-`max`/inner-`min` fold pair, and
/// `directed_le_iff` gives its certifying meaning (`≤ r` iff every point of `a`
/// is within `r` of some point of `b`); `directed_triangle` /
/// `directed_eq_zero_iff` are the directed half of the metric laws.
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

/// A **checkable certificate** for a directed Hausdorff distance: the pair of
/// points that *attains* it.
///
/// [`directed_hausdorff`] reports only the scalar `maxₚ∈ₐ minₑ∈ᵦ d(p, q)`. This
/// records *where* that maximum is realized — the worst point of the source set
/// (`source_index`) and its nearest neighbour in the target set
/// (`target_index`), whose separation is exactly [`distance`](Self::distance).
/// A caller holding the original slices can re-evaluate
/// `metric.distance(a[source_index], b[target_index])` and confirm it equals the
/// reported directed distance — the witness is *verifiable*, not just reported.
///
/// Proof-backed by `riemann` `Hausdorff/SetDistance.directed_exists_witness`:
/// the max-min is attained at concrete `p ∈ a`, `q ∈ b`, so the value is a
/// realizable certificate (the worst point and its nearest neighbour), not
/// merely a bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HausdorffWitness {
    /// Index into the *source* slice of the point furthest from the target set.
    pub source_index: usize,
    /// Index into the *target* slice of the source point's nearest neighbour.
    pub target_index: usize,
    /// The realized distance between the two indexed points — equal to
    /// [`directed_hausdorff`] over the same arguments.
    pub distance: f64,
}

/// Like [`directed_hausdorff`], but returns the attaining [`HausdorffWitness`]
/// — the worst point of `a`, its nearest neighbour in `b`, and their distance —
/// rather than just the scalar. Returns `None` if either set is empty.
///
/// `witness.distance == directed_hausdorff(metric, a, b).unwrap()`, and the pair
/// `(a[witness.source_index], b[witness.target_index])` realizes it. Ties resolve
/// to the lowest index (the first point of `a` attaining the max-min, and its
/// first nearest neighbour in `b`).
///
/// This is the certificate side of `riemann`
/// `Hausdorff/SetDistance.directed_exists_witness` — useful for *explaining* a
/// shape mismatch (which vertex is the outlier, and what it is closest to), not
/// only measuring it.
#[must_use]
pub fn directed_hausdorff_witness<M, P>(metric: &M, a: &[P], b: &[P]) -> Option<HausdorffWitness>
where
    M: MetricSpace<P>,
    P: Copy,
{
    if a.is_empty() || b.is_empty() {
        return None;
    }
    let mut best = HausdorffWitness {
        source_index: 0,
        target_index: 0,
        distance: f64::NEG_INFINITY,
    };
    for (i, &p) in a.iter().enumerate() {
        // The nearest neighbour of `p` in `b` (lowest index on ties).
        let (nearest_j, nearest_d) =
            b.iter()
                .enumerate()
                .fold((0usize, f64::INFINITY), |(bj, bd), (j, &q)| {
                    let d = metric.distance(p, q);
                    if d < bd { (j, d) } else { (bj, bd) }
                });
        // Keep the point of `a` whose nearest neighbour is farthest (the max-min).
        if nearest_d > best.distance {
            best = HausdorffWitness {
                source_index: i,
                target_index: nearest_j,
                distance: nearest_d,
            };
        }
    }
    Some(best)
}

/// (Undirected) Hausdorff distance between two point sets under `metric`: the
/// maximum of the two directed distances. Returns `None` if either set is empty.
///
/// Proof-backed by `riemann` `Hausdorff/SetDistance.hausdorffDist` — the
/// `max(directed a b, directed b a)` here — which is proved a genuine metric on
/// finite nonempty point sets: `hausdorffDist_comm` (symmetric),
/// `hausdorffDist_nonneg`, `hausdorffDist_triangle`, and `hausdorffDist_eq_zero_iff`
/// (vanishes *iff* the sets are equal).
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
///
/// Proof-backed by `riemann` `Hausdorff/Frechet.frechet_isLeast`: the DP value
/// this evaluates is the *least* leash over every monotone coupling of the two
/// polylines — i.e. it really is the discrete Fréchet distance, not just an
/// upper bound. See [`FrechetWorkspace::distance`] for the recurrence.
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
    ///
    /// The rolling row below is the `O(n)`-space evaluation of `riemann`
    /// `Hausdorff/Frechet.frechetDP` — the `max(min(left, up, diag), d)`
    /// recurrence (`frechetDP_succ_succ`), whose value `frechet_isLeast` proves is
    /// the least coupling cost.
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
                // `frechetDP_succ_succ`: min over the three predecessors (up, left,
                // diag), then the new pair's leash — the cheapest way to extend a
                // coupling to (i, j).
                let best_prev = row[j].min(row[j - 1]).min(diag);
                row[j] = best_prev.max(metric.distance(pi, q[j]));
                diag = up;
            }
        }
        Some(row[n - 1])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::r2::{Coord, Euclidean};

    fn pts(coords: &[(f64, f64)]) -> Vec<Coord> {
        coords.iter().map(|&(x, y)| Coord::new(x, y)).collect()
    }

    #[test]
    fn witness_attains_the_directed_distance() {
        let a = pts(&[(0.0, 0.0), (10.0, 0.0)]);
        let b = pts(&[(0.0, 1.0), (10.0, 2.0)]);

        let scalar = directed_hausdorff(&Euclidean, &a, &b).expect("non-empty");
        let w = directed_hausdorff_witness(&Euclidean, &a, &b).expect("non-empty");

        // The worst point of `a` is (10,0); its nearest neighbour in `b` is (10,2).
        assert_eq!(w.source_index, 1);
        assert_eq!(w.target_index, 1);
        assert_eq!(w.distance, 2.0);
        // The witness distance is the directed Hausdorff distance...
        assert_eq!(w.distance, scalar);
        // ...and the certificate is verifiable: re-evaluating the indexed pair
        // reproduces it (the `directed_exists_witness` property).
        assert_eq!(
            MetricSpace::distance(&Euclidean, a[w.source_index], b[w.target_index]),
            w.distance
        );
    }

    #[test]
    fn witness_is_asymmetric_like_the_scalar() {
        let a = pts(&[(0.0, 0.0), (10.0, 0.0)]);
        let b = pts(&[(0.0, 1.0), (10.0, 2.0)]);

        let ab = directed_hausdorff_witness(&Euclidean, &a, &b).expect("non-empty");
        let ba = directed_hausdorff_witness(&Euclidean, &b, &a).expect("non-empty");
        // b→a worst point is (10,2), nearest in a is (10,0), distance 2.
        assert_eq!(ab.distance, 2.0);
        assert_eq!(ba.distance, 2.0);
        // Indices are into the respective source/target slices.
        assert_eq!(ba.source_index, 1);
        assert_eq!(ba.target_index, 1);
    }

    #[test]
    fn witness_ties_resolve_to_lowest_index() {
        // Two points of `a` are equidistant (1.0) from their nearest neighbours;
        // the lowest source index wins, and its lowest nearest-neighbour index.
        let a = pts(&[(0.0, 1.0), (5.0, 1.0)]);
        let b = pts(&[(0.0, 0.0), (5.0, 0.0)]);
        let w = directed_hausdorff_witness(&Euclidean, &a, &b).expect("non-empty");
        assert_eq!(w.source_index, 0);
        assert_eq!(w.target_index, 0);
        assert_eq!(w.distance, 1.0);
    }

    #[test]
    fn witness_of_a_subset_is_zero() {
        // a ⊆ b ⟹ directed distance 0 (every point of `a` sits on a point of `b`).
        let a = pts(&[(1.0, 1.0)]);
        let b = pts(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)]);
        let w = directed_hausdorff_witness(&Euclidean, &a, &b).expect("non-empty");
        assert_eq!(w.distance, 0.0);
        assert_eq!(w.target_index, 1);
    }

    #[test]
    fn witness_none_on_empty() {
        let a = pts(&[(0.0, 0.0)]);
        let empty: Vec<Coord> = Vec::new();
        assert_eq!(directed_hausdorff_witness(&Euclidean, &a, &empty), None);
        assert_eq!(directed_hausdorff_witness(&Euclidean, &empty, &a), None);
    }
}
