//! Webs — spatially-embedded graphs, generic over the embedding space.
//!
//! A *web* is a graph whose nodes carry coordinates in a space `S ∈ {R2, S2}`: an
//! **R2 web** (planar — local or projected) or an **S2 web** (on the sphere —
//! Earth). It is the carrier type for spatial-network analysis. The *graph* half
//! of a web's character is ramsey's domain (degree, clustering, spectral, …); the
//! *spatial* half lives here, computed in the web's own metric — Euclidean for R2,
//! great-circle for S2.
//!
//! [`SpatialMetric`] abstracts the space (a point type plus a distance), so the
//! [`spatial_fingerprint`] is written once and instantiated at [`R2`] and [`S2`] —
//! the same generic-over-the-space pattern the `space::r2` / `space::s2` geometry
//! traits use. ramsey's graph fingerprint ⊕ this spatial fingerprint is the full
//! web fingerprint.

use crate::earth::{haversine_m, GeoPoint};

/// A space a web can be embedded in: a point type plus a distance between points.
pub trait SpatialMetric {
    /// The coordinate type of a node in this space.
    type Point: Copy;
    /// Distance between two embedded points, in the space's natural unit
    /// (Euclidean units for [`R2`], meters for [`S2`]).
    fn distance(a: Self::Point, b: Self::Point) -> f64;
}

/// The Euclidean plane ℝ². Points are `(x, y)`.
#[derive(Clone, Copy, Debug, Default)]
pub struct R2;

impl SpatialMetric for R2 {
    type Point = (f64, f64);
    #[inline]
    fn distance(a: (f64, f64), b: (f64, f64)) -> f64 {
        let dx = a.0 - b.0;
        let dy = a.1 - b.1;
        (dx * dx + dy * dy).sqrt()
    }
}

/// The sphere S² (Earth). Points are [`GeoPoint`]s; distance is the great-circle
/// (haversine) kernel from `crate::earth` — the one riemann proves correct.
#[derive(Clone, Copy, Debug, Default)]
pub struct S2;

impl SpatialMetric for S2 {
    type Point = GeoPoint;
    #[inline]
    fn distance(a: GeoPoint, b: GeoPoint) -> f64 {
        haversine_m(a, b)
    }
}

/// A spatially-embedded graph over space `S`: node coordinates plus undirected
/// edges as `(node, node)` index pairs.
pub struct Web<S: SpatialMetric> {
    /// Node coordinates in space `S`.
    pub nodes: Vec<S::Point>,
    /// Undirected edges as index pairs into [`Web::nodes`].
    pub edges: Vec<(usize, usize)>,
}

impl<S: SpatialMetric> Web<S> {
    /// Build a web from node coordinates and undirected edge index pairs.
    #[must_use]
    pub fn new(nodes: Vec<S::Point>, edges: Vec<(usize, usize)>) -> Self {
        Self { nodes, edges }
    }

    /// Number of nodes.
    #[must_use]
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges.
    #[must_use]
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }
}

/// The spatial half of a web's fingerprint — signals ramsey's graph fingerprint
/// cannot see, computed in the web's own metric.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialFingerprint {
    /// Node count.
    pub num_nodes: usize,
    /// Edge count.
    pub num_edges: usize,
    /// Mean nearest-neighbour distance.
    pub nn_mean: f64,
    /// Coefficient of variation of nearest-neighbour distances — the dimensionless
    /// dispersion of the node point process: ≈0 for a regular lattice, ≈0.5 for a
    /// Poisson field, larger for clustered points. (Space-agnostic; the
    /// area-normalised Clark–Evans index is the space-specific refinement.)
    pub nn_cv: f64,
    /// Mean edge length.
    pub edge_mean: f64,
    /// Coefficient of variation of edge lengths.
    pub edge_cv: f64,
}

/// Mean and coefficient of variation (`std / mean`) of a sample; `(0, 0)` for an
/// empty sample or zero mean.
fn mean_cv(xs: &[f64]) -> (f64, f64) {
    if xs.is_empty() {
        return (0.0, 0.0);
    }
    let n = xs.len() as f64;
    let mean = xs.iter().sum::<f64>() / n;
    if mean == 0.0 {
        return (0.0, 0.0);
    }
    let var = xs.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / n;
    (mean, var.sqrt() / mean)
}

/// Compute the spatial fingerprint of a web: nearest-neighbour dispersion (the
/// node point process) and the edge-length distribution, in the web's metric.
/// Generic over the space `S`, so the same code fingerprints R2 and S2 webs.
///
/// Nearest neighbours are found by brute force (`O(n²)`); for large webs swap in
/// the `hausdorff-index` grid.
#[must_use]
pub fn spatial_fingerprint<S: SpatialMetric>(web: &Web<S>) -> SpatialFingerprint {
    let n = web.nodes.len();
    let mut nn = Vec::with_capacity(n);
    for i in 0..n {
        let mut best = f64::INFINITY;
        for j in 0..n {
            if i == j {
                continue;
            }
            let d = S::distance(web.nodes[i], web.nodes[j]);
            if d < best {
                best = d;
            }
        }
        if best.is_finite() {
            nn.push(best);
        }
    }
    let edge_len: Vec<f64> = web
        .edges
        .iter()
        .map(|&(u, v)| S::distance(web.nodes[u], web.nodes[v]))
        .collect();
    let (nn_mean, nn_cv) = mean_cv(&nn);
    let (edge_mean, edge_cv) = mean_cv(&edge_len);
    SpatialFingerprint {
        num_nodes: n,
        num_edges: web.edges.len(),
        nn_mean,
        nn_cv,
        edge_mean,
        edge_cv,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid_web(k: usize) -> Web<R2> {
        let mut nodes = Vec::new();
        for y in 0..k {
            for x in 0..k {
                nodes.push((x as f64, y as f64));
            }
        }
        let idx = |x: usize, y: usize| y * k + x;
        let mut edges = Vec::new();
        for y in 0..k {
            for x in 0..k {
                if x + 1 < k {
                    edges.push((idx(x, y), idx(x + 1, y)));
                }
                if y + 1 < k {
                    edges.push((idx(x, y), idx(x, y + 1)));
                }
            }
        }
        Web::new(nodes, edges)
    }

    fn lcg(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((*state >> 33) as f64) / ((1u64 << 31) as f64)
    }

    fn random_web(n: usize, span: f64, seed: u64) -> Web<R2> {
        let mut s = seed;
        let mut nodes = Vec::new();
        for _ in 0..n {
            let x = lcg(&mut s) * span;
            let y = lcg(&mut s) * span;
            nodes.push((x, y));
        }
        Web::new(nodes, Vec::new())
    }

    #[test]
    fn grid_is_a_regular_point_process() {
        let fp = spatial_fingerprint(&grid_web(20));
        // Every nearest-neighbour distance is exactly 1, every edge length is 1.
        assert!((fp.nn_mean - 1.0).abs() < 1e-9);
        assert!(fp.nn_cv < 1e-9, "grid nn_cv should be ~0, got {}", fp.nn_cv);
        assert!(fp.edge_cv < 1e-9, "grid edge_cv should be ~0, got {}", fp.edge_cv);
    }

    #[test]
    fn random_is_more_dispersed_than_grid() {
        let grid = spatial_fingerprint(&grid_web(20)); // 400 nodes, regular
        let rand = spatial_fingerprint(&random_web(400, 19.0, 42)); // Poisson-ish
        assert!(grid.nn_cv < 0.05, "grid nn_cv {}", grid.nn_cv);
        assert!(rand.nn_cv > 0.2, "random nn_cv {} should be substantial", rand.nn_cv);
        assert!(rand.nn_cv > grid.nn_cv);
    }

    #[test]
    fn s2_web_fingerprints_with_great_circle_metric() {
        let nodes = vec![
            GeoPoint::from_lat_lon(40.00, -74.00),
            GeoPoint::from_lat_lon(40.01, -74.00),
            GeoPoint::from_lat_lon(40.00, -74.01),
            GeoPoint::from_lat_lon(40.02, -74.02),
        ];
        let edges = vec![(0, 1), (0, 2), (1, 3)];
        let fp = spatial_fingerprint(&Web::<S2>::new(nodes, edges));
        assert_eq!(fp.num_nodes, 4);
        // Distances are in meters (great-circle), so they are finite and positive.
        assert!(fp.nn_mean > 0.0 && fp.nn_mean.is_finite());
        assert!(fp.edge_mean > 0.0 && fp.edge_mean.is_finite());
    }
}
