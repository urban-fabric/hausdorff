//! Set- and curve-distances between whole geometries.
//!
//! The [`crate::algebra::metric`] measures are defined over point slices in an
//! arbitrary [`MetricSpace`]; here they are lifted onto the [`Geometry`] enum by
//! flattening each geometry's coordinate sequence ([`CoordsIter`]). Pair them
//! with [`Euclidean`](super::Euclidean) for planar units or
//! [`Haversine`](super::Haversine) for great-circle metres.
//!
//! Distance to an empty geometry is undefined, so these return `Option<f64>` and
//! yield `None` for an empty operand — surfacing emptiness rather than letting
//! `0.0` read as "identical". The underlying [`crate::algebra`] kernels also
//! return `Option`, so the same semantics hold at every layer.

use crate::algebra::{
    HausdorffWitness, MetricSpace, directed_hausdorff,
    directed_hausdorff_witness as directed_hausdorff_witness_points,
    frechet_distance as frechet_points, hausdorff_distance as hausdorff_points,
};

use super::{Coord, CoordsIter, Geometry};

/// (Undirected) Hausdorff distance between the coordinate sets of two geometries
/// under `metric`: the greater of the two directed distances, i.e. the largest
/// distance from a vertex of either geometry to the nearest vertex of the other.
/// Symmetric and order-independent. Returns `None` when either geometry has no
/// coordinates.
#[must_use]
pub fn hausdorff_distance<M: MetricSpace<Coord, Scalar = f64>>(
    metric: &M,
    a: &Geometry,
    b: &Geometry,
) -> Option<f64> {
    let _span = crate::profile::span("hausdorff_distance");
    let pa: Vec<Coord> = a.coords_iter().collect();
    let pb: Vec<Coord> = b.coords_iter().collect();
    hausdorff_points(metric, &pa, &pb)
}

/// Directed Hausdorff distance from `a` to `b` under `metric`: the greatest, over
/// every vertex of `a`, of the distance to the nearest vertex of `b`. Asymmetric.
/// Returns `None` when either geometry has no coordinates.
#[must_use]
pub fn directed_hausdorff_distance<M: MetricSpace<Coord, Scalar = f64>>(
    metric: &M,
    a: &Geometry,
    b: &Geometry,
) -> Option<f64> {
    let _span = crate::profile::span("directed_hausdorff_distance");
    let pa: Vec<Coord> = a.coords_iter().collect();
    let pb: Vec<Coord> = b.coords_iter().collect();
    directed_hausdorff(metric, &pa, &pb)
}

/// Like [`directed_hausdorff_distance`], but returns the attaining
/// [`HausdorffWitness`] — the vertex of `a` farthest from `b`, its nearest
/// neighbour in `b`, and that distance — rather than just the scalar. The
/// indices refer to each geometry's flattened coordinate order
/// ([`CoordsIter`]). Returns `None` when either geometry has no coordinates.
///
/// Use this to *explain* a directed shape mismatch — which vertex is the
/// outlier and what it sits closest to — not only measure it. The certificate
/// is verifiable (re-evaluate the indexed pair); see `riemann`
/// `Hausdorff/SetDistance.directed_exists_witness`.
#[must_use]
pub fn directed_hausdorff_witness<M: MetricSpace<Coord, Scalar = f64>>(
    metric: &M,
    a: &Geometry,
    b: &Geometry,
) -> Option<HausdorffWitness> {
    let _span = crate::profile::span("directed_hausdorff_witness");
    let pa: Vec<Coord> = a.coords_iter().collect();
    let pb: Vec<Coord> = b.coords_iter().collect();
    directed_hausdorff_witness_points(metric, &pa, &pb)
}

/// Discrete Fréchet distance between the coordinate *sequences* of two geometries
/// under `metric` — the "dog walker on a leash" curve similarity. Order-sensitive
/// (the vertex order of a `LineString` is its path direction), unlike
/// [`hausdorff_distance`]. Returns `None` when either geometry has no coordinates.
#[must_use]
pub fn frechet_distance<M: MetricSpace<Coord, Scalar = f64>>(
    metric: &M,
    a: &Geometry,
    b: &Geometry,
) -> Option<f64> {
    let _span = crate::profile::span("frechet_distance");
    let pa: Vec<Coord> = a.coords_iter().collect();
    let pb: Vec<Coord> = b.coords_iter().collect();
    frechet_points(metric, &pa, &pb)
}

#[cfg(test)]
mod tests {
    use super::super::{Euclidean, Haversine, LineString, Point};
    use super::*;

    fn line(coords: &[(f64, f64)]) -> Geometry {
        Geometry::LineString(LineString::new(
            coords.iter().map(|&(x, y)| Coord::new(x, y)).collect(),
        ))
    }

    #[test]
    fn identical_geometries_have_zero_distance() {
        let a = line(&[(0.0, 0.0), (0.001, 0.0), (0.002, 0.0)]);
        assert_eq!(hausdorff_distance(&Haversine, &a, &a), Some(0.0));
        assert_eq!(frechet_distance(&Haversine, &a, &a), Some(0.0));
    }

    #[test]
    fn parallel_lines_offset_north_measure_the_offset() {
        let a = line(&[(0.0, 0.0), (0.001, 0.0), (0.002, 0.0)]);
        let b = line(&[(0.0, 0.001), (0.001, 0.001), (0.002, 0.001)]);
        let h = hausdorff_distance(&Haversine, &a, &b).expect("non-empty");
        let f = frechet_distance(&Haversine, &a, &b).expect("non-empty");
        assert!((110.0..112.5).contains(&h), "hausdorff {h}");
        assert!((110.0..112.5).contains(&f), "frechet {f}");
    }

    #[test]
    fn witness_identifies_the_outlier_vertex() {
        // a's second vertex (5,0) is the outlier: its nearest neighbour in b is
        // (5,1) at distance 1; a's first vertex sits exactly on b's first vertex.
        let a = line(&[(0.0, 0.0), (5.0, 0.0)]);
        let b = line(&[(0.0, 0.0), (5.0, 1.0)]);

        let w = directed_hausdorff_witness(&Euclidean, &a, &b).expect("non-empty");
        assert_eq!(w.source_index, 1);
        assert_eq!(w.target_index, 1);
        assert_eq!(w.distance, 1.0);
        // The witness distance is the directed Hausdorff distance.
        assert_eq!(
            Some(w.distance),
            directed_hausdorff_distance(&Euclidean, &a, &b)
        );
    }

    #[test]
    fn witness_none_on_empty_geometry() {
        let a = line(&[(0.0, 0.0), (1.0, 1.0)]);
        let empty = Geometry::MultiPoint(super::super::MultiPoint::new(Vec::<Point>::new()));
        assert_eq!(directed_hausdorff_witness(&Euclidean, &a, &empty), None);
        assert_eq!(directed_hausdorff_witness(&Euclidean, &empty, &a), None);
    }

    #[test]
    fn empty_geometry_yields_none() {
        let a = line(&[(0.0, 0.0), (1.0, 1.0)]);
        let empty = Geometry::MultiPoint(super::super::MultiPoint::new(Vec::<Point>::new()));
        assert_eq!(hausdorff_distance(&Haversine, &a, &empty), None);
        assert_eq!(directed_hausdorff_distance(&Haversine, &empty, &a), None);
        assert_eq!(frechet_distance(&Haversine, &empty, &empty), None);
    }
}
