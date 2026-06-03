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
    MetricSpace, directed_hausdorff, frechet_distance as frechet_points,
    hausdorff_distance as hausdorff_points,
};

use super::{Coord, CoordsIter, Geometry};

/// (Undirected) Hausdorff distance between the coordinate sets of two geometries
/// under `metric`: the greater of the two directed distances, i.e. the largest
/// distance from a vertex of either geometry to the nearest vertex of the other.
/// Symmetric and order-independent. Returns `None` when either geometry has no
/// coordinates.
#[must_use]
pub fn hausdorff_distance<M: MetricSpace<Coord>>(
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
pub fn directed_hausdorff_distance<M: MetricSpace<Coord>>(
    metric: &M,
    a: &Geometry,
    b: &Geometry,
) -> Option<f64> {
    let _span = crate::profile::span("directed_hausdorff_distance");
    let pa: Vec<Coord> = a.coords_iter().collect();
    let pb: Vec<Coord> = b.coords_iter().collect();
    directed_hausdorff(metric, &pa, &pb)
}

/// Discrete Fréchet distance between the coordinate *sequences* of two geometries
/// under `metric` — the "dog walker on a leash" curve similarity. Order-sensitive
/// (the vertex order of a `LineString` is its path direction), unlike
/// [`hausdorff_distance`]. Returns `None` when either geometry has no coordinates.
#[must_use]
pub fn frechet_distance<M: MetricSpace<Coord>>(
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
    use super::super::{Haversine, LineString, Point};
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
    fn empty_geometry_yields_none() {
        let a = line(&[(0.0, 0.0), (1.0, 1.0)]);
        let empty = Geometry::MultiPoint(super::super::MultiPoint::new(Vec::<Point>::new()));
        assert_eq!(hausdorff_distance(&Haversine, &a, &empty), None);
        assert_eq!(directed_hausdorff_distance(&Haversine, &empty, &a), None);
        assert_eq!(frechet_distance(&Haversine, &empty, &empty), None);
    }
}
