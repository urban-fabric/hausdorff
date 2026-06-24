use super::*;
use crate::space::r2::{Coord, Euclidean};

#[derive(Clone)]
struct CloneOnlyPoint(f64);

struct CloneOnlyLineMetric;

impl MetricSpace<CloneOnlyPoint> for CloneOnlyLineMetric {
    type Scalar = f64;

    fn distance(&self, a: CloneOnlyPoint, b: CloneOnlyPoint) -> Self::Scalar {
        (a.0 - b.0).abs()
    }
}

struct F32LineMetric;

impl MetricSpace<CloneOnlyPoint> for F32LineMetric {
    type Scalar = f32;

    fn distance(&self, a: CloneOnlyPoint, b: CloneOnlyPoint) -> Self::Scalar {
        (a.0 - b.0).abs() as f32
    }
}

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
fn set_and_curve_distances_accept_clone_only_points() {
    let a = vec![CloneOnlyPoint(0.0), CloneOnlyPoint(3.0)];
    let b = vec![CloneOnlyPoint(1.0), CloneOnlyPoint(4.0)];

    assert_eq!(directed_hausdorff(&CloneOnlyLineMetric, &a, &b), Some(1.0));
    let witness = directed_hausdorff_witness(&CloneOnlyLineMetric, &a, &b)
        .expect("non-empty clone-only sets");
    assert_eq!(witness.source_index, 0);
    assert_eq!(witness.target_index, 0);
    assert_eq!(witness.distance, 1.0);
    assert_eq!(hausdorff_distance(&CloneOnlyLineMetric, &a, &b), Some(1.0));
    assert_eq!(frechet_distance(&CloneOnlyLineMetric, &a, &b), Some(1.0));
}

#[test]
fn set_and_curve_distances_return_metric_scalar() {
    let a = vec![CloneOnlyPoint(0.0), CloneOnlyPoint(3.0)];
    let b = vec![CloneOnlyPoint(1.0), CloneOnlyPoint(4.0)];

    let directed: Option<f32> = directed_hausdorff(&F32LineMetric, &a, &b);
    let symmetric: Option<f32> = hausdorff_distance(&F32LineMetric, &a, &b);
    let frechet: Option<f32> = frechet_distance(&F32LineMetric, &a, &b);
    let witness: HausdorffWitness<f32> =
        directed_hausdorff_witness(&F32LineMetric, &a, &b).expect("non-empty");

    assert_eq!(directed, Some(1.0));
    assert_eq!(symmetric, Some(1.0));
    assert_eq!(frechet, Some(1.0));
    assert_eq!(witness.distance, 1.0);
}

#[test]
fn witness_is_asymmetric_like_the_scalar() {
    let a = pts(&[(0.0, 0.0), (10.0, 0.0)]);
    let b = pts(&[(0.0, 1.0), (10.0, 2.0)]);

    let ab = directed_hausdorff_witness(&Euclidean, &a, &b).expect("non-empty");
    let ba = directed_hausdorff_witness(&Euclidean, &b, &a).expect("non-empty");
    // b->a worst point is (10,2), nearest in a is (10,0), distance 2.
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
    // a subset of b means directed distance 0: every point of `a` sits on a point of `b`.
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
