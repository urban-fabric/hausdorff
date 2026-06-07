//! Tests for the cross-cutting `algebra` abstractions: the metric-generic
//! set/curve distances (here over the planar `Euclidean` metric, complementing
//! the `Haversine` coverage in the earth tests) and the reusable Fréchet
//! workspace.

use hausdorff::algebra::{
    FrechetWorkspace, directed_hausdorff, frechet_distance, hausdorff_distance,
};
use hausdorff::{Coord, Euclidean};

fn c(x: f64, y: f64) -> Coord {
    Coord::new(x, y)
}

#[test]
fn euclidean_hausdorff_known_values() {
    let a = [c(0.0, 0.0), c(1.0, 0.0), c(2.0, 0.0)];
    let b = [c(0.0, 3.0), c(1.0, 3.0), c(2.0, 3.0)];
    assert!((directed_hausdorff(&Euclidean, &a, &b) - 3.0).abs() < 1e-12);
    assert!((hausdorff_distance(&Euclidean, &a, &b) - 3.0).abs() < 1e-12);

    assert_eq!(hausdorff_distance(&Euclidean, &a, &[]), 0.0);

    let a2 = [c(0.0, 0.0), c(10.0, 0.0)];
    let b2 = [c(0.0, 0.0)];
    assert!((directed_hausdorff(&Euclidean, &a2, &b2) - 10.0).abs() < 1e-12);
    assert_eq!(directed_hausdorff(&Euclidean, &b2, &a2), 0.0);
}

#[test]
fn euclidean_frechet_known_values() {
    let p = [c(0.0, 0.0), c(1.0, 0.0), c(2.0, 0.0)];
    let q = [c(0.0, 1.0), c(1.0, 1.0), c(2.0, 1.0)];
    assert_eq!(frechet_distance(&Euclidean, &p, &p), 0.0);
    assert!((frechet_distance(&Euclidean, &p, &q) - 1.0).abs() < 1e-12);
}

#[test]
fn frechet_workspace_matches_free_fn_across_reuse() {
    let curves: [(Vec<Coord>, Vec<Coord>); 3] = [
        (
            vec![c(0.0, 0.0), c(1.0, 1.0), c(2.0, 0.0)],
            vec![c(0.0, 0.5), c(2.0, 0.5)],
        ),
        (
            vec![c(0.0, 0.0), c(5.0, 0.0)],
            vec![c(0.0, 2.0), c(5.0, 3.0)],
        ),
        (vec![c(0.0, 0.0)], vec![c(3.0, 4.0)]),
    ];

    let mut ws = FrechetWorkspace::new();
    for (p, q) in &curves {
        assert_eq!(
            ws.distance(&Euclidean, p, q),
            frechet_distance(&Euclidean, p, q)
        );
    }

    let big: Vec<Coord> = (0..50).map(|i| c(f64::from(i), 0.0)).collect();
    let small = [c(0.0, 2.0), c(49.0, 2.0)];
    assert_eq!(
        ws.distance(&Euclidean, &small, &big),
        frechet_distance(&Euclidean, &small, &big)
    );
    assert_eq!(
        ws.distance(&Euclidean, &small, &small),
        frechet_distance(&Euclidean, &small, &small)
    );
}
