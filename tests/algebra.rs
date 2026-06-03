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
    // Two horizontal point rows offset vertically by 3 → Hausdorff distance 3.
    let a = [c(0.0, 0.0), c(1.0, 0.0), c(2.0, 0.0)];
    let b = [c(0.0, 3.0), c(1.0, 3.0), c(2.0, 3.0)];
    assert!((directed_hausdorff(&Euclidean, &a, &b).unwrap() - 3.0).abs() < 1e-12);
    assert!((hausdorff_distance(&Euclidean, &a, &b).unwrap() - 3.0).abs() < 1e-12);

    // empty operand → None
    assert_eq!(hausdorff_distance(&Euclidean, &a, &[]), None);
    assert_eq!(directed_hausdorff(&Euclidean, &a, &[]), None);
    assert_eq!(directed_hausdorff(&Euclidean, &[], &a), None);
    assert_eq!(directed_hausdorff::<_, Coord>(&Euclidean, &[], &[]), None);

    // asymmetry: an extra far point in `a` raises directed(a, b) but not (b, a)
    let a2 = [c(0.0, 0.0), c(10.0, 0.0)];
    let b2 = [c(0.0, 0.0)];
    assert!((directed_hausdorff(&Euclidean, &a2, &b2).unwrap() - 10.0).abs() < 1e-12);
    assert_eq!(directed_hausdorff(&Euclidean, &b2, &a2), Some(0.0));
}

#[test]
fn euclidean_frechet_known_values() {
    let p = [c(0.0, 0.0), c(1.0, 0.0), c(2.0, 0.0)];
    let q = [c(0.0, 1.0), c(1.0, 1.0), c(2.0, 1.0)];
    // identical paths → Some(0); parallel paths offset by 1 → Some(1)
    assert_eq!(frechet_distance(&Euclidean, &p, &p), Some(0.0));
    assert!((frechet_distance(&Euclidean, &p, &q).unwrap() - 1.0).abs() < 1e-12);

    // empty path → None
    assert_eq!(frechet_distance(&Euclidean, &p, &[]), None);
    assert_eq!(frechet_distance(&Euclidean, &[], &q), None);
    assert_eq!(frechet_distance::<_, Coord>(&Euclidean, &[], &[]), None);
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

    // one workspace reused across comparisons of varying width must match the
    // allocate-internally free function exactly, every time.
    let mut ws = FrechetWorkspace::new();
    for (p, q) in &curves {
        assert_eq!(
            ws.distance(&Euclidean, p, q),
            frechet_distance(&Euclidean, p, q)
        );
    }

    // sizing the buffer up (width 50) then back down must not corrupt results.
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
