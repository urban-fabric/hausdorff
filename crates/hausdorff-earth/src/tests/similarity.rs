//! Great-circle similarity measures: directed/symmetric Hausdorff and Frechet
//! distance between point sets and paths, measured on the curved Earth rather
//! than in a projected plane.
//!
//! The properties are the ones that make a metric usable: symmetry, zero on
//! equal sets, a known great-circle gap on hand-checkable inputs, `None` on
//! empty input, and agreement between the unit-typed wrappers and the raw
//! `_m` kernels.

use crate::earth::*;

fn pt(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

#[test]
fn hausdorff_is_symmetric_and_zero_on_equal_sets() {
    let a = [pt(0.0, 0.0), pt(0.0, 1.0), pt(1.0, 1.0)];
    assert_eq!(hausdorff_distance_m(&a, &a), Some(0.0));

    let b = [pt(0.0, 0.0), pt(0.0, 1.0)];
    let ab = hausdorff_distance_m(&a, &b).expect("non-empty");
    let ba = hausdorff_distance_m(&b, &a).expect("non-empty");
    assert!((ab - ba).abs() < 1e-9);
}

#[test]
fn hausdorff_equals_known_great_circle_gap() {
    let a = [pt(0.0, 0.0), pt(0.0, 1.0)];
    let b = [pt(0.0, 0.0)];
    let expected = haversine_m(pt(0.0, 1.0), pt(0.0, 0.0));
    assert!((hausdorff_distance_m(&a, &b).expect("non-empty") - expected).abs() < 1e-6);
}

#[test]
fn empty_sets_are_none() {
    let a = [pt(1.0, 2.0)];
    assert_eq!(hausdorff_distance_m(&a, &[]), None);
    assert_eq!(frechet_distance_m(&[], &a), None);
}

#[test]
fn frechet_zero_on_identical_path() {
    let path = [pt(0.0, 0.0), pt(0.0, 1.0), pt(1.0, 1.0)];
    assert!(frechet_distance_m(&path, &path).expect("non-empty") < 1e-6);
}

#[test]
fn frechet_of_parallel_offset_lines_is_the_offset() {
    let lower = [pt(0.0, 0.0), pt(0.0, 1.0), pt(0.0, 2.0)];
    let upper = [pt(1.0, 0.0), pt(1.0, 1.0), pt(1.0, 2.0)];
    let gap = haversine_m(pt(0.0, 0.0), pt(1.0, 0.0));
    assert!((frechet_distance_m(&lower, &upper).expect("non-empty") - gap).abs() < 1.0);
}

#[test]
fn frechet_dominates_endpoint_distance() {
    let p = [pt(0.0, 0.0), pt(2.0, 3.0)];
    let q = [pt(1.0, 1.0), pt(3.0, 5.0)];
    let d = frechet_distance_m(&p, &q).expect("non-empty");
    assert!(d >= haversine_m(p[0], q[0]) - 1e-6);
    assert!(d >= haversine_m(p[1], q[1]) - 1e-6);
}

#[test]
fn typed_wrappers_match_raw() {
    let a = [pt(0.0, 0.0), pt(0.0, 1.0)];
    let b = [pt(0.5, 0.5)];
    assert_eq!(
        hausdorff_distance(&a, &b).expect("non-empty").value(),
        hausdorff_distance_m(&a, &b).expect("non-empty")
    );
    assert_eq!(
        frechet_distance(&a, &b).expect("non-empty").value(),
        frechet_distance_m(&a, &b).expect("non-empty")
    );
}
