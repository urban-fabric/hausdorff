use super::*;
use crate::space::r2::{Area, Validation};

fn poly(rings: &[&[(f64, f64)]]) -> Polygon {
    let ring = |r: &[(f64, f64)]| LineString::from(r.to_vec());
    Polygon::new(ring(rings[0]), rings[1..].iter().map(|r| ring(r)).collect())
}

fn area(mp: &MultiPolygon) -> f64 {
    mp.unsigned_area()
}

#[test]
fn bowtie_splits_into_two_components() {
    // The boundary crosses itself at (5, 5); each lobe is a triangle of area 25.
    let bowtie = poly(&[&[(0., 0.), (0., 10.), (10., 0.), (10., 10.), (0., 0.)]]);
    assert!(!bowtie.is_valid(), "the fixture should start invalid");

    let repaired = bowtie.make_valid().expect("bowtie is repairable");
    assert!(repaired.is_valid());
    assert_eq!(repaired.0.len(), 2, "one component per lobe");
    assert!((area(&repaired) - 50.0).abs() < 1e-9, "{}", area(&repaired));
}

#[test]
fn dangling_spur_cancels_against_itself() {
    // Out to (15, 5) and straight back: covered twice, so it is not boundary.
    let spur = poly(&[&[
        (0., 0.),
        (10., 0.),
        (15., 5.),
        (10., 0.),
        (10., 10.),
        (0., 10.),
        (0., 0.),
    ]]);
    let repaired = spur.make_valid().expect("spur is repairable");
    assert!(repaired.is_valid());
    assert!((area(&repaired) - 100.0).abs() < 1e-9, "{}", area(&repaired));
}

#[test]
fn hole_is_subtracted() {
    let with_hole = poly(&[
        &[(0., 0.), (10., 0.), (10., 10.), (0., 10.), (0., 0.)],
        &[(2., 2.), (2., 4.), (4., 4.), (4., 2.), (2., 2.)],
    ]);
    let repaired = with_hole.make_valid().expect("repairable");
    assert!(repaired.is_valid());
    assert!((area(&repaired) - 96.0).abs() < 1e-9, "{}", area(&repaired));
}

#[test]
fn valid_input_keeps_its_area() {
    let square = poly(&[&[(0., 0.), (4., 0.), (4., 4.), (0., 4.), (0., 0.)]]);
    let repaired = square.make_valid().expect("repairable");
    assert!(repaired.is_valid());
    assert_eq!(repaired.0.len(), 1);
    assert!((area(&repaired) - 16.0).abs() < 1e-9, "{}", area(&repaired));
}

#[test]
fn overlapping_components_cancel_where_they_double_up() {
    // Two unit squares overlapping in a 1x2 strip. Under odd-even the overlap
    // is covered twice, so it drops out and the result is the symmetric
    // difference: 6 + 6 - 2 * 2 = 8.
    let two = MultiPolygon::new(vec![
        poly(&[&[(0., 0.), (3., 0.), (3., 2.), (0., 2.), (0., 0.)]]),
        poly(&[&[(2., 0.), (5., 0.), (5., 2.), (2., 2.), (2., 0.)]]),
    ]);
    let repaired = two.make_valid().expect("repairable");
    assert!(repaired.is_valid());
    assert!((area(&repaired) - 8.0).abs() < 1e-9, "{}", area(&repaired));
}

#[test]
fn non_finite_coordinates_are_rejected() {
    let bad = poly(&[&[(0., 0.), (f64::NAN, 0.), (1., 1.), (0., 0.)]]);
    assert_eq!(bad.make_valid(), Err(MakeValidError::NonFiniteCoordinate));
}

#[test]
fn figure_eight_touching_at_a_point_keeps_both_lobes() {
    // Two squares meeting at exactly one vertex — the pinch case. Both lobes
    // are area 1, and neither may be swallowed as a hole of the other.
    let pinched = poly(&[&[
        (0., 0.),
        (1., 0.),
        (1., 1.),
        (0., 1.),
        (0., 0.),
        (-1., 0.),
        (-1., -1.),
        (0., -1.),
        (0., 0.),
    ]]);
    let repaired = pinched.make_valid().expect("repairable");
    assert!(repaired.is_valid());
    assert!((area(&repaired) - 2.0).abs() < 1e-9, "{}", area(&repaired));
}
