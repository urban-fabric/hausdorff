//! Conformance for the native Martinez-Rueda BooleanOps against geo's BooleanOps
//! (areas of the results must match).

// NOTE: hausdorff::r2::bool_ops is a WORK IN PROGRESS — the Martinez-Rueda
// clipper is correct for disjoint/nested polygons but not yet for crossing /
// overlapping ones (those cases are #[ignore]d). compat keeps geo's BooleanOps
// as the blessed path until this is complete.
use geo::BooleanOps as _; // reference, for method resolution on geo::MultiPolygon
use hausdorff::r2::bool_ops::BooleanOps;
use hausdorff::{Area, LineString, MultiPolygon, Polygon};

fn hpoly(coords: &[(f64, f64)]) -> Polygon {
    Polygon::new(LineString::from(coords.to_vec()), vec![])
}
fn gpoly(coords: &[(f64, f64)]) -> geo::Polygon<f64> {
    geo::Polygon::new(
        geo::LineString(coords.iter().map(|&(x, y)| geo::Coord { x, y }).collect()),
        vec![],
    )
}

fn area(mp: &MultiPolygon) -> f64 {
    mp.unsigned_area()
}
fn garea(mp: &geo::MultiPolygon<f64>) -> f64 {
    use geo::Area;
    mp.unsigned_area()
}

fn check(a: &[(f64, f64)], b: &[(f64, f64)]) {
    let (ha, hb) = (hpoly(a), hpoly(b));
    let (ga, gb) = (gpoly(a), gpoly(b));
    let gma = geo::MultiPolygon(vec![ga.clone()]);
    let gmb = geo::MultiPolygon(vec![gb.clone()]);

    let eps = 1e-6;
    assert!(
        (area(&ha.intersection(&hb)) - garea(&gma.intersection(&gmb))).abs() < eps,
        "intersection area mismatch: {} vs {}",
        area(&ha.intersection(&hb)),
        garea(&gma.intersection(&gmb))
    );
    assert!(
        (area(&ha.union(&hb)) - garea(&gma.union(&gmb))).abs() < eps,
        "union area mismatch: {} vs {}",
        area(&ha.union(&hb)),
        garea(&gma.union(&gmb))
    );
    assert!(
        (area(&ha.difference(&hb)) - garea(&gma.difference(&gmb))).abs() < eps,
        "difference area mismatch: {} vs {}",
        area(&ha.difference(&hb)),
        garea(&gma.difference(&gmb))
    );
}

#[test]
#[ignore = "Martinez clipper WIP: crossing cases"]
fn overlapping_squares() {
    check(
        &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        &[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)],
    );
}

#[test]
fn disjoint_and_nested() {
    // disjoint
    check(
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        &[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0)],
    );
    // nested (b inside a)
    check(
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        &[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0)],
    );
}

#[test]
#[ignore = "Martinez clipper WIP: crossing cases"]
fn concave_overlap() {
    check(
        &[(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (3.0, 3.0), (0.0, 6.0)],
        &[(1.0, 1.0), (5.0, 1.0), (5.0, 5.0), (1.0, 5.0)],
    );
}
