use super::*;

const EPS: f64 = 1e-12;

fn assert_coords_eq(a: &LineString, b: &[(f64, f64)]) {
    assert_eq!(a.0.len(), b.len());
    for (c, &(x, y)) in a.0.iter().zip(b) {
        assert!(
            (c.x - x).abs() < EPS && (c.y - y).abs() < EPS,
            "({}, {}) vs ({x}, {y})",
            c.x,
            c.y
        );
    }
}

#[test]
fn identity_and_is_identity() {
    let t = AffineTransform::identity();
    assert!(t.is_identity());
    assert!(!t.translated(1.0, 2.0).is_identity());
    assert!(t.translated(1.0, 2.0).translated(-1.0, -2.0).is_identity());
    assert_eq!(t.apply(Coord::new(3.0, -4.0)), Coord::new(3.0, -4.0));
}

#[test]
fn compose_runs_in_order() {
    let scale = AffineTransform::scale(2.0, 2.0, Coord::new(0.0, 0.0));
    let translate = AffineTransform::translate(10.0, 0.0);
    let scale_then_translate = scale.compose(&translate);
    assert_eq!(
        scale_then_translate.apply(Coord::new(1.0, 1.0)),
        Coord::new(12.0, 2.0)
    );
    let translate_then_scale = translate.compose(&scale);
    assert_eq!(
        translate_then_scale.apply(Coord::new(1.0, 1.0)),
        Coord::new(22.0, 2.0)
    );
    let many = AffineTransform::identity().compose_many(&[scale, translate]);
    assert_eq!(many.apply(Coord::new(1.0, 1.0)), Coord::new(12.0, 2.0));
}

#[test]
fn inverse_round_trips() {
    let t = AffineTransform::translate(3.0, -7.0)
        .rotated(31.0, Coord::new(1.0, 2.0))
        .scaled(2.0, 0.5, Coord::new(-3.0, 4.0));
    let inv = t.inverse().unwrap();
    let round = t.compose(&inv);
    let c = Coord::new(5.5, -2.25);
    let rc = round.apply(c);
    assert!((rc.x - c.x).abs() < EPS && (rc.y - c.y).abs() < EPS);
    assert!(
        AffineTransform::scale(0.0, 1.0, Coord::new(0.0, 0.0))
            .inverse()
            .is_none()
    );
}

#[test]
fn translate_line_string() {
    let ls = LineString::from(vec![(0.0, 0.0), (5.0, 5.0), (10.0, 10.0)]);
    let translated = ls.translate(1.5, 3.5);
    assert_coords_eq(&translated, &[(1.5, 3.5), (6.5, 8.5), (11.5, 13.5)]);
    let mut in_place = ls;
    in_place.translate_mut(1.5, 3.5);
    assert_eq!(in_place, translated);
}

#[test]
fn rotate_point_about_itself_is_fixed() {
    let p = Point::new(1.0, 5.0);
    assert_eq!(p.rotate_around_centroid(30.0), p);
    let rotated = p.rotate_around_point(90.0, Point::new(0.0, 0.0));
    assert!((rotated.x() + 5.0).abs() < EPS && (rotated.y() - 1.0).abs() < EPS);
}

#[test]
fn scale_about_bounding_box_centre() {
    let ls = LineString::from(vec![(0.0, 0.0), (10.0, 10.0)]);
    assert_coords_eq(&ls.scale(2.0), &[(-5.0, -5.0), (15.0, 15.0)]);
    assert_coords_eq(&ls.scale_xy(2.0, 4.0), &[(-5.0, -15.0), (15.0, 25.0)]);
    assert_coords_eq(
        &ls.scale_around_point(2.0, 4.0, Coord::new(100.0, 100.0)),
        &[(-100.0, -300.0), (-80.0, -260.0)],
    );
}

#[test]
fn skew_zero_is_identity() {
    let poly = Polygon::new(
        LineString::from(vec![(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        vec![],
    );
    assert_eq!(poly.skew(0.0), poly);
    let skewed = poly.skew(30.0);
    assert!((skewed.exterior().0[0].x - -2.886751345948129).abs() < 1e-9);
}

#[test]
fn mul_is_compose_and_a_group() {
    let a = AffineTransform::scale(2.0, 3.0, Coord::new(0.0, 0.0));
    let b = AffineTransform::translate(10.0, -4.0);
    let p = Coord::new(1.5, -2.0);
    // `a * b` applies a then b — identical to `b.apply(a.apply(p))`.
    assert_eq!((a * b).apply(p), b.apply(a.apply(p)));
    assert_eq!((a * b).apply(p), a.compose(&b).apply(p));
    // By-reference Mul agrees (called explicitly: the values are `Copy`, so
    // the `&a * &b` sugar would trip `clippy::op_ref`).
    assert_eq!(std::ops::Mul::mul(&a, &b).apply(p), (a * b).apply(p));
    // Identity is the unit.
    assert_eq!((a * AffineTransform::identity()).apply(p), a.apply(p));
    assert_eq!((AffineTransform::identity() * a).apply(p), a.apply(p));
    // MulAssign matches Mul.
    let mut acc = a;
    acc *= b;
    assert_eq!(acc.apply(p), (a * b).apply(p));
    // Inverse round-trips to the identity.
    let inv = a.inverse().unwrap();
    let round = (a * inv).apply(p);
    assert!((round.x - p.x).abs() < EPS && (round.y - p.y).abs() < EPS);
}

#[test]
fn is_rigid_distinguishes_isometries() {
    // A rotation is rigid (orthonormal columns).
    assert!(AffineTransform::rotate(37.0, Coord::new(2.0, -1.0)).is_rigid());
    // A pure translation is rigid.
    assert!(AffineTransform::translate(5.0, -3.0).is_rigid());
    // The identity is rigid.
    assert!(AffineTransform::identity().is_rigid());
    // A reflection (negative uniform scale) is rigid (|det| = 1).
    assert!(AffineTransform::scale(-1.0, 1.0, Coord::new(0.0, 0.0)).is_rigid());
    // A non-uniform scale (x by 2) is NOT rigid.
    assert!(!AffineTransform::scale(2.0, 1.0, Coord::new(0.0, 0.0)).is_rigid());
    // A uniform scale ≠ 1 is not rigid (columns not unit length).
    assert!(!AffineTransform::scale(2.0, 2.0, Coord::new(0.0, 0.0)).is_rigid());
    // A shear is not rigid (columns not orthogonal).
    assert!(!AffineTransform::skew(30.0, 0.0, Coord::new(0.0, 0.0)).is_rigid());
    // A rotation composed with a translation stays rigid.
    let rt =
        AffineTransform::rotate(20.0, Coord::new(0.0, 0.0)) * AffineTransform::translate(7.0, 2.0);
    assert!(rt.is_rigid());
}

#[test]
fn empty_geometries_pass_through() {
    let empty = LineString::new(vec![]);
    assert_eq!(empty.scale(2.0), empty);
    assert_eq!(empty.skew(30.0), empty);
    assert_eq!(empty.rotate_around_centroid(45.0), empty);
    assert_eq!(empty.rotate_around_center(45.0), empty);
}
