use super::*;

#[test]
fn near_collinear_scaled_point_keeps_exact_sign() {
    let scale = 1_000_000_000_000_000_000_i128;
    let a = ExactCoord::new(0, 0);
    let b = ExactCoord::new(3 * scale, 3 * scale);
    let c = ExactCoord::new(scale, scale + 1);

    assert_eq!(orientation_determinant(a, b, c).unwrap(), 3 * scale);
    assert_eq!(orient(a, b, c).unwrap(), OrientationSign::CounterClockwise);
}

#[test]
fn incircle_classifies_inside_boundary_outside() {
    let a = ExactCoord::new(0, 0);
    let b = ExactCoord::new(2, 0);
    let c = ExactCoord::new(0, 2);

    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(1, 1)).unwrap(),
        InCirclePosition::Inside
    );
    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(2, 2)).unwrap(),
        InCirclePosition::OnCircle
    );
    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(3, 3)).unwrap(),
        InCirclePosition::Outside
    );
}

#[test]
fn incircle_position_is_independent_of_triangle_winding() {
    let a = ExactCoord::new(0, 0);
    let b = ExactCoord::new(2, 0);
    let c = ExactCoord::new(0, 2);
    let query = ExactCoord::new(1, 1);

    assert_eq!(
        in_circle_position(a, b, c, query).unwrap(),
        in_circle_position(c, b, a, query).unwrap()
    );
}

#[test]
fn segments_classify_cross_touch_overlap_and_disjoint() {
    assert_eq!(
        segment_intersection(
            ExactCoord::new(0, 0),
            ExactCoord::new(4, 4),
            ExactCoord::new(0, 4),
            ExactCoord::new(4, 0),
        )
        .unwrap(),
        ExactSegmentIntersection::ProperCross
    );
    assert!(
        segments_properly_cross(
            ExactCoord::new(0, 0),
            ExactCoord::new(4, 4),
            ExactCoord::new(0, 4),
            ExactCoord::new(4, 0),
        )
        .unwrap()
    );
    assert_eq!(
        segment_intersection(
            ExactCoord::new(0, 0),
            ExactCoord::new(2, 0),
            ExactCoord::new(2, 0),
            ExactCoord::new(2, 3),
        )
        .unwrap(),
        ExactSegmentIntersection::Touch
    );
    assert_eq!(
        segment_intersection(
            ExactCoord::new(0, 0),
            ExactCoord::new(4, 0),
            ExactCoord::new(2, 0),
            ExactCoord::new(6, 0),
        )
        .unwrap(),
        ExactSegmentIntersection::CollinearOverlap
    );
    assert_eq!(
        segment_intersection(
            ExactCoord::new(0, 0),
            ExactCoord::new(1, 0),
            ExactCoord::new(2, 0),
            ExactCoord::new(3, 0),
        )
        .unwrap(),
        ExactSegmentIntersection::Disjoint
    );
}

#[test]
fn point_in_triangle_classifies_boundary_and_winding() {
    let a = ExactCoord::new(0, 0);
    let b = ExactCoord::new(4, 0);
    let c = ExactCoord::new(0, 4);

    assert_eq!(
        point_in_triangle(a, b, c, ExactCoord::new(1, 1)).unwrap(),
        ExactTrianglePosition::Inside
    );
    assert_eq!(
        point_in_triangle(a, b, c, ExactCoord::new(2, 0)).unwrap(),
        ExactTrianglePosition::OnBoundary
    );
    assert_eq!(
        point_in_triangle(a, b, c, ExactCoord::new(3, 3)).unwrap(),
        ExactTrianglePosition::Outside
    );
    assert_eq!(
        point_in_triangle(c, b, a, ExactCoord::new(1, 1)).unwrap(),
        ExactTrianglePosition::Inside
    );
}
