use super::*;

#[test]
fn orient2d_signs_and_masks() {
    assert_eq!(orient2d(0.0, 0.0, 1.0, 0.0, 0.0, 1.0), 1.0);
    assert_eq!(orient2d(0.0, 0.0, 0.0, 1.0, 1.0, 0.0), -1.0);
    assert!(turns_positive2d::<f64>(0.0, 0.0, 1.0, 0.0, 0.0, 1.0));
    assert!(!turns_positive2d::<f64>(0.0, 0.0, 0.0, 1.0, 1.0, 0.0));
    assert!(collinear2d::<f64>(0.0, 0.0, 2.0, 2.0, 4.0, 4.0, 0.0));
}

#[test]
fn triangle_kernel_uses_mask_algebra() {
    assert!(in_triangle2d::<f64>(0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 1.0, 1.0));
    assert!(!in_triangle2d::<f64>(
        0.0, 0.0, 4.0, 0.0, 0.0, 4.0, -1.0, 1.0
    ));
    assert!(!in_triangle2d::<f64>(
        0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 2.0, 0.0
    ));

    assert!(in_triangle2d::<f64>(0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 1.0, 1.0));
}

#[test]
fn rect_circle_and_ellipse_predicates_are_mask_valued() {
    assert!(rect_contains_xy::<f64>(0.0, 0.0, 2.0, 3.0, 1.0, 3.0));
    assert!(!rect_contains_xy::<f64>(0.0, 0.0, 2.0, 3.0, 3.0, 1.0));
    assert!(rects_overlap_xy::<f64>(
        0.0, 0.0, 2.0, 2.0, 2.0, 1.0, 3.0, 3.0
    ));
    assert!(!rects_overlap_xy::<f64>(
        0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0
    ));

    assert!(circle_contains_xy::<f64>(0.0, 0.0, 2.0, 1.0, 1.0));
    assert!(!circle_contains_xy::<f64>(0.0, 0.0, 2.0, 3.0, 0.0));

    assert!(ellipse_contains_xy::<f64>(0.0, 0.0, 4.0, 2.0, 2.0, 0.5));
    assert!(!ellipse_contains_xy::<f64>(0.0, 0.0, 4.0, 2.0, 5.0, 0.0));
    assert!(!ellipse_contains_xy::<f64>(0.0, 0.0, 0.0, 2.0, 0.0, 0.0));
}

#[test]
fn segment_predicates_are_mask_valued() {
    assert!(point_on_segment2d::<f64>(0.0, 0.0, 4.0, 0.0, 2.0, 0.0));
    assert!(!point_on_segment2d::<f64>(0.0, 0.0, 4.0, 0.0, 5.0, 0.0));
    assert!(!point_on_segment2d::<f64>(0.0, 0.0, 4.0, 0.0, 2.0, 1.0));

    assert!(segments_intersect2d::<f64>(
        0.0, 0.0, 4.0, 4.0, 0.0, 4.0, 4.0, 0.0
    ));
    assert!(segments_intersect2d::<f64>(
        0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 6.0, 0.0
    ));
    assert!(!segments_intersect2d::<f64>(
        0.0, 0.0, 4.0, 0.0, 5.0, 0.0, 6.0, 0.0
    ));
    assert!(!segments_intersect2d::<f64>(
        0.0, 0.0, 4.0, 0.0, 0.0, 1.0, 4.0, 1.0
    ));
}

#[test]
fn distance_kernels_match_known_values() {
    assert_eq!(distance_squared2d(0.0, 0.0, 3.0, 4.0), 25.0);
    assert_eq!(
        rect_distance_squared_xy(0.0, 0.0, 1.0, 1.0, 3.0, 4.0, 5.0, 6.0),
        13.0
    );
    assert_eq!(
        point_segment_distance_squared2d(2.0, 3.0, 0.0, 0.0, 4.0, 0.0),
        9.0
    );
    assert_eq!(
        point_segment_distance_squared2d(3.0, 4.0, 1.0, 1.0, 1.0, 1.0),
        13.0
    );
    assert_eq!(
        segment_distance_squared2d(0.0, 0.0, 4.0, 0.0, 0.0, 1.0, 4.0, 1.0),
        1.0
    );
    assert_eq!(
        segment_distance_squared2d(0.0, 0.0, 4.0, 4.0, 0.0, 4.0, 4.0, 0.0),
        0.0
    );
}

#[test]
fn kernels_run_for_f32() {
    let orient = orient2d(0.0_f32, 0.0, 1.0, 0.0, 0.0, 1.0);
    assert_eq!(orient, 1.0);
    assert!(rect_contains_xy(0.0_f32, 0.0, 2.0, 2.0, 1.0, 1.0));
    assert!(ellipse_contains_xy(0.0_f32, 0.0, 3.0, 2.0, 1.0, 1.0));
    assert!(segments_intersect2d(
        0.0_f32, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0
    ));
    assert_eq!(
        point_segment_distance_squared2d(2.0_f32, 3.0, 0.0, 0.0, 4.0, 0.0),
        9.0
    );
}
