use hausdorff::space::r2::*;

#[test]
fn dimensions_by_type() {
    assert_eq!(
        Point::new(0.0, 0.0).dimensions(),
        Dimensions::ZeroDimensional
    );
    assert_eq!(
        LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]).dimensions(),
        Dimensions::OneDimensional
    );
    assert_eq!(LineString::default().dimensions(), Dimensions::Empty);
    assert!(LineString::default().is_empty());
    assert_eq!(
        Polygon::new(
            LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]),
            vec![]
        )
        .dimensions(),
        Dimensions::TwoDimensional
    );
}
