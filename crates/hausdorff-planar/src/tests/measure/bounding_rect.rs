use crate::space::r2::*;

#[test]
fn polygon_bounding_rect() {
    let p = Polygon::new(
        LineString::from(vec![(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (1.0, 5.0)]),
        vec![],
    );
    let r = p.bounding_rect().unwrap();
    assert_eq!(r.min(), Coord::new(0.0, 0.0));
    assert_eq!(r.max(), Coord::new(3.0, 5.0));
}

#[test]
fn point_bounding_rect_is_unconditional() {
    let r: Rect = Point::new(2.0, 4.0).bounding_rect();
    assert_eq!(r.min(), Coord::new(2.0, 4.0));
    assert_eq!(r.max(), Coord::new(2.0, 4.0));
}

#[test]
fn empty_line_string_has_no_rect() {
    assert!(LineString::default().bounding_rect().is_none());
}
