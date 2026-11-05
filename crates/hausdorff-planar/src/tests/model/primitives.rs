use crate::space::r2::*;

#[test]
fn line_length() {
    let l = Line::new((0.0, 0.0), (3.0, 4.0));
    assert_eq!(l.length(), 5.0);
    assert_eq!(l.delta(), Coord::new(3.0, 4.0));
}

#[test]
fn rect_normalizes_corners() {
    let r = Rect::new((2.0, 5.0), (-1.0, 1.0));
    assert_eq!(r.min(), Coord::new(-1.0, 1.0));
    assert_eq!(r.max(), Coord::new(2.0, 5.0));
    assert_eq!(r.width(), 3.0);
    assert_eq!(r.height(), 4.0);
    assert!(r.contains_coord(Coord::new(0.0, 2.0)));
    assert!(!r.contains_coord(Coord::new(9.0, 9.0)));
}

#[test]
fn rect_to_polygon_is_closed_quad() {
    let r = Rect::new((0.0, 0.0), (1.0, 1.0));
    let p = r.to_polygon();
    assert_eq!(p.exterior().coords_count(), 5);
}
