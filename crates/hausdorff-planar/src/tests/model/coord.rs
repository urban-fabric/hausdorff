use crate::space::r2::*;

#[test]
fn coord_ops_and_cross() {
    let a = Coord::new(1.0, 2.0);
    let b = Coord::new(3.0, 4.0);
    assert_eq!(a + b, Coord::new(4.0, 6.0));
    assert_eq!(b - a, Coord::new(2.0, 2.0));
    assert_eq!(a * 2.0, Coord::new(2.0, 4.0));
    let o = Coord::new(0.0, 0.0);
    assert!(Coord::new(1.0, 0.0).cross(o, Coord::new(0.0, 1.0)) > 0.0);
}

#[test]
fn point_coord_interop() {
    let p = Point::new(5.0, 6.0);
    assert_eq!(p.x_y(), (5.0, 6.0));
    let c: Coord = p.into();
    assert_eq!(c, Coord::new(5.0, 6.0));
    assert_eq!(Point::from(c), p);
}
