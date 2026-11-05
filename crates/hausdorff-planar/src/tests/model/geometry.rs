use crate::space::r2::*;

#[test]
fn from_conversions() {
    let g: Geometry = Point::new(1.0, 2.0).into();
    assert!(matches!(g, Geometry::Point(_)));
    let g: Geometry = Polygon::default().into();
    assert!(matches!(g, Geometry::Polygon(_)));
}
