use crate::space::r2::*;

fn poly(coords: &[(f64, f64)]) -> Polygon {
    Polygon::new(
        LineString::new(coords.iter().map(|&(x, y)| Coord::new(x, y)).collect()),
        vec![],
    )
}

#[test]
fn unit_square_area() {
    let p = poly(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
    assert_eq!(p.unsigned_area(), 1.0);
    assert_eq!(p.signed_area(), 1.0);
}

#[test]
fn clockwise_is_negative_signed() {
    let p = poly(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]);
    assert_eq!(p.signed_area(), -1.0);
    assert_eq!(p.unsigned_area(), 1.0);
}

#[test]
fn hole_subtracts_area() {
    let p = Polygon::new(
        LineString::new(
            [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
                .iter()
                .map(|&(x, y)| Coord::new(x, y))
                .collect(),
        ),
        vec![LineString::new(
            [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]
                .iter()
                .map(|&(x, y)| Coord::new(x, y))
                .collect(),
        )],
    );
    assert_eq!(p.unsigned_area(), 16.0 - 1.0);
}

#[test]
fn rect_and_triangle() {
    assert_eq!(Rect::new((0.0, 0.0), (2.0, 3.0)).unsigned_area(), 6.0);
    let t = Triangle::new(
        Coord::new(0.0, 0.0),
        Coord::new(4.0, 0.0),
        Coord::new(0.0, 3.0),
    );
    assert_eq!(t.unsigned_area(), 6.0);
}
