use hausdorff::space::r2::*;

fn square() -> Polygon {
    Polygon::new(
        LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
        vec![LineString::from(vec![
            (1.0, 1.0),
            (3.0, 1.0),
            (3.0, 3.0),
            (1.0, 3.0),
        ])],
    )
}

#[test]
fn point_in_polygon_respecting_hole() {
    let s = square();
    assert!(s.contains(&Point::new(0.5, 0.5)));
    assert!(!s.contains(&Point::new(2.0, 2.0)));
    assert!(!s.contains(&Point::new(9.0, 9.0)));
}

#[test]
fn multipolygon_and_rect() {
    let mp = MultiPolygon::new(vec![square()]);
    assert!(mp.contains(&Coord::new(0.5, 0.5)));
    let r = Rect::new((0.0, 0.0), (2.0, 2.0));
    assert!(r.contains(&Point::new(1.0, 1.0)));
    assert!(!r.contains(&Point::new(2.0, 2.0)));
}
