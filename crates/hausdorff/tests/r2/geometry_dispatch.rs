use hausdorff::{
    Area, BoundingRect, Centroid, ConvexHull, Geometry, LineString, Point, Polygon, Simplify,
};

fn square() -> Geometry {
    Geometry::Polygon(Polygon::new(
        LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
        vec![],
    ))
}

#[test]
fn algorithms_dispatch_over_geometry_enum() {
    let g = square();
    assert_eq!(g.unsigned_area(), 16.0);
    let r = g.bounding_rect().unwrap();
    assert_eq!(
        (r.min().x, r.min().y, r.max().x, r.max().y),
        (0.0, 0.0, 4.0, 4.0)
    );
    let c = g.centroid().unwrap();
    assert!((c.x() - 2.0).abs() < 1e-9 && (c.y() - 2.0).abs() < 1e-9);
    assert_eq!(g.convex_hull().unsigned_area(), 16.0);

    let p = Geometry::Point(Point::new(3.0, 5.0));
    assert_eq!(p.unsigned_area(), 0.0);
    assert_eq!(
        p.bounding_rect().unwrap().min(),
        hausdorff::Coord::new(3.0, 5.0)
    );

    assert!(matches!(p.simplify(1.0), Geometry::Point(_)));
    let wiggly = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 0.01), (2.0, 0.0)]));
    let Geometry::LineString(s) = wiggly.simplify(0.1) else {
        panic!()
    };
    assert_eq!(s.coords_count(), 2);
}
