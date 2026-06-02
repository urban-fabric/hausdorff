//! Geometry-enum trait dispatch added for galois: `Length<Geometry>` (boundary
//! length / perimeter), `Intersects for Rect`, and `Orient for Geometry`.

use hausdorff::space::r2::{Direction, Euclidean, Intersects, Length, Orient, Winding, WindingOrder};
use hausdorff::{Coord, Geometry, LineString, Point, Polygon, Rect};

#[test]
fn length_for_geometry_sums_boundaries() {
    // Polygon → perimeter (a 2×2 square's exterior is 8).
    let square = Geometry::Polygon(Polygon::new(
        LineString::from(vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]),
        vec![],
    ));
    assert!((Euclidean.length(&square) - 8.0).abs() < 1e-9);

    // LineString → path length (3 + 4 = 7).
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (3.0, 0.0), (3.0, 4.0)]));
    assert!((Euclidean.length(&line) - 7.0).abs() < 1e-9);

    // Puntal → 0.
    assert_eq!(Euclidean.length(&Geometry::Point(Point::new(1.0, 1.0))), 0.0);
}

#[test]
fn rect_intersects_is_box_overlap() {
    let a = Rect::new(Coord::new(0.0, 0.0), Coord::new(2.0, 2.0));
    let overlap = Rect::new(Coord::new(1.0, 1.0), Coord::new(3.0, 3.0));
    let touch = Rect::new(Coord::new(2.0, 0.0), Coord::new(4.0, 2.0)); // shares edge x=2
    let disjoint = Rect::new(Coord::new(5.0, 5.0), Coord::new(6.0, 6.0));
    assert!(a.intersects(&overlap));
    assert!(a.intersects(&touch)); // touching counts
    assert!(!a.intersects(&disjoint));
}

#[test]
fn orient_for_geometry_reorients_areal_only() {
    // A clockwise exterior ring; Default orientation wants it counter-clockwise.
    let cw = Polygon::new(
        LineString::from(vec![(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)]),
        vec![],
    );
    match Geometry::Polygon(cw).orient(Direction::Default) {
        Geometry::Polygon(p) => {
            assert_eq!(p.exterior().winding_order(), Some(WindingOrder::CounterClockwise));
        }
        _ => panic!("polygon stays a polygon"),
    }

    // A non-areal geometry passes through unchanged.
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]));
    match line.orient(Direction::Default) {
        Geometry::LineString(ls) => assert_eq!(ls.0.len(), 2),
        _ => panic!("line stays a line"),
    }
}
