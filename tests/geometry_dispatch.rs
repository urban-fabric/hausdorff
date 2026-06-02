//! Geometry-enum trait dispatch added for galois: `Length<Geometry>` (boundary
//! length / perimeter), `Intersects for Rect`, and `Orient for Geometry`.

use hausdorff::space::r2::{
    Area, Centroid, Direction, Euclidean, Intersects, Length, Orient, Simplify, Winding,
    WindingOrder,
};
use hausdorff::{Coord, Geometry, LineString, MultiLineString, Point, Polygon, Rect};

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

// ---- refinements (galois-correct semantics for collections) ----------------

#[test]
fn collection_area_sums_absolute_parts() {
    // Two unit squares of opposite winding: `|Σ signed|` would cancel to 0, but
    // the unsigned area must be 2.
    let ccw = Geometry::Polygon(Polygon::new(
        LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]),
        vec![],
    ));
    let cw = Geometry::Polygon(Polygon::new(
        LineString::from(vec![(2.0, 0.0), (2.0, 1.0), (3.0, 1.0), (3.0, 0.0), (2.0, 0.0)]),
        vec![],
    ));
    let gc = Geometry::GeometryCollection(vec![ccw, cw]);
    assert!((gc.unsigned_area() - 2.0).abs() < 1e-9, "got {}", gc.unsigned_area());
}

#[test]
fn collection_centroid_is_mean_of_all_coords() {
    let gc = Geometry::GeometryCollection(vec![
        Geometry::Point(Point::new(0.0, 0.0)),
        Geometry::Point(Point::new(4.0, 0.0)),
        Geometry::Point(Point::new(2.0, 6.0)),
    ]);
    let c = gc.centroid().unwrap();
    assert!((c.0.x - 2.0).abs() < 1e-9 && (c.0.y - 2.0).abs() < 1e-9, "got {:?}", c.0);
}

#[test]
fn multilinestring_centroid_is_length_weighted() {
    // A long line on y=0 (length 10) and a short line on y=10 (length 1): the
    // centroid's y is pulled toward 0, not the unweighted midpoint 5.
    let long = LineString::from(vec![(0.0, 0.0), (10.0, 0.0)]);
    let short = LineString::from(vec![(0.0, 10.0), (1.0, 10.0)]);
    let mls = Geometry::MultiLineString(MultiLineString::new(vec![long, short]));
    let c = mls.centroid().unwrap();
    assert!((c.0.y - 10.0 / 11.0).abs() < 1e-9, "weighted y={}, not the mean 5", c.0.y);
}

#[test]
fn simplify_recurses_into_collection() {
    // A line with a redundant collinear midpoint, wrapped in a collection: the
    // member must be simplified (the old dispatch cloned collections untouched).
    let line = LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
    let gc = Geometry::GeometryCollection(vec![Geometry::LineString(line)]);
    match gc.simplify(0.1) {
        Geometry::GeometryCollection(members) => match &members[0] {
            Geometry::LineString(ls) => assert_eq!(ls.0.len(), 2, "midpoint should drop"),
            _ => panic!("member stays a line"),
        },
        _ => panic!("stays a collection"),
    }
}
