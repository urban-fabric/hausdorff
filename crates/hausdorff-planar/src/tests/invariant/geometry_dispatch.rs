//! Geometry-enum trait dispatch added for galois: `Length<Geometry>` (linear
//! length, curve parts only) and `Perimeter<Geometry>` (areal boundary),
//! `Intersects for Rect`, and `Orient for Geometry`.

use crate::space::r2::{
    Area, Centroid, Direction, Euclidean, Intersects, Length, Orient, Perimeter, Simplify, Winding,
    WindingOrder,
};
use crate::space::r2::{
    BoundingRect, ConvexHull, Coord, Geometry, LineString, MultiLineString, MultiPolygon, Point,
    Polygon, Rect,
};

/// A 4x4 square as a `Geometry`, for the accessor sweep below.
fn square() -> Geometry {
    Geometry::Polygon(Polygon::new(
        LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
        vec![],
    ))
}

#[test]
fn length_is_linear_and_perimeter_is_areal() {
    let square = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    assert_eq!(Euclidean.length(&square), 0.0);
    assert!((Euclidean.perimeter(&square) - 8.0).abs() < 1e-9);

    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (3.0, 0.0), (3.0, 4.0)]));
    assert!((Euclidean.length(&line) - 7.0).abs() < 1e-9);
    assert_eq!(Euclidean.perimeter(&line), 0.0);

    let pt = Geometry::Point(Point::new(1.0, 1.0));
    assert_eq!(Euclidean.length(&pt), 0.0);
    assert_eq!(Euclidean.perimeter(&pt), 0.0);

    let gc = Geometry::GeometryCollection(vec![
        Geometry::LineString(LineString::from(vec![(0.0, 0.0), (3.0, 0.0), (3.0, 4.0)])),
        Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (2.0, 0.0),
                (2.0, 2.0),
                (0.0, 2.0),
                (0.0, 0.0),
            ]),
            vec![],
        )),
    ]);
    assert!((Euclidean.length(&gc) - 7.0).abs() < 1e-9);
    assert!((Euclidean.perimeter(&gc) - 8.0).abs() < 1e-9);
}

#[test]
fn rect_intersects_is_box_overlap() {
    let a = Rect::new(Coord::new(0.0, 0.0), Coord::new(2.0, 2.0));
    let overlap = Rect::new(Coord::new(1.0, 1.0), Coord::new(3.0, 3.0));
    let touch = Rect::new(Coord::new(2.0, 0.0), Coord::new(4.0, 2.0));
    let disjoint = Rect::new(Coord::new(5.0, 5.0), Coord::new(6.0, 6.0));
    assert!(a.intersects(&overlap));
    assert!(a.intersects(&touch));
    assert!(!a.intersects(&disjoint));
}

#[test]
fn orient_for_geometry_reorients_areal_only() {
    let cw = Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (0.0, 2.0),
            (2.0, 2.0),
            (2.0, 0.0),
            (0.0, 0.0),
        ]),
        vec![],
    );
    match Geometry::Polygon(cw).orient(Direction::Default) {
        Geometry::Polygon(p) => {
            assert_eq!(
                p.exterior().winding_order(),
                Some(WindingOrder::CounterClockwise)
            );
        }
        _ => panic!("polygon stays a polygon"),
    }

    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]));
    match line.orient(Direction::Default) {
        Geometry::LineString(ls) => assert_eq!(ls.0.len(), 2),
        _ => panic!("line stays a line"),
    }
}

#[test]
fn collection_area_sums_absolute_parts() {
    let ccw = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    let cw = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (2.0, 0.0),
            (2.0, 1.0),
            (3.0, 1.0),
            (3.0, 0.0),
            (2.0, 0.0),
        ]),
        vec![],
    ));
    let gc = Geometry::GeometryCollection(vec![ccw, cw]);
    assert!(
        (gc.unsigned_area() - 2.0).abs() < 1e-9,
        "got {}",
        gc.unsigned_area()
    );
}

#[test]
fn collection_centroid_is_mean_of_all_coords() {
    let gc = Geometry::GeometryCollection(vec![
        Geometry::Point(Point::new(0.0, 0.0)),
        Geometry::Point(Point::new(4.0, 0.0)),
        Geometry::Point(Point::new(2.0, 6.0)),
    ]);
    let c = gc.centroid().unwrap();
    assert!(
        (c.0.x - 2.0).abs() < 1e-9 && (c.0.y - 2.0).abs() < 1e-9,
        "got {:?}",
        c.0
    );
}

#[test]
fn multilinestring_centroid_is_length_weighted() {
    let long = LineString::from(vec![(0.0, 0.0), (10.0, 0.0)]);
    let short = LineString::from(vec![(0.0, 10.0), (1.0, 10.0)]);
    let mls = Geometry::MultiLineString(MultiLineString::new(vec![long, short]));
    let c = mls.centroid().unwrap();
    assert!(
        (c.0.y - 10.0 / 11.0).abs() < 1e-9,
        "weighted y={}, not the mean 5",
        c.0.y
    );
}

#[test]
fn simplify_recurses_into_collection() {
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

#[test]
fn geometry_accessors_and_classifiers() {
    let poly = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]),
        vec![LineString::from(vec![
            (1.0, 1.0),
            (2.0, 1.0),
            (2.0, 2.0),
            (1.0, 2.0),
            (1.0, 1.0),
        ])],
    ));
    assert_eq!(poly.exterior_ring().unwrap().0.len(), 5);
    assert_eq!(poly.num_interior_rings(), Some(1));
    assert_eq!(poly.interior_ring(0).unwrap().0.len(), 5);
    assert!(poly.interior_ring(1).is_none());
    assert_eq!(poly.geometry_type(), "Polygon");
    assert_eq!(poly.num_geometries(), 1);
    assert_eq!(poly.is_polygon_ccw(), Some(true));
    assert_eq!(poly.is_polygon_cw(), Some(false));
    assert_eq!(poly.as_multi_polygon().unwrap().0.len(), 1);
    assert!(poly.point_n(0).is_none());

    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]));
    assert_eq!(line.point_n(0), Some(Point::new(0.0, 0.0)));
    assert_eq!(line.point_n(2), Some(Point::new(2.0, 0.0)));
    assert!(line.point_n(3).is_none());
    assert_eq!(line.start_point(), Some(Point::new(0.0, 0.0)));
    assert_eq!(line.end_point(), Some(Point::new(2.0, 0.0)));
    assert_eq!(line.is_closed(), Some(false));
    assert_eq!(line.exterior_ring(), None);
    assert!(line.as_multi_polygon().is_none());

    let ring = Geometry::LineString(LineString::from(vec![
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 0.0),
    ]));
    assert_eq!(ring.is_closed(), Some(true));

    let mp = Geometry::MultiPolygon(MultiPolygon::new(vec![
        Polygon::new(
            LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]),
            vec![],
        ),
        Polygon::new(
            LineString::from(vec![(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 2.0)]),
            vec![],
        ),
    ]));
    assert_eq!(mp.geometry_type(), "MultiPolygon");
    assert_eq!(mp.num_geometries(), 2);
    assert!(matches!(mp.geometry_n(1), Some(Geometry::Polygon(_))));
    assert!(mp.geometry_n(2).is_none());
    assert_eq!(mp.as_multi_polygon().unwrap().0.len(), 2);

    let pt = Geometry::Point(Point::new(5.0, 5.0));
    assert_eq!(pt.geometry_type(), "Point");
    assert_eq!(pt.num_geometries(), 1);
    assert!(matches!(pt.geometry_n(0), Some(Geometry::Point(_))));
    assert!(pt.geometry_n(1).is_none());
    assert!(pt.is_closed().is_none());
    assert!(pt.is_polygon_cw().is_none());
}

// Folded in from `r2/geometry_dispatch.rs`, which held this single test and a
// second `square()` helper. Two files named the same thing, covering the same
// question, is how the duplicate went unnoticed.
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
        crate::space::r2::Coord::new(3.0, 5.0)
    );

    assert!(matches!(p.simplify(1.0), Geometry::Point(_)));
    let wiggly = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 0.01), (2.0, 0.0)]));
    let Geometry::LineString(s) = wiggly.simplify(0.1) else {
        panic!()
    };
    assert_eq!(s.coords_count(), 2);
}
