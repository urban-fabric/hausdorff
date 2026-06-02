//! Public-API tests for `hausdorff::space::r2`, relocated out of `src/`.

mod coord {
    use hausdorff::space::r2::*;

    #[test]
    fn coord_ops_and_cross() {
        let a = Coord::new(1.0, 2.0);
        let b = Coord::new(3.0, 4.0);
        assert_eq!(a + b, Coord::new(4.0, 6.0));
        assert_eq!(b - a, Coord::new(2.0, 2.0));
        assert_eq!(a * 2.0, Coord::new(2.0, 4.0));
        // CCW turn (origin -> +x -> +y) is positive.
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
}

mod primitives {
    use hausdorff::space::r2::*;

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
        assert_eq!(p.exterior().coords_count(), 5); // 4 corners + closing
    }
}

mod line_string {
    use hausdorff::space::r2::*;

    #[test]
    fn closing_a_ring() {
        let mut ls = LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]);
        assert!(!ls.is_closed());
        ls.close();
        assert!(ls.is_closed());
        assert_eq!(ls.coords_count(), 4);
        // already closed -> no-op
        ls.close();
        assert_eq!(ls.coords_count(), 4);
    }

    #[test]
    fn lines_iterates_segments() {
        let ls = LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]);
        let segs: Vec<_> = ls.lines().collect();
        assert_eq!(segs.len(), 2);
        assert_eq!(segs[0], Line::new((0.0, 0.0), (1.0, 0.0)));
    }
}

mod polygon {
    use hausdorff::space::r2::*;

    #[test]
    fn polygon_closes_rings() {
        let ext = LineString::new(vec![
            Coord::new(0.0, 0.0),
            Coord::new(1.0, 0.0),
            Coord::new(1.0, 1.0),
        ]);
        let p = Polygon::new(ext, vec![]);
        assert!(p.exterior().is_closed());
        assert_eq!(p.exterior().coords_count(), 4);
        assert!(p.interiors().is_empty());
    }

    #[test]
    fn exterior_mut_recloses() {
        let mut p = Polygon::new(
            LineString::new(vec![
                Coord::new(0.0, 0.0),
                Coord::new(2.0, 0.0),
                Coord::new(1.0, 2.0),
            ]),
            vec![],
        );
        p.exterior_mut(|ring| {
            ring.0.reverse();
        });
        assert!(p.exterior().is_closed());
    }
}

mod geometry {
    use hausdorff::space::r2::*;

    #[test]
    fn from_conversions() {
        let g: Geometry = Point::new(1.0, 2.0).into();
        assert!(matches!(g, Geometry::Point(_)));
        let g: Geometry = Polygon::default().into();
        assert!(matches!(g, Geometry::Polygon(_)));
    }
}

mod winding {
    use hausdorff::space::r2::*;

    fn ls(coords: &[(f64, f64)]) -> LineString {
        LineString::new(coords.iter().map(|&(x, y)| Coord::new(x, y)).collect())
    }

    #[test]
    fn winding_detection() {
        // CCW square
        let ccw = ls(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]);
        assert_eq!(ccw.winding_order(), Some(WindingOrder::CounterClockwise));
        assert!(ccw.is_ccw());

        let mut rev = ccw.0.clone();
        rev.reverse();
        assert_eq!(
            LineString::new(rev).winding_order(),
            Some(WindingOrder::Clockwise)
        );
    }

    #[test]
    fn orient_default_makes_exterior_ccw() {
        // Start with a CW exterior.
        let cw = Polygon::new(
            ls(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]),
            vec![],
        );
        assert!(cw.exterior().is_cw());
        let oriented = cw.orient(Direction::Default);
        assert!(oriented.exterior().is_ccw());
    }
}

mod area {
    use hausdorff::space::r2::*;

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
        assert_eq!(p.signed_area(), 1.0); // CCW -> positive
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
}

mod bounding_rect {
    use hausdorff::space::r2::*;

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
}

mod dimensions {
    use hausdorff::space::r2::*;

    #[test]
    fn dimensions_by_type() {
        assert_eq!(
            Point::new(0.0, 0.0).dimensions(),
            Dimensions::ZeroDimensional
        );
        assert_eq!(
            LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]).dimensions(),
            Dimensions::OneDimensional
        );
        assert_eq!(LineString::default().dimensions(), Dimensions::Empty);
        assert!(LineString::default().is_empty());
        assert_eq!(
            Polygon::new(
                LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]),
                vec![]
            )
            .dimensions(),
            Dimensions::TwoDimensional
        );
    }
}

mod centroid {
    use hausdorff::space::r2::*;

    fn poly(coords: &[(f64, f64)]) -> Polygon {
        Polygon::new(LineString::from(coords.to_vec()), vec![])
    }

    #[test]
    fn square_centroid_is_center() {
        let p = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
        let c = p.centroid().unwrap();
        assert!((c.x() - 1.0).abs() < 1e-12);
        assert!((c.y() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn centroid_is_winding_independent() {
        let ccw = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
        let cw = poly(&[(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)]);
        assert_eq!(ccw.centroid(), cw.centroid());
    }

    #[test]
    fn multipoint_centroid_is_mean() {
        let mp = MultiPoint::from(vec![(0.0, 0.0), (2.0, 0.0), (1.0, 3.0)]);
        let c = mp.centroid().unwrap();
        assert!((c.x() - 1.0).abs() < 1e-12);
        assert!((c.y() - 1.0).abs() < 1e-12);
    }
}

mod contains {
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
        assert!(s.contains(&Point::new(0.5, 0.5))); // inside exterior, outside hole
        assert!(!s.contains(&Point::new(2.0, 2.0))); // inside the hole
        assert!(!s.contains(&Point::new(9.0, 9.0))); // outside entirely
    }

    #[test]
    fn multipolygon_and_rect() {
        let mp = MultiPolygon::new(vec![square()]);
        assert!(mp.contains(&Coord::new(0.5, 0.5)));
        let r = Rect::new((0.0, 0.0), (2.0, 2.0));
        assert!(r.contains(&Point::new(1.0, 1.0)));
        assert!(!r.contains(&Point::new(2.0, 2.0))); // on boundary -> not contained
    }
}

mod convex_hull {
    use hausdorff::space::r2::*;

    #[test]
    fn hull_of_square_with_interior_point() {
        // Four corners plus a point inside; the hull is just the square.
        let mp = MultiPoint::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (2.0, 2.0), // interior — should be dropped
        ]);
        let hull = mp.convex_hull();
        assert_eq!(hull.unsigned_area(), 16.0);
        assert!(hull.exterior().is_closed());
        // 4 distinct corners + closing point
        assert_eq!(hull.exterior().coords_count(), 5);
    }

    #[test]
    fn hull_is_counter_clockwise() {
        let mp = MultiPoint::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        let hull = mp.convex_hull();
        assert!(hull.signed_area() > 0.0); // CCW
    }
}

mod macros {
    use hausdorff::space::r2::*;
    use hausdorff::{coord, line_string, point, polygon};

    #[test]
    fn coord_and_point_macros() {
        assert_eq!(coord! { x: 1.0, y: 2.0 }, Coord::new(1.0, 2.0));
        assert_eq!(point! { x: 1.0, y: 2.0 }, Point::new(1.0, 2.0));
    }

    #[test]
    fn line_string_macro_tuples() {
        let ls = line_string![(x: 0.0, y: 0.0), (x: 1.0, y: 1.0)];
        assert_eq!(ls.coords_count(), 2);
        assert_eq!(ls.0[1], Coord::new(1.0, 1.0));
    }

    #[test]
    fn polygon_macro_exterior_only() {
        let p = polygon![(x: 0.0, y: 0.0), (x: 1.0, y: 0.0), (x: 0.0, y: 1.0)];
        assert!(p.exterior().is_closed());
        assert_eq!(p.exterior().coords_count(), 4); // closed
        assert!(p.interiors().is_empty());
    }

    #[test]
    fn polygon_macro_with_hole() {
        let p = polygon!(
            exterior: [(x: 0.0, y: 0.0), (x: 4.0, y: 0.0), (x: 4.0, y: 4.0), (x: 0.0, y: 4.0)],
            interiors: [[(x: 1.0, y: 1.0), (x: 2.0, y: 1.0), (x: 2.0, y: 2.0), (x: 1.0, y: 2.0)]]
        );
        assert_eq!(p.interiors().len(), 1);
        assert!(p.interiors()[0].is_closed());
    }
}

mod geometry_dispatch {
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
        assert_eq!((r.min().x, r.min().y, r.max().x, r.max().y), (0.0, 0.0, 4.0, 4.0));
        let c = g.centroid().unwrap();
        assert!((c.x() - 2.0).abs() < 1e-9 && (c.y() - 2.0).abs() < 1e-9);
        assert_eq!(g.convex_hull().unsigned_area(), 16.0);

        // a point has zero area and bounds itself
        let p = Geometry::Point(Point::new(3.0, 5.0));
        assert_eq!(p.unsigned_area(), 0.0);
        assert_eq!(p.bounding_rect().unwrap().min(), hausdorff::Coord::new(3.0, 5.0));

        // simplify is a no-op on a point, real on a line string
        assert!(matches!(p.simplify(1.0), Geometry::Point(_)));
        let wiggly = Geometry::LineString(LineString::from(vec![
            (0.0, 0.0), (1.0, 0.01), (2.0, 0.0),
        ]));
        let Geometry::LineString(s) = wiggly.simplify(0.1) else { panic!() };
        assert_eq!(s.coords_count(), 2);
    }
}
