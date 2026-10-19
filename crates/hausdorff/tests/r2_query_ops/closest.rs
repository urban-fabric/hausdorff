use super::*;

#[test]
fn closest_point_on_line_matches_reference() {
    let segments = [
        ((0.0, 0.0), (4.0, 0.0)),
        ((0.0, 0.0), (4.0, 4.0)),
        ((-2.0, 1.0), (3.0, -1.0)),
        ((1.0, 1.0), (1.0, 1.0)),
    ];
    for (start, end) in segments {
        let h = Line::new(Coord::new(start.0, start.1), Coord::new(end.0, end.1));
        for (x, y) in query_points() {
            assert_closest_eq(
                h.closest_point(&Point::new(x, y)),
                ref_closest_segment((x, y), start, end),
                &format!("line {start:?}-{end:?} vs ({x}, {y})"),
            );
        }
    }
}

#[test]
fn closest_point_on_point_and_multipoint_matches_reference() {
    let h_pt = Point::new(2.0, 2.0);
    let r_pt = (2.0, 2.0);
    let pts = vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
    let h_mp = MultiPoint::from(pts.clone());
    for (x, y) in query_points() {
        let q = (x, y);
        let p = Point::new(x, y);
        // Single point: intersection iff it is the query, else that point.
        let r_single = if q == r_pt {
            RefClosest::Intersection(r_pt.0, r_pt.1)
        } else {
            RefClosest::SinglePoint(r_pt.0, r_pt.1)
        };
        assert_closest_eq(
            h_pt.closest_point(&p),
            r_single,
            &format!("point vs ({x}, {y})"),
        );
        // Multipoint: fold over each member's single-point result.
        let r_mp = ref_fold(
            pts.iter().map(|&pt| {
                if q == pt {
                    RefClosest::Intersection(pt.0, pt.1)
                } else {
                    RefClosest::SinglePoint(pt.0, pt.1)
                }
            }),
            q,
        );
        assert_closest_eq(
            h_mp.closest_point(&p),
            r_mp,
            &format!("multipoint vs ({x}, {y})"),
        );
    }
    assert_closest_eq(
        MultiPoint::new(vec![]).closest_point(&h_pt),
        RefClosest::Indeterminate,
        "empty multipoint",
    );
}

#[test]
fn closest_point_on_line_string_matches_reference() {
    let paths = [
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        vec![(-3.0, -2.0), (1.0, 3.0), (5.0, -1.0), (8.0, 6.0)],
        vec![],
    ];
    for path in &paths {
        let h = h_line_string(path);
        for (x, y) in query_points() {
            assert_closest_eq(
                h.closest_point(&Point::new(x, y)),
                ref_closest_path((x, y), path),
                &format!("linestring {path:?} vs ({x}, {y})"),
            );
        }
    }
}

#[test]
fn closest_point_on_polygon_matches_reference() {
    let exterior = vec![(0.0, 0.0), (8.0, 0.0), (8.0, 8.0), (0.0, 8.0)];
    let hole = vec![(3.0, 3.0), (5.0, 3.0), (5.0, 5.0), (3.0, 5.0)];
    let cases: Vec<PolyFixture> = vec![
        (h_polygon(&exterior, &[]), exterior.clone(), vec![]),
        (
            h_polygon(&exterior, std::slice::from_ref(&hole)),
            exterior.clone(),
            vec![hole.clone()],
        ),
        (h_polygon(&[], &[]), vec![], vec![]),
    ];
    let queries = [
        (1.0, 1.0),
        (0.0, 0.0),
        (4.0, 0.0),
        (8.0, 8.0),
        (4.0, 3.5),
        (3.0, 4.0),
        (-2.0, 4.0),
        (12.0, 12.0),
        (4.0, -5.0),
    ];
    for (hp, ext, holes) in &cases {
        for (x, y) in queries {
            assert_closest_eq(
                hp.closest_point(&Point::new(x, y)),
                ref_closest_polygon((x, y), ext, holes),
                &format!("polygon vs ({x}, {y})"),
            );
        }
    }
}

#[test]
fn closest_point_on_multi_geometries_matches_reference() {
    let part_a = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
    let part_b = vec![(6.0, 6.0), (9.0, 6.0), (9.0, 9.0), (6.0, 9.0)];
    let h_mp = hausdorff::MultiPolygon::new(vec![h_polygon(&part_a, &[]), h_polygon(&part_b, &[])]);

    let line_a = vec![(0.0, 0.0), (4.0, 0.0)];
    let line_b = vec![(0.0, 6.0), (4.0, 6.0)];
    let h_mls =
        hausdorff::MultiLineString::new(vec![h_line_string(&line_a), h_line_string(&line_b)]);

    for (x, y) in query_points() {
        let q = (x, y);
        let p = Point::new(x, y);
        let r_mp = ref_fold(
            [
                ref_closest_polygon(q, &part_a, &[]),
                ref_closest_polygon(q, &part_b, &[]),
            ],
            q,
        );
        assert_closest_eq(
            h_mp.closest_point(&p),
            r_mp,
            &format!("multipolygon vs ({x}, {y})"),
        );
        let r_mls = ref_fold(
            [ref_closest_path(q, &line_a), ref_closest_path(q, &line_b)],
            q,
        );
        assert_closest_eq(
            h_mls.closest_point(&p),
            r_mls,
            &format!("multilinestring vs ({x}, {y})"),
        );
    }
    assert_closest_eq(
        hausdorff::MultiPolygon::new(vec![]).closest_point(&Point::new(1.0, 1.0)),
        RefClosest::Indeterminate,
        "empty multipolygon",
    );
}

#[test]
fn closest_point_on_rect_triangle_and_geometry_matches_reference() {
    let h_rect = Rect::new(Coord::new(1.0, 1.0), Coord::new(5.0, 3.0));
    let rect_ring = vec![(1.0, 1.0), (5.0, 1.0), (5.0, 3.0), (1.0, 3.0), (1.0, 1.0)];
    let h_tri = Triangle::new(
        Coord::new(0.0, 0.0),
        Coord::new(6.0, 0.0),
        Coord::new(0.0, 6.0),
    );
    let tri_ring = vec![(0.0, 0.0), (6.0, 0.0), (0.0, 6.0), (0.0, 0.0)];
    for (x, y) in query_points() {
        let q = (x, y);
        let p = Point::new(x, y);
        // Rect and Triangle are filled areas: inside-or-on is an intersection,
        // otherwise the nearest point on the boundary ring.
        let r_rect = ref_closest_polygon(q, &rect_ring, &[]);
        assert_closest_eq(
            h_rect.closest_point(&p),
            r_rect,
            &format!("rect vs ({x}, {y})"),
        );
        let r_tri = ref_closest_polygon(q, &tri_ring, &[]);
        assert_closest_eq(
            h_tri.closest_point(&p),
            r_tri,
            &format!("triangle vs ({x}, {y})"),
        );
        assert_closest_eq(
            Geometry::Triangle(h_tri).closest_point(&p),
            r_tri,
            &format!("geometry vs ({x}, {y})"),
        );
    }
}

// ---------------------------------------------------------------------------
// Independent linear-referencing reference oracle.
// ---------------------------------------------------------------------------
