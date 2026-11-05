use crate::space::r2::{
    Closest, ClosestPoint, Coord, Geometry, Line, MultiPoint, Point, Rect, Triangle,
};
use crate::tests::support::{
    PolyFixture, assert_points_eq, closed, dist2, line_string, polygon, query_points,
};

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
        let h = line_string(path);
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
        (polygon(&exterior, &[]), exterior.clone(), vec![]),
        (
            polygon(&exterior, std::slice::from_ref(&hole)),
            exterior.clone(),
            vec![hole.clone()],
        ),
        (polygon(&[], &[]), vec![], vec![]),
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
    let h_mp =
        crate::space::r2::MultiPolygon::new(vec![polygon(&part_a, &[]), polygon(&part_b, &[])]);

    let line_a = vec![(0.0, 0.0), (4.0, 0.0)];
    let line_b = vec![(0.0, 6.0), (4.0, 6.0)];
    let h_mls =
        crate::space::r2::MultiLineString::new(vec![line_string(&line_a), line_string(&line_b)]);

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
        crate::space::r2::MultiPolygon::new(vec![]).closest_point(&Point::new(1.0, 1.0)),
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

// ---- reference oracle ----
//
// A self-contained closest-point implementation: closed-form point-to-segment
// projection plus a ray-cast point-in-polygon. It lives here rather than in
// `support` because no other suite can use it, and keeping it beside the tests it
// checks is what makes the comparison auditable.

#[derive(Clone, Copy, Debug)]
enum RefClosest {
    Intersection(f64, f64),
    SinglePoint(f64, f64),
    Indeterminate,
}

fn on_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> bool {
    let cross = (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0);
    cross == 0.0
        && p.0 <= a.0.max(b.0)
        && p.0 >= a.0.min(b.0)
        && p.1 <= a.1.max(b.1)
        && p.1 >= a.1.min(b.1)
}

/// Closest point of a single (possibly degenerate) segment to `p`.
fn ref_closest_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> RefClosest {
    let (vx, vy) = (b.0 - a.0, b.1 - a.1);
    let v_sq = vx * vx + vy * vy;
    if v_sq == 0.0 {
        return RefClosest::Indeterminate;
    }
    let t = ((p.0 - a.0) * vx + (p.1 - a.1) * vy) / v_sq;
    if t < 0.0 {
        return RefClosest::SinglePoint(a.0, a.1);
    }
    if t > 1.0 {
        return RefClosest::SinglePoint(b.0, b.1);
    }
    let c = (a.0 + vx * t, a.1 + vy * t);
    if on_segment(p, a, b) {
        RefClosest::Intersection(c.0, c.1)
    } else {
        RefClosest::SinglePoint(c.0, c.1)
    }
}

/// Mirrors `Closest::best_of_two`: intersection wins outright, indeterminate
/// loses to anything, otherwise the nearer single point wins (left on ties).
fn ref_best(left: RefClosest, right: RefClosest, p: (f64, f64)) -> RefClosest {
    let l = match left {
        RefClosest::Indeterminate => return right,
        RefClosest::Intersection(..) => return left,
        RefClosest::SinglePoint(x, y) => (x, y),
    };
    let r = match right {
        RefClosest::Indeterminate => return left,
        RefClosest::Intersection(..) => return right,
        RefClosest::SinglePoint(x, y) => (x, y),
    };
    if dist2(l, p) <= dist2(r, p) {
        left
    } else {
        right
    }
}

/// Fold sub-results, short-circuiting on the first intersection.
fn ref_fold(items: impl IntoIterator<Item = RefClosest>, p: (f64, f64)) -> RefClosest {
    let mut best = RefClosest::Indeterminate;
    for item in items {
        best = ref_best(item, best, p);
        if matches!(best, RefClosest::Intersection(..)) {
            return best;
        }
    }
    best
}

/// Closest result over a polyline's segments.
fn ref_closest_path(p: (f64, f64), path: &[(f64, f64)]) -> RefClosest {
    ref_fold(
        path.windows(2).map(|w| ref_closest_segment(p, w[0], w[1])),
        p,
    )
}

/// Ray-cast (even-odd) test: is `p` strictly inside the ring? Boundary handling
/// is done separately, so a point exactly on an edge may report either way here.
fn point_in_ring(p: (f64, f64), ring: &[(f64, f64)]) -> bool {
    let mut inside = false;
    for w in ring.windows(2) {
        let (x1, y1) = w[0];
        let (x2, y2) = w[1];
        let crosses = (y1 > p.1) != (y2 > p.1);
        if crosses {
            let x_int = x1 + (p.1 - y1) / (y2 - y1) * (x2 - x1);
            if p.0 < x_int {
                inside = !inside;
            }
        }
    }
    inside
}

/// Append the closing vertex if `ring` is not already closed.
fn on_ring(p: (f64, f64), ring: &[(f64, f64)]) -> bool {
    ring.windows(2).any(|w| on_segment(p, w[0], w[1]))
}

/// Closest result for a polygon (exterior + holes), matching the kernel: a
/// point inside the filled area or on any ring is an intersection, otherwise
/// the nearest point on any ring.
fn ref_closest_polygon(
    p: (f64, f64),
    exterior: &[(f64, f64)],
    holes: &[Vec<(f64, f64)>],
) -> RefClosest {
    if exterior.is_empty() {
        return ref_closest_path(p, exterior);
    }
    let ext = closed(exterior);
    let hole_rings: Vec<Vec<(f64, f64)>> = holes.iter().map(|h| closed(h)).collect();
    let interior = point_in_ring(p, &ext) && !hole_rings.iter().any(|h| point_in_ring(p, h));
    let boundary = on_ring(p, &ext) || hole_rings.iter().any(|h| on_ring(p, h));
    if interior || boundary {
        return RefClosest::Intersection(p.0, p.1);
    }
    // Kernel order: holes first, then exterior (only matters for SinglePoint
    // tie-breaks, which the squared-distance comparison resolves identically).
    let mut items: Vec<RefClosest> = Vec::new();
    for h in &hole_rings {
        items.push(ref_closest_path(p, h));
    }
    items.push(ref_closest_path(p, &ext));
    ref_fold(items, p)
}

fn assert_closest_eq(h: Closest, r: RefClosest, ctx: &str) {
    match (h, r) {
        (Closest::Intersection(hp), RefClosest::Intersection(x, y)) => {
            assert_points_eq(hp, (x, y), ctx)
        }
        (Closest::SinglePoint(hp), RefClosest::SinglePoint(x, y)) => {
            assert_points_eq(hp, (x, y), ctx)
        }
        (Closest::Indeterminate, RefClosest::Indeterminate) => {}
        (h, r) => panic!("{ctx}: variant mismatch {h:?} vs ref {r:?}"),
    }
}
