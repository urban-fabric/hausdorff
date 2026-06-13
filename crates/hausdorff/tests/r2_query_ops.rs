//! Conformance tests for the planar query operations added alongside the core
//! algorithm surface — `ClosestPoint`, linear referencing
//! (`LineInterpolatePoint` / `LineLocatePoint`), the affine transform stack
//! (`AffineTransform` + `AffineOps` / `Translate` / `Rotate` / `Scale` /
//! `Skew`), and `Extremes`.
//!
//! These tests have **no third-party geometry dependency**. Every hausdorff
//! result is cross-checked against a small, self-contained reference oracle
//! written inline below: closed-form point-to-segment projection and ray-cast
//! point-in-polygon for `ClosestPoint`, a cumulative-length walk for the linear
//! referencing operations, and direct 3×3 matrix algebra for the affine stack.
//! Where a convenience transform pivots about a geometry's centroid or
//! bounding-box centre, we anchor on hausdorff's own (independently verified)
//! `Centroid` / `BoundingRect` kernels and apply the matrix with our own math —
//! a legitimate cross-subsystem check.

use hausdorff::space::r2::{
    AffineOps, AffineTransform, Centroid, Closest, ClosestPoint, Extremes, LineInterpolatePoint,
    LineLocatePoint, Rotate, Scale, Skew, Translate,
};
use hausdorff::{Coord, Geometry, Line, LineString, MultiPoint, Point, Polygon, Rect, Triangle};

const EPS: f64 = 1e-9;

fn h_line_string(coords: &[(f64, f64)]) -> LineString {
    LineString::new(coords.iter().map(|&(x, y)| Coord::new(x, y)).collect())
}

fn h_polygon(exterior: &[(f64, f64)], holes: &[Vec<(f64, f64)>]) -> Polygon {
    Polygon::new(
        h_line_string(exterior),
        holes.iter().map(|h| h_line_string(h)).collect(),
    )
}

fn assert_points_eq(h: Point, r: (f64, f64), ctx: &str) {
    assert!(
        (h.x() - r.0).abs() < EPS && (h.y() - r.1).abs() < EPS,
        "{ctx}: ({}, {}) vs ref ({}, {})",
        h.x(),
        h.y(),
        r.0,
        r.1
    );
}

// ---------------------------------------------------------------------------
// Independent `ClosestPoint` reference oracle.
//
// `RefClosest` mirrors the documented three-way result: the query lies *on* the
// geometry (`Intersection`), a unique nearest spot exists (`SinglePoint`), or
// the geometry is degenerate / empty (`Indeterminate`). Everything below is
// computed from first principles — segment projection and ray-cast
// point-in-ring — with no dependence on hausdorff's own closest-point code.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
enum RefClosest {
    Intersection(f64, f64),
    SinglePoint(f64, f64),
    Indeterminate,
}

fn dist2(a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dy) = (a.0 - b.0, a.1 - b.1);
    dx * dx + dy * dy
}

/// True iff `p` lies on the bounded segment `a`–`b` (collinear, within bbox),
/// matching the kernel's exact-zero cross-product test.
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
fn closed(ring: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut v = ring.to_vec();
    if let (Some(&first), Some(&last)) = (v.first(), v.last()) {
        if first != last {
            v.push(first);
        }
    }
    v
}

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
    let interior =
        point_in_ring(p, &ext) && !hole_rings.iter().any(|h| point_in_ring(p, h));
    let boundary =
        on_ring(p, &ext) || hole_rings.iter().any(|h| on_ring(p, h));
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

fn query_points() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0),
        (4.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (5.0, 5.0),
        (-3.0, 1.0),
        (10.0, -10.0),
        (0.5, 3.75),
        (4.0, 4.0),
        (100.0, 0.0),
    ]
}

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
    let cases: Vec<(Polygon, Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>)> = vec![
        (h_polygon(&exterior, &[]), exterior.clone(), vec![]),
        (
            h_polygon(&exterior, &[hole.clone()]),
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
    let rect_ring = vec![
        (1.0, 1.0),
        (5.0, 1.0),
        (5.0, 3.0),
        (1.0, 3.0),
        (1.0, 1.0),
    ];
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

fn seg_len(a: (f64, f64), b: (f64, f64)) -> f64 {
    dist2(a, b).sqrt()
}

fn path_len(path: &[(f64, f64)]) -> f64 {
    path.windows(2).map(|w| seg_len(w[0], w[1])).sum()
}

/// Reference for `Line::line_interpolate_point` / `LineString` over a polyline.
/// Returns `None` for `NaN` fractions or non-finite results, clamps out-of-range
/// fractions to the endpoints, and returns `None` for empty paths in range.
fn ref_interpolate(path: &[(f64, f64)], fraction: f64) -> Option<(f64, f64)> {
    if fraction.is_nan() {
        return None;
    }
    let f = fraction.clamp(0.0, 1.0);
    let total = path_len(path);
    let target = total * f;
    let mut cum = 0.0;
    for w in path.windows(2) {
        let (a, b) = (w[0], w[1]);
        let len = seg_len(a, b);
        if cum + len >= target {
            // Mirror the kernel exactly: a zero-length segment yields a NaN
            // segment-fraction (0/0), and interpolating that gives None.
            let sf = (target - cum) / len;
            return ref_interpolate_seg(a, b, sf);
        }
        cum += len;
    }
    None
}

fn ref_interpolate_seg(a: (f64, f64), b: (f64, f64), fraction: f64) -> Option<(f64, f64)> {
    if fraction.is_nan() {
        return None;
    }
    let f = fraction.clamp(0.0, 1.0);
    let r = (a.0 + (b.0 - a.0) * f, a.1 + (b.1 - a.1) * f);
    (r.0.is_finite() && r.1.is_finite()).then_some(r)
}

/// Perpendicular distance from `p` to the bounded segment `a`–`b`.
fn point_seg_dist(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (abx, aby) = (b.0 - a.0, b.1 - a.1);
    let len2 = abx * abx + aby * aby;
    let proj = if len2 == 0.0 {
        a
    } else {
        let t = ((p.0 - a.0) * abx + (p.1 - a.1) * aby) / len2;
        let t = t.clamp(0.0, 1.0);
        (a.0 + abx * t, a.1 + aby * t)
    };
    seg_len(p, proj)
}

/// Reference for `Line::line_locate_point`: fraction along the segment of the
/// clamped projection; `None` on non-finite input, `Some(0)` for zero-length.
fn ref_locate_seg(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> Option<f64> {
    let (vx, vy) = (b.0 - a.0, b.1 - a.1);
    let v_sq = vx * vx + vy * vy;
    if v_sq == 0.0 {
        return Some(0.0);
    }
    let l = ((p.0 - a.0) * vx + (p.1 - a.1) * vy) / v_sq;
    l.is_finite().then(|| l.clamp(0.0, 1.0))
}

/// Reference for `LineString::line_locate_point`: the fraction of total length
/// of the closest location along the polyline.
fn ref_locate_path(path: &[(f64, f64)], p: (f64, f64)) -> Option<f64> {
    let total = path_len(path);
    if total == 0.0 {
        return Some(0.0);
    }
    let mut cum = 0.0;
    let mut closest = f64::INFINITY;
    let mut fraction = 0.0;
    for w in path.windows(2) {
        let (a, b) = (w[0], w[1]);
        let d = point_seg_dist(p, a, b);
        let len = seg_len(a, b);
        let sf = ref_locate_seg(a, b, p)?;
        if d < closest {
            closest = d;
            fraction = (cum + sf * len) / total;
        }
        cum += len;
    }
    Some(fraction)
}

#[test]
fn line_interpolate_point_matches_reference() {
    let fractions = [
        -1.0,
        0.0,
        0.1,
        0.25,
        0.5,
        0.75,
        0.9,
        1.0,
        2.0,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NAN,
    ];
    let seg = ((-1.0, 0.0), (3.0, 2.0));
    let h_seg = Line::new(Coord::new(seg.0 .0, seg.0 .1), Coord::new(seg.1 .0, seg.1 .1));
    let paths = [
        vec![(-1.0, 0.0), (0.0, 0.0), (0.0, 1.0)],
        vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)],
        vec![(1.0, 1.0), (1.0, 1.0)],
        vec![],
    ];
    for &f in &fractions {
        let h = h_seg.line_interpolate_point(f);
        let r = ref_interpolate_seg(seg.0, seg.1, f);
        match (h, r) {
            (Some(hp), Some(rp)) => assert_points_eq(hp, rp, &format!("line interpolate {f}")),
            (None, None) => {}
            (h, r) => panic!("line interpolate {f}: {h:?} vs ref {r:?}"),
        }
        for path in &paths {
            let h = h_line_string(path).line_interpolate_point(f);
            // A two-point zero-length path yields a finite point; an empty path
            // yields None.
            let r = if path.len() < 2 {
                None
            } else {
                ref_interpolate(path, f)
            };
            match (h, r) {
                (Some(hp), Some(rp)) => {
                    assert_points_eq(hp, rp, &format!("ls interpolate {path:?} at {f}"))
                }
                (None, None) => {}
                (h, r) => panic!("ls interpolate {path:?} at {f}: {h:?} vs ref {r:?}"),
            }
        }
    }
}

#[test]
fn line_locate_point_matches_reference() {
    let seg = ((-1.0, 0.0), (3.0, 2.0));
    let h_seg = Line::new(Coord::new(seg.0 .0, seg.0 .1), Coord::new(seg.1 .0, seg.1 .1));
    let zero = ((1.0, 1.0), (1.0, 1.0));
    let h_zero = Line::new(Coord::new(zero.0 .0, zero.0 .1), Coord::new(zero.1 .0, zero.1 .1));
    let path = vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)];
    let h_ls = h_line_string(&path);

    let mut queries = query_points();
    queries.push((f64::NAN, 1.0));
    for (x, y) in queries {
        let q = (x, y);
        let p = Point::new(x, y);
        for (name, h, r) in [
            (
                "segment",
                h_seg.line_locate_point(&p),
                ref_locate_seg(seg.0, seg.1, q),
            ),
            (
                "zero",
                h_zero.line_locate_point(&p),
                ref_locate_seg(zero.0, zero.1, q),
            ),
            (
                "linestring",
                h_ls.line_locate_point(&p),
                ref_locate_path(&path, q),
            ),
        ] {
            match (h, r) {
                (Some(hf), Some(rf)) => assert!(
                    (hf - rf).abs() < EPS,
                    "{name} locate ({x}, {y}): {hf} vs ref {rf}"
                ),
                (None, None) => {}
                (h, r) => panic!("{name} locate ({x}, {y}): {h:?} vs ref {r:?}"),
            }
        }
    }

    // An empty (zero-length) line string locates everything at fraction 0.
    assert_eq!(
        h_line_string(&[]).line_locate_point(&Point::new(1.0, 1.0)),
        ref_locate_path(&[], (1.0, 1.0)),
    );
    assert_eq!(
        h_line_string(&[]).line_locate_point(&Point::new(1.0, 1.0)),
        Some(0.0),
    );
}

#[test]
fn locate_interpolate_round_trip_is_consistent() {
    let path = vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)];
    let h_ls = h_line_string(&path);
    for i in 0..=20 {
        let f = i as f64 / 20.0;
        let hp = h_ls.line_interpolate_point(f).unwrap();
        let rp = ref_interpolate(&path, f).unwrap();
        assert_points_eq(hp, rp, &format!("round-trip interpolate {f}"));
        // Locating the interpolated point must recover the original fraction.
        let hf = h_ls.line_locate_point(&hp).unwrap();
        let rf = ref_locate_path(&path, (hp.x(), hp.y())).unwrap();
        assert!(
            (hf - rf).abs() < EPS && (hf - f).abs() < 1e-12,
            "round-trip locate {f}: {hf} vs ref {rf}"
        );
    }
}

// ---------------------------------------------------------------------------
// Independent affine reference oracle: a 3×3 (six-parameter) matrix built and
// applied with our own algebra, mirroring the documented constructors.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
struct RefAffine {
    a: f64,
    b: f64,
    xoff: f64,
    d: f64,
    e: f64,
    yoff: f64,
}

impl RefAffine {
    fn identity() -> Self {
        Self { a: 1.0, b: 0.0, xoff: 0.0, d: 0.0, e: 1.0, yoff: 0.0 }
    }
    fn translate(x: f64, y: f64) -> Self {
        Self { a: 1.0, b: 0.0, xoff: x, d: 0.0, e: 1.0, yoff: y }
    }
    fn scale(xf: f64, yf: f64, origin: (f64, f64)) -> Self {
        let (x0, y0) = origin;
        Self {
            a: xf,
            b: 0.0,
            xoff: x0 - x0 * xf,
            d: 0.0,
            e: yf,
            yoff: y0 - y0 * yf,
        }
    }
    fn rotate(degrees: f64, origin: (f64, f64)) -> Self {
        let (s, c) = degrees.to_radians().sin_cos();
        let (x0, y0) = origin;
        Self {
            a: c,
            b: -s,
            xoff: x0 - x0 * c + y0 * s,
            d: s,
            e: c,
            yoff: y0 - x0 * s - y0 * c,
        }
    }
    fn skew(xs: f64, ys: f64, origin: (f64, f64)) -> Self {
        let (x0, y0) = origin;
        let mut tanx = xs.to_radians().tan();
        let mut tany = ys.to_radians().tan();
        if tanx.abs() < 2.5e-16 {
            tanx = 0.0;
        }
        if tany.abs() < 2.5e-16 {
            tany = 0.0;
        }
        Self {
            a: 1.0,
            b: tanx,
            xoff: -y0 * tanx,
            d: tany,
            e: 1.0,
            yoff: -x0 * tany,
        }
    }
    /// `self` then `other`: the matrix product `other · self`.
    fn then(self, other: RefAffine) -> RefAffine {
        RefAffine {
            a: other.a * self.a + other.b * self.d,
            b: other.a * self.b + other.b * self.e,
            xoff: other.a * self.xoff + other.b * self.yoff + other.xoff,
            d: other.d * self.a + other.e * self.d,
            e: other.d * self.b + other.e * self.e,
            yoff: other.d * self.xoff + other.e * self.yoff + other.yoff,
        }
    }
    fn apply(&self, p: (f64, f64)) -> (f64, f64) {
        (
            self.a * p.0 + self.b * p.1 + self.xoff,
            self.d * p.0 + self.e * p.1 + self.yoff,
        )
    }
    fn params(&self) -> [f64; 6] {
        [self.a, self.b, self.xoff, self.d, self.e, self.yoff]
    }
}

fn assert_line_strings_eq_ref(h: &LineString, r: &[(f64, f64)], ctx: &str) {
    assert_eq!(h.0.len(), r.len(), "{ctx}: coord count");
    for (hc, &(rx, ry)) in h.0.iter().zip(r.iter()) {
        assert!(
            (hc.x - rx).abs() < EPS && (hc.y - ry).abs() < EPS,
            "{ctx}: ({}, {}) vs ref ({rx}, {ry})",
            hc.x,
            hc.y,
        );
    }
}

fn assert_polygons_eq_ref(
    h: &Polygon,
    ext: &[(f64, f64)],
    holes: &[Vec<(f64, f64)>],
    ctx: &str,
) {
    assert_line_strings_eq_ref(h.exterior(), ext, ctx);
    assert_eq!(h.interiors().len(), holes.len(), "{ctx}: holes");
    for (hr, gr) in h.interiors().iter().zip(holes) {
        assert_line_strings_eq_ref(hr, gr, ctx);
    }
}

/// Hausdorff stores rings closed; produce the closed form a path takes once
/// inside a `Polygon`, then map each coord through `t`.
fn apply_ring(t: &RefAffine, ring: &[(f64, f64)]) -> Vec<(f64, f64)> {
    closed(ring).iter().map(|&c| t.apply(c)).collect()
}

fn apply_path(t: &RefAffine, path: &[(f64, f64)]) -> Vec<(f64, f64)> {
    path.iter().map(|&c| t.apply(c)).collect()
}

// The fixtures used across the affine tests, kept identical to the originals.
fn poly_fixture() -> (Polygon, Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>) {
    let exterior = vec![(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
    let hole = vec![(3.0, 3.0), (5.0, 3.0), (5.0, 5.0), (3.0, 5.0)];
    (
        h_polygon(&exterior, &[hole.clone()]),
        exterior,
        vec![hole],
    )
}

fn ls_fixture() -> (LineString, Vec<(f64, f64)>) {
    let path = vec![(0.0, 0.0), (5.0, 1.0), (7.0, 6.0), (-2.0, 4.0)];
    (h_line_string(&path), path)
}

#[test]
fn affine_transform_application_matches_reference() {
    let transforms: Vec<(AffineTransform, RefAffine)> = vec![
        (AffineTransform::identity(), RefAffine::identity()),
        (
            AffineTransform::translate(3.5, -2.25),
            RefAffine::translate(3.5, -2.25),
        ),
        (
            AffineTransform::scale(2.0, -0.5, Coord::new(1.0, 2.0)),
            RefAffine::scale(2.0, -0.5, (1.0, 2.0)),
        ),
        (
            AffineTransform::rotate(37.0, Coord::new(-1.0, 4.0)),
            RefAffine::rotate(37.0, (-1.0, 4.0)),
        ),
        (
            AffineTransform::skew(15.0, -30.0, Coord::new(2.0, 2.0)),
            RefAffine::skew(15.0, -30.0, (2.0, 2.0)),
        ),
        (
            AffineTransform::translate(1.0, 2.0)
                .rotated(45.0, Coord::new(0.0, 0.0))
                .scaled(3.0, 0.25, Coord::new(5.0, 5.0))
                .skewed(10.0, 20.0, Coord::new(-2.0, 3.0)),
            RefAffine::translate(1.0, 2.0)
                .then(RefAffine::rotate(45.0, (0.0, 0.0)))
                .then(RefAffine::scale(3.0, 0.25, (5.0, 5.0)))
                .then(RefAffine::skew(10.0, 20.0, (-2.0, 3.0))),
        ),
    ];
    let (h_poly, ext, holes) = poly_fixture();
    let (h_ls, ls_path) = ls_fixture();
    for (i, (ht, rt)) in transforms.iter().enumerate() {
        assert_eq!(
            [ht.a(), ht.b(), ht.xoff(), ht.d(), ht.e(), ht.yoff()],
            rt.params(),
            "transform {i}: matrix"
        );
        let ref_ext = apply_ring(rt, &ext);
        let ref_holes: Vec<Vec<(f64, f64)>> = holes.iter().map(|h| apply_ring(rt, h)).collect();
        assert_polygons_eq_ref(
            &h_poly.affine_transform(ht),
            &ref_ext,
            &ref_holes,
            &format!("transform {i} on polygon"),
        );
        assert_line_strings_eq_ref(
            &h_ls.affine_transform(ht),
            &apply_path(rt, &ls_path),
            &format!("transform {i} on linestring"),
        );
        let mut h_mut = h_ls.clone();
        h_mut.affine_transform_mut(ht);
        assert_eq!(h_mut, h_ls.affine_transform(ht), "transform {i} mut");
    }
}

/// Centroid of a hausdorff geometry, taken from its own verified `Centroid`
/// kernel — the anchor for the centroid-pivot convenience transforms.
fn centroid_of(g: &impl Centroid) -> (f64, f64) {
    let c = g.centroid().expect("non-empty fixture has a centroid");
    (c.x(), c.y())
}

/// Bounding-box centre, computed independently from the geometry's coords.
fn bbox_center(coords: &[(f64, f64)]) -> (f64, f64) {
    let (mut minx, mut miny) = (f64::INFINITY, f64::INFINITY);
    let (mut maxx, mut maxy) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(x, y) in coords {
        minx = minx.min(x);
        miny = miny.min(y);
        maxx = maxx.max(x);
        maxy = maxy.max(y);
    }
    ((minx + maxx) / 2.0, (miny + maxy) / 2.0)
}

#[test]
fn translate_rotate_scale_skew_match_reference() {
    let (h_poly, ext, holes) = poly_fixture();
    let (h_ls, ls_path) = ls_fixture();

    // The polygon's bounding box is over its closed exterior ring; the
    // linestring's over its own coords. Centroids come from hausdorff's kernel.
    let poly_bbc = bbox_center(&closed(&ext));
    let poly_centroid = centroid_of(&h_poly);
    let ls_centroid = centroid_of(&h_ls);

    let apply_poly = |t: &RefAffine| -> (Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>) {
        (
            apply_ring(t, &ext),
            holes.iter().map(|h| apply_ring(t, h)).collect(),
        )
    };

    // translate
    {
        let t = RefAffine::translate(17.0, -3.0);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.translate(17.0, -3.0), &re, &rh, "translate polygon");
    }
    {
        let t = RefAffine::translate(-0.5, 9.0);
        assert_line_strings_eq_ref(
            &h_ls.translate(-0.5, 9.0),
            &apply_path(&t, &ls_path),
            "translate linestring",
        );
    }

    // rotate
    {
        let t = RefAffine::rotate(-45.0, (10.0, 0.0));
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(
            &h_poly.rotate_around_point(-45.0, Point::new(10.0, 0.0)),
            &re,
            &rh,
            "rotate around point",
        );
    }
    {
        let t = RefAffine::rotate(33.0, poly_centroid);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(
            &h_poly.rotate_around_centroid(33.0),
            &re,
            &rh,
            "rotate around centroid",
        );
    }
    {
        let t = RefAffine::rotate(-120.0, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(
            &h_poly.rotate_around_center(-120.0),
            &re,
            &rh,
            "rotate around center",
        );
    }
    {
        let t = RefAffine::rotate(75.0, ls_centroid);
        assert_line_strings_eq_ref(
            &h_ls.rotate_around_centroid(75.0),
            &apply_path(&t, &ls_path),
            "rotate linestring around centroid",
        );
    }

    // scale (about bounding-box centre)
    {
        let t = RefAffine::scale(2.5, 2.5, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.scale(2.5), &re, &rh, "scale");
    }
    {
        let t = RefAffine::scale(2.0, -4.0, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.scale_xy(2.0, -4.0), &re, &rh, "scale_xy");
    }
    {
        let t = RefAffine::scale(2.0, 4.0, (100.0, 100.0));
        assert_line_strings_eq_ref(
            &h_ls.scale_around_point(2.0, 4.0, Coord::new(100.0, 100.0)),
            &apply_path(&t, &ls_path),
            "scale_around_point",
        );
    }

    // skew (about bounding-box centre)
    {
        let t = RefAffine::skew(30.0, 30.0, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.skew(30.0), &re, &rh, "skew");
    }
    {
        let t = RefAffine::skew(30.0, 12.0, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.skew_xy(30.0, 12.0), &re, &rh, "skew_xy");
    }
    {
        let t = RefAffine::skew(45.0, 10.0, (2.0, 2.0));
        assert_line_strings_eq_ref(
            &h_ls.skew_around_point(45.0, 10.0, Coord::new(2.0, 2.0)),
            &apply_path(&t, &ls_path),
            "skew_around_point",
        );
    }

    // mut variants must agree with their immutable counterparts.
    let mut h_mut = h_poly.clone();
    h_mut.rotate_around_centroid_mut(33.0);
    assert_eq!(h_mut, h_poly.rotate_around_centroid(33.0), "rotate mut");
    let mut h_mut = h_ls.clone();
    h_mut.scale_mut(2.5);
    assert_eq!(h_mut, h_ls.scale(2.5), "scale mut");
    let mut h_mut = h_ls.clone();
    h_mut.skew_xy_mut(30.0, 12.0);
    assert_eq!(h_mut, h_ls.skew_xy(30.0, 12.0), "skew mut");
    let mut h_mut = h_ls.clone();
    h_mut.translate_mut(1.0, 2.0);
    assert_eq!(h_mut, h_ls.translate(1.0, 2.0), "translate mut");
}

// ---------------------------------------------------------------------------
// Independent extremes reference: a strict-comparison scan over the exterior
// coordinate walk (interior rings excluded, first coord wins ties).
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
struct RefExtreme {
    index: usize,
    coord: (f64, f64),
}

#[derive(Clone, Copy, Debug)]
struct RefExtremes {
    x_min: RefExtreme,
    y_min: RefExtreme,
    x_max: RefExtreme,
    y_max: RefExtreme,
}

/// Compute extremes over a flat coordinate walk, exactly mirroring the kernel:
/// the *first* coordinate attaining an extreme value wins (strict `<` / `>`).
fn ref_extremes(coords: &[(f64, f64)]) -> Option<RefExtremes> {
    let &first = coords.first()?;
    let mut out = RefExtremes {
        x_min: RefExtreme { index: 0, coord: first },
        y_min: RefExtreme { index: 0, coord: first },
        x_max: RefExtreme { index: 0, coord: first },
        y_max: RefExtreme { index: 0, coord: first },
    };
    for (index, &coord) in coords.iter().enumerate().skip(1) {
        if coord.0 < out.x_min.coord.0 {
            out.x_min = RefExtreme { index, coord };
        }
        if coord.1 < out.y_min.coord.1 {
            out.y_min = RefExtreme { index, coord };
        }
        if coord.0 > out.x_max.coord.0 {
            out.x_max = RefExtreme { index, coord };
        }
        if coord.1 > out.y_max.coord.1 {
            out.y_max = RefExtreme { index, coord };
        }
    }
    Some(out)
}

#[test]
fn extremes_match_reference() {
    let diamond = vec![(1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0), (1.0, 0.0)];
    let irregular = vec![
        (1.0, -3.1),
        (3.7, 2.7),
        (0.9, 7.6),
        (-4.8, 6.7),
        (-7.5, 0.9),
    ];

    // The exterior coordinate walk a geometry exposes to `Extremes`:
    //  - a polygon: its closed exterior ring;
    //  - a linestring: its coords as given;
    //  - a multipoint: its points in order.
    let h_cases: Vec<Geometry> = vec![
        Geometry::Polygon(h_polygon(&diamond, &[])),
        Geometry::Polygon(h_polygon(&irregular, &[])),
        Geometry::LineString(h_line_string(&irregular)),
        Geometry::MultiPoint(MultiPoint::from(irregular.clone())),
    ];
    let r_walks: Vec<Vec<(f64, f64)>> = vec![
        closed(&diamond),
        closed(&irregular),
        irregular.clone(),
        irregular.clone(),
    ];

    for (h, walk) in h_cases.iter().zip(&r_walks) {
        let ho = h.extremes().unwrap();
        let ro = ref_extremes(walk).unwrap();
        for (name, he, re) in [
            ("x_min", ho.x_min, ro.x_min),
            ("y_min", ho.y_min, ro.y_min),
            ("x_max", ho.x_max, ro.x_max),
            ("y_max", ho.y_max, ro.y_max),
        ] {
            assert_eq!(he.index, re.index, "{name} index on {h:?}");
            assert!(
                (he.coord.x - re.coord.0).abs() < EPS && (he.coord.y - re.coord.1).abs() < EPS,
                "{name} coord on {h:?}: {:?} vs ref {:?}",
                he.coord,
                re.coord
            );
        }
    }

    assert!(MultiPoint::new(vec![]).extremes().is_none());
}

#[test]
fn extremes_ignore_interior_rings() {
    let exterior = vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
    let hole = vec![(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)];
    let h = h_polygon(&exterior, &[hole.clone()]);
    // The reference walks only the (closed) exterior ring, never the hole.
    let ro = ref_extremes(&closed(&exterior)).unwrap();
    let ho = h.extremes().unwrap();
    assert_eq!(ho.x_max.index, ro.x_max.index);
    assert_eq!(ho.x_max.coord.x, ro.x_max.coord.0);
    assert!(ho.x_max.index < h.exterior().coords_count());
}
