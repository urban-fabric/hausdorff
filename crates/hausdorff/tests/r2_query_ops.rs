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

#[path = "r2_query_ops/closest.rs"]
mod closest;

#[path = "r2_query_ops/interpolate.rs"]
mod interpolate;

#[path = "r2_query_ops/affine.rs"]
mod affine;

#[path = "r2_query_ops/extremes.rs"]
mod extremes;
