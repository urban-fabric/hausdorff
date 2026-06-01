//! Native concave-hull (chi-shape) invariants. We don't match geo's
//! `concaveman` output vertex-for-vertex (different heuristic), so we validate
//! the geometric guarantees a concave hull must satisfy.

use hausdorff::{Area, ConcaveHull, ConvexHull, Coord, MultiPoint, Point};

fn mp(pts: &[(f64, f64)]) -> MultiPoint {
    MultiPoint::new(pts.iter().map(|&(x, y)| Point::new(x, y)).collect())
}

fn orient(a: Coord, b: Coord, c: Coord) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// Two segments share a strictly-interior point (no shared endpoints).
fn proper_cross(a: Coord, b: Coord, c: Coord, d: Coord) -> bool {
    let (d1, d2, d3, d4) = (orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b));
    d1 != 0.0
        && d2 != 0.0
        && d3 != 0.0
        && d4 != 0.0
        && (d1 > 0.0) != (d2 > 0.0)
        && (d3 > 0.0) != (d4 > 0.0)
}

/// A polygon ring has no self-intersections among non-adjacent edges.
fn is_simple(ring: &[Coord]) -> bool {
    let v: Vec<Coord> = ring[..ring.len() - 1].to_vec(); // drop closing dup
    let n = v.len();
    for i in 0..n {
        let (a, b) = (v[i], v[(i + 1) % n]);
        for j in (i + 1)..n {
            // skip the two edges adjacent to edge i (they share a vertex)
            if j == i || (j + 1) % n == i || (i + 1) % n == j {
                continue;
            }
            let (c, d) = (v[j], v[(j + 1) % n]);
            if proper_cross(a, b, c, d) {
                return false;
            }
        }
    }
    true
}

fn point_on_seg(p: Coord, a: Coord, b: Coord) -> bool {
    orient(a, b, p) == 0.0
        && p.x <= a.x.max(b.x)
        && p.x >= a.x.min(b.x)
        && p.y <= a.y.max(b.y)
        && p.y >= a.y.min(b.y)
}

/// A point is covered if it is strictly inside or lies on a ring edge.
fn covered(ring: &[Coord], p: Coord) -> bool {
    // ray cast for strict interior
    let n = ring.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (pi, pj) = (ring[i], ring[j]);
        if (pi.y > p.y) != (pj.y > p.y)
            && p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x
        {
            inside = !inside;
        }
        j = i;
    }
    if inside {
        return true;
    }
    ring.windows(2).any(|w| point_on_seg(p, w[0], w[1]))
}

#[test]
fn covers_all_points_and_is_simple() {
    let pts = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (5.0, 1.0),
        (3.0, 4.0),
        (7.0, 6.0),
        (2.0, 8.0),
    ];
    let hull = mp(&pts).concave_hull(2.0);
    let ring = &hull.exterior().0;

    assert_eq!(ring.first(), ring.last(), "ring must be closed");
    assert!(is_simple(ring), "concave hull must be a simple polygon");
    for &(x, y) in &pts {
        assert!(covered(ring, Coord::new(x, y)), "point ({x},{y}) not covered");
    }
}

#[test]
fn never_exceeds_convex_hull_area() {
    let pts = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (5.0, 1.0),
        (4.0, 5.0),
    ];
    let convex = mp(&pts).convex_hull().unsigned_area();
    let concave = mp(&pts).concave_hull(2.0).unsigned_area();
    assert!(concave > 0.0);
    assert!(
        concave <= convex + 1e-9,
        "concave area {concave} exceeds convex area {convex}"
    );
}

#[test]
fn large_concavity_degenerates_to_convex_hull() {
    let pts = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (5.0, 1.0)];
    let convex = mp(&pts).convex_hull().unsigned_area();
    let concave = mp(&pts).concave_hull(1.0e6).unsigned_area();
    assert!((concave - convex).abs() < 1e-9, "huge concavity should equal convex hull");
}

#[test]
fn small_concavity_digs_in() {
    // square corners + a point just inside the bottom edge: a low concavity
    // should pull the long bottom edge in to that point, shrinking the area.
    let pts = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (5.0, 1.0)];
    let convex = mp(&pts).convex_hull().unsigned_area();
    let concave = mp(&pts).concave_hull(2.0).unsigned_area();
    assert!(
        concave < convex - 1.0,
        "small concavity should dig in: concave {concave} vs convex {convex}"
    );
}
