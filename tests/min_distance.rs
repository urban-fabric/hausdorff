//! Conformance for the native minimum-distance kernel (the `ST_Distance` /
//! `ST_DWithin` substrate) against geo's Euclidean geometry-to-geometry distance.

use geo::{Distance as _, Euclidean};
use hausdorff::space::r2::{MinDistanceWorkspace, euclidean_distance, euclidean_within};
use hausdorff::{Coord, Geometry, Line, LineString, Point, Polygon};

mod common;
use common::to_geo_geometry;

fn pt(x: f64, y: f64) -> Geometry {
    Geometry::Point(Point::new(x, y))
}
fn line(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    Geometry::Line(Line::new(Coord::new(x0, y0), Coord::new(x1, y1)))
}
fn ls(coords: &[(f64, f64)]) -> Geometry {
    Geometry::LineString(LineString::from(coords.to_vec()))
}
fn poly(coords: &[(f64, f64)]) -> Geometry {
    Geometry::Polygon(Polygon::new(LineString::from(coords.to_vec()), vec![]))
}
fn sq(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    poly(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
}
fn holed(shell: &[(f64, f64)], hole: &[(f64, f64)]) -> Geometry {
    Geometry::Polygon(Polygon::new(
        LineString::from(shell.to_vec()),
        vec![LineString::from(hole.to_vec())],
    ))
}

/// Compare the native min-distance to geo, and check symmetry + DWithin.
fn check(a: &Geometry, b: &Geometry) {
    let ours = euclidean_distance(a, b);
    let (ga, gb) = (to_geo_geometry(a), to_geo_geometry(b));
    let theirs = Euclidean.distance(&ga, &gb);
    assert!(
        (ours - theirs).abs() <= 1e-9 * theirs.max(1.0),
        "distance native={ours} geo={theirs}\n  A={a:?}\n  B={b:?}"
    );
    // symmetric
    assert!(
        (euclidean_distance(b, a) - ours).abs() < 1e-9,
        "asymmetric distance"
    );
    // DWithin is consistent with the computed distance. (`ours` is √(min_sq);
    // squaring it back inside `within` can round a hair under the exact min, so
    // test just above the distance.)
    assert!(
        euclidean_within(a, b, ours + 1e-9 * ours.max(1.0)),
        "within its own distance"
    );
    assert!(
        !euclidean_within(a, b, -1.0),
        "never within a negative radius"
    );
    if ours > 1e-6 {
        assert!(
            !euclidean_within(a, b, ours - 1e-6),
            "not within less than the gap"
        );
    }
}

#[test]
fn point_point() {
    check(&pt(0.0, 0.0), &pt(3.0, 4.0)); // 5
    check(&pt(1.0, 1.0), &pt(1.0, 1.0)); // 0
}

#[test]
fn point_to_segment_and_linestring() {
    check(&pt(0.0, 0.0), &line(2.0, -1.0, 2.0, 1.0)); // perpendicular → 2
    check(&pt(5.0, 0.0), &line(0.0, 0.0, 3.0, 0.0)); // beyond endpoint → 2
    check(&pt(1.0, 3.0), &ls(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]));
    check(&pt(2.0, 1.0), &ls(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)])); // on the LS → 0
}

#[test]
fn line_line() {
    check(&line(0.0, 0.0, 1.0, 1.0), &line(2.0, 0.0, 3.0, 1.0)); // disjoint
    check(&line(0.0, 0.0, 2.0, 2.0), &line(0.0, 2.0, 2.0, 0.0)); // crossing → 0
    check(&line(0.0, 0.0, 4.0, 0.0), &line(1.0, 3.0, 3.0, 3.0)); // parallel → 3
}

#[test]
fn point_polygon_containment() {
    let s = sq(0.0, 0.0, 4.0, 4.0);
    check(&pt(2.0, 2.0), &s); // inside → 0
    check(&pt(6.0, 2.0), &s); // outside → 2 from the right edge
    check(&pt(-3.0, -4.0), &s); // outside corner
}

#[test]
fn polygon_polygon() {
    check(&sq(0.0, 0.0, 2.0, 2.0), &sq(5.0, 0.0, 7.0, 2.0)); // disjoint, gap 3
    check(&sq(0.0, 0.0, 2.0, 2.0), &sq(2.0, 0.0, 4.0, 2.0)); // shared edge → 0
    check(&sq(0.0, 0.0, 10.0, 10.0), &sq(3.0, 3.0, 7.0, 7.0)); // nested → 0
    check(&sq(0.0, 0.0, 2.0, 2.0), &sq(3.0, 3.0, 5.0, 5.0)); // diagonal gap
}

#[test]
fn polygon_with_hole() {
    let donut = holed(
        &[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ],
        &[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)],
    );
    check(&donut, &pt(5.0, 5.0)); // in the hole → distance to hole ring (2)
    check(&donut, &pt(1.0, 1.0)); // in the solid ring → 0
    check(&donut, &sq(4.0, 4.0, 6.0, 6.0)); // square inside the hole
}

// ---- randomized stress against geo -----------------------------------------

struct Lcg(u64);
impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
    fn coord(&mut self, n: i64) -> f64 {
        f64::from(self.next_u32() % (n as u32)) - f64::from(n as u32) / 2.0
    }
}

#[test]
fn large_polygons_use_index_and_match_geo() {
    // Rings with >64 segments each push the pair count past INDEX_THRESHOLD, so
    // the spatial-index path runs — it must agree with geo (and the brute path).
    let ring = |cx: f64, cy: f64, r: f64, n: usize| -> Vec<(f64, f64)> {
        (0..=n)
            .map(|i| {
                let t = i as f64 / n as f64 * std::f64::consts::TAU;
                (cx + r * t.cos(), cy + r * t.sin())
            })
            .collect()
    };
    let mk = |c: Vec<(f64, f64)>| Geometry::Polygon(Polygon::new(LineString::from(c), vec![]));

    check(
        &mk(ring(0.0, 0.0, 10.0, 200)),
        &mk(ring(40.0, 5.0, 8.0, 200)),
    ); // disjoint
    check(
        &mk(ring(0.0, 0.0, 20.0, 200)),
        &mk(ring(2.0, 1.0, 5.0, 150)),
    ); // nested → 0
    check(
        &mk(ring(0.0, 0.0, 10.0, 200)),
        &mk(ring(15.0, 0.0, 10.0, 180)),
    ); // overlap → 0
    check(
        &mk(ring(0.0, 0.0, 10.0, 128)),
        &mk(ring(0.0, 30.0, 10.0, 128)),
    ); // stacked, gap 10
}

#[test]
fn index_path_scattered_points_match_geo() {
    // Two overlapping clouds of many points: the cost model picks the spatial
    // index (overlapping boxes + large), which must still agree with geo.
    let mut rng = Lcg(0x1357_9BDF_2468_ACE0);
    let cloud = |rng: &mut Lcg, n: usize| -> Geometry {
        Geometry::MultiPoint(hausdorff::MultiPoint::from(
            (0..n)
                .map(|_| (rng.coord(100), rng.coord(100)))
                .collect::<Vec<_>>(),
        ))
    };
    for _ in 0..20 {
        let a = cloud(&mut rng, 150);
        let b = cloud(&mut rng, 150);
        check(&a, &b);
    }
}

/// A regular `n`-gon (always convex) — closed ring, arbitrary rotation. Drives
/// the convex GJK fast path at the vertex counts the brute scan struggles with.
fn ngon(cx: f64, cy: f64, r: f64, n: usize, rot: f64) -> Geometry {
    let coords: Vec<(f64, f64)> = (0..=n)
        .map(|i| {
            let t = rot + i as f64 / n as f64 * std::f64::consts::TAU;
            (cx + r * t.cos(), cy + r * t.sin())
        })
        .collect();
    poly(&coords)
}

#[test]
fn convex_fast_path_matches_geo() {
    // Random convex polygons (regular n-gons, varied size/rotation/position):
    // disjoint pairs take the O(n+m) GJK path, touching/overlapping pairs the
    // intersects short-circuit — both must agree with geo, and with the squares
    // and triangles already covered by `polygon_polygon` / the stress test.
    let mut rng = Lcg(0x00C0_FFEE_1234_5678);
    for _ in 0..400 {
        let n1 = 3 + (rng.next_u32() % 32) as usize;
        let n2 = 3 + (rng.next_u32() % 32) as usize;
        let r1 = 1.0 + rng.coord(8).abs();
        let r2 = 1.0 + rng.coord(8).abs();
        let a = ngon(rng.coord(40), rng.coord(40), r1, n1, rng.coord(7));
        let b = ngon(rng.coord(40), rng.coord(40), r2, n2, rng.coord(7));
        check(&a, &b);
    }
    // Large convex rings well separated: pure GJK, no intersects involvement.
    check(
        &ngon(0.0, 0.0, 5.0, 256, 0.1),
        &ngon(30.0, 4.0, 4.0, 512, 0.7),
    );
    check(
        &ngon(0.0, 0.0, 5.0, 512, 0.0),
        &ngon(0.0, 25.0, 5.0, 200, 1.3),
    ); // stacked
}

#[test]
fn workspace_matches_free_functions() {
    // A reused workspace must return exactly what the one-shot free functions do,
    // across a heterogeneous batch (small, convex, point-cloud) of queries.
    let mut ws = MinDistanceWorkspace::new();
    let cases: Vec<(Geometry, Geometry)> = vec![
        (pt(0.0, 0.0), pt(3.0, 4.0)),
        (sq(0.0, 0.0, 2.0, 2.0), sq(5.0, 0.0, 7.0, 2.0)),
        (ngon(0.0, 0.0, 5.0, 64, 0.0), ngon(20.0, 3.0, 4.0, 64, 0.0)),
        (line(0.0, 0.0, 4.0, 0.0), pt(2.0, 3.0)),
        (sq(0.0, 0.0, 10.0, 10.0), sq(3.0, 3.0, 7.0, 7.0)), // nested → 0
    ];
    // Two passes over the same workspace to confirm buffer reuse is stateless.
    for _ in 0..2 {
        for (a, b) in &cases {
            assert_eq!(
                ws.distance(a, b),
                euclidean_distance(a, b),
                "distance reuse mismatch"
            );
            let d = euclidean_distance(a, b);
            assert_eq!(ws.within(a, b, d + 0.5), euclidean_within(a, b, d + 0.5));
            assert_eq!(
                ws.within(a, b, (d - 0.5).max(0.0)),
                euclidean_within(a, b, (d - 0.5).max(0.0))
            );
        }
    }
}

#[test]
fn stress_random_geometries_match_geo() {
    let mut rng = Lcg(0x51ED_2317_BEEF_0042);
    for _ in 0..600 {
        let mk = |rng: &mut Lcg, kind: u32| -> Geometry {
            match kind % 4 {
                0 => pt(rng.coord(20), rng.coord(20)),
                1 => line(rng.coord(20), rng.coord(20), rng.coord(20), rng.coord(20)),
                2 => {
                    let (x, y) = (rng.coord(16), rng.coord(16));
                    let (w, h) = (2.0 + rng.coord(8).abs(), 2.0 + rng.coord(8).abs());
                    sq(x, y, x + w, y + h)
                }
                _ => poly(&[
                    (rng.coord(20), rng.coord(20)),
                    (rng.coord(20), rng.coord(20)),
                    (rng.coord(20), rng.coord(20)),
                ]),
            }
        };
        let ka = rng.next_u32();
        let a = mk(&mut rng, ka);
        let kb = rng.next_u32();
        let b = mk(&mut rng, kb);
        let ours = euclidean_distance(&a, &b);
        let theirs = Euclidean.distance(&to_geo_geometry(&a), &to_geo_geometry(&b));
        assert!(
            (ours - theirs).abs() <= 1e-7 * theirs.max(1.0),
            "native={ours} geo={theirs}\n  A={a:?}\n  B={b:?}"
        );
    }
}
