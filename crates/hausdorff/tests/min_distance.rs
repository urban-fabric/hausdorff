//! Conformance for the native minimum-distance kernel (the `ST_Distance` /
//! `ST_DWithin` substrate) against an independent brute-force segment oracle,
//! plus self-contained symmetry / DWithin invariants. No third-party crate.

use hausdorff::space::r2::Intersects as _;
use hausdorff::space::r2::{MinDistanceWorkspace, euclidean_distance, euclidean_within};
use hausdorff::{Coord, Geometry, Line, LineString, Point, Polygon};

/// A boundary segment as its two endpoints.
type Segment = ((f64, f64), (f64, f64));

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

// ---- independent brute-force min-distance oracle (dependency-free) ----

/// Distance from point `p` to segment `ab`.
fn point_seg(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (vx, vy) = (b.0 - a.0, b.1 - a.1);
    let (wx, wy) = (p.0 - a.0, p.1 - a.1);
    let len2 = vx * vx + vy * vy;
    let t = if len2 <= 0.0 {
        0.0
    } else {
        ((wx * vx + wy * vy) / len2).clamp(0.0, 1.0)
    };
    let (cx, cy) = (a.0 + t * vx, a.1 + t * vy);
    ((p.0 - cx).powi(2) + (p.1 - cy).powi(2)).sqrt()
}

/// Distance between two segments known not to cross (endpoints-to-segment min).
fn seg_seg(p1: (f64, f64), p2: (f64, f64), q1: (f64, f64), q2: (f64, f64)) -> f64 {
    point_seg(p1, q1, q2)
        .min(point_seg(p2, q1, q2))
        .min(point_seg(q1, p1, p2))
        .min(point_seg(q2, p1, p2))
}

/// Decompose a geometry into its boundary segments and isolated points.
fn primitives(g: &Geometry) -> (Vec<Segment>, Vec<(f64, f64)>) {
    let mut segs = Vec::new();
    let mut pts = Vec::new();
    fn push_ring(segs: &mut Vec<Segment>, coords: &[Coord]) {
        for w in coords.windows(2) {
            segs.push(((w[0].x, w[0].y), (w[1].x, w[1].y)));
        }
    }
    match g {
        Geometry::Point(p) => pts.push((p.x(), p.y())),
        Geometry::MultiPoint(mp) => pts.extend(mp.iter().map(|p| (p.x(), p.y()))),
        Geometry::Line(l) => segs.push(((l.start.x, l.start.y), (l.end.x, l.end.y))),
        Geometry::LineString(l) => push_ring(&mut segs, &l.0),
        Geometry::Polygon(poly) => {
            push_ring(&mut segs, &poly.exterior().0);
            for h in poly.interiors() {
                push_ring(&mut segs, &h.0);
            }
        }
        _ => {}
    }
    (segs, pts)
}

/// Independent O(n·m) brute-force min distance: 0 if the geometries intersect
/// (via the separate Intersects kernel), else the min over all boundary
/// primitive pairs.
fn brute_distance(a: &Geometry, b: &Geometry) -> f64 {
    if a.intersects(b) {
        return 0.0;
    }
    let (sa, pa) = primitives(a);
    let (sb, pb) = primitives(b);
    let mut best = f64::INFINITY;
    for &(a1, a2) in &sa {
        for &(b1, b2) in &sb {
            best = best.min(seg_seg(a1, a2, b1, b2));
        }
        for &p in &pb {
            best = best.min(point_seg(p, a1, a2));
        }
    }
    for &p in &pa {
        for &(b1, b2) in &sb {
            best = best.min(point_seg(p, b1, b2));
        }
        for &q in &pb {
            best = best.min(((p.0 - q.0).powi(2) + (p.1 - q.1).powi(2)).sqrt());
        }
    }
    best
}

/// Compare the native min-distance to the brute oracle, plus symmetry + DWithin.
fn check(a: &Geometry, b: &Geometry) {
    let ours = euclidean_distance(a, b);
    let reference = brute_distance(a, b);
    assert!(
        (ours - reference).abs() <= 1e-9 * reference.max(1.0),
        "distance native={ours} brute={reference}\n  A={a:?}\n  B={b:?}"
    );
    assert!(
        (euclidean_distance(b, a) - ours).abs() < 1e-9,
        "asymmetric distance"
    );
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
    check(&pt(0.0, 0.0), &pt(3.0, 4.0));
    check(&pt(1.0, 1.0), &pt(1.0, 1.0));
}

#[test]
fn point_to_segment_and_linestring() {
    check(&pt(0.0, 0.0), &line(2.0, -1.0, 2.0, 1.0));
    check(&pt(5.0, 0.0), &line(0.0, 0.0, 3.0, 0.0));
    check(&pt(1.0, 3.0), &ls(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]));
    check(&pt(2.0, 1.0), &ls(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]));
}

#[test]
fn line_line() {
    check(&line(0.0, 0.0, 1.0, 1.0), &line(2.0, 0.0, 3.0, 1.0));
    check(&line(0.0, 0.0, 2.0, 2.0), &line(0.0, 2.0, 2.0, 0.0));
    check(&line(0.0, 0.0, 4.0, 0.0), &line(1.0, 3.0, 3.0, 3.0));
}

#[test]
fn point_polygon_containment() {
    let s = sq(0.0, 0.0, 4.0, 4.0);
    check(&pt(2.0, 2.0), &s);
    check(&pt(6.0, 2.0), &s);
    check(&pt(-3.0, -4.0), &s);
}

#[test]
fn polygon_polygon() {
    check(&sq(0.0, 0.0, 2.0, 2.0), &sq(5.0, 0.0, 7.0, 2.0));
    check(&sq(0.0, 0.0, 2.0, 2.0), &sq(2.0, 0.0, 4.0, 2.0));
    check(&sq(0.0, 0.0, 10.0, 10.0), &sq(3.0, 3.0, 7.0, 7.0));
    check(&sq(0.0, 0.0, 2.0, 2.0), &sq(3.0, 3.0, 5.0, 5.0));
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
    check(&donut, &pt(5.0, 5.0));
    check(&donut, &pt(1.0, 1.0));
    check(&donut, &sq(4.0, 4.0, 6.0, 6.0));
}

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
    );
    check(
        &mk(ring(0.0, 0.0, 20.0, 200)),
        &mk(ring(2.0, 1.0, 5.0, 150)),
    );
    check(
        &mk(ring(0.0, 0.0, 10.0, 200)),
        &mk(ring(15.0, 0.0, 10.0, 180)),
    );
    check(
        &mk(ring(0.0, 0.0, 10.0, 128)),
        &mk(ring(0.0, 30.0, 10.0, 128)),
    );
}

#[test]
fn index_path_scattered_points_match_geo() {
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
    check(
        &ngon(0.0, 0.0, 5.0, 256, 0.1),
        &ngon(30.0, 4.0, 4.0, 512, 0.7),
    );
    check(
        &ngon(0.0, 0.0, 5.0, 512, 0.0),
        &ngon(0.0, 25.0, 5.0, 200, 1.3),
    );
}

#[test]
fn workspace_matches_free_functions() {
    let mut ws = MinDistanceWorkspace::new();
    let cases: Vec<(Geometry, Geometry)> = vec![
        (pt(0.0, 0.0), pt(3.0, 4.0)),
        (sq(0.0, 0.0, 2.0, 2.0), sq(5.0, 0.0, 7.0, 2.0)),
        (ngon(0.0, 0.0, 5.0, 64, 0.0), ngon(20.0, 3.0, 4.0, 64, 0.0)),
        (line(0.0, 0.0, 4.0, 0.0), pt(2.0, 3.0)),
        (sq(0.0, 0.0, 10.0, 10.0), sq(3.0, 3.0, 7.0, 7.0)),
    ];
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
        let reference = brute_distance(&a, &b);
        assert!(
            (ours - reference).abs() <= 1e-7 * reference.max(1.0),
            "native={ours} brute={reference}\n  A={a:?}\n  B={b:?}"
        );
    }
}
