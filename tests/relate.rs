//! Conformance for the native DE-9IM relate against geo's `Relate`. Every OGC
//! predicate must agree with geo across a battery of geometry-type pairs.

use geo::Relate;
use hausdorff::space::r2::relate::{geometry_dimension, relate};
use hausdorff::{Coord, Geometry, Line, LineString, MultiPolygon, Point, Polygon};

mod common;
use common::to_geo_geometry;

fn poly(coords: &[(f64, f64)]) -> Geometry {
    Geometry::Polygon(Polygon::new(LineString::from(coords.to_vec()), vec![]))
}
fn sq(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    poly(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
}
fn pt(x: f64, y: f64) -> Geometry {
    Geometry::Point(Point::new(x, y))
}
fn line(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    Geometry::Line(Line::new(Coord::new(x0, y0), Coord::new(x1, y1)))
}
fn ls(coords: &[(f64, f64)]) -> Geometry {
    Geometry::LineString(LineString::from(coords.to_vec()))
}

/// Compare every predicate of the native relate to geo's for one ordered pair.
fn check(a: &Geometry, b: &Geometry) {
    let m = relate(a, b);
    let (da, db) = (geometry_dimension(a), geometry_dimension(b));
    let (ga, gb) = (to_geo_geometry(a), to_geo_geometry(b));
    let g = ga.relate(&gb);

    let cases = [
        ("disjoint", m.is_disjoint(), g.is_disjoint()),
        ("intersects", m.is_intersects(), g.is_intersects()),
        ("within", m.is_within(), g.is_within()),
        ("contains", m.is_contains(), g.is_contains()),
        ("covers", m.is_covers(), g.is_covers()),
        ("coveredby", m.is_coveredby(), g.is_coveredby()),
        ("equals", m.is_equal_topo(), g.is_equal_topo()),
        ("touches", m.is_touches(), g.is_touches()),
        ("crosses", m.is_crosses(da, db), g.is_crosses()),
        ("overlaps", m.is_overlaps(da, db), g.is_overlaps()),
    ];
    for (name, got, want) in cases {
        assert_eq!(
            got, want,
            "{name}: native={got} geo={want}\n  A={a:?}\n  B={b:?}\n  IM={m:?}"
        );
    }
}

/// Check both orderings.
fn check2(a: &Geometry, b: &Geometry) {
    check(a, b);
    check(b, a);
}

#[test]
fn polygon_polygon() {
    let a = sq(0.0, 0.0, 4.0, 4.0);
    check2(&a, &sq(2.0, 2.0, 6.0, 6.0)); // overlap
    check2(&a, &sq(1.0, 1.0, 2.0, 2.0)); // nested
    check2(&a, &sq(10.0, 10.0, 12.0, 12.0)); // disjoint
    check2(&a, &sq(4.0, 0.0, 8.0, 4.0)); // edge-touch
    check2(&a, &sq(4.0, 4.0, 8.0, 8.0)); // corner-touch
    check2(&a, &sq(0.0, 0.0, 4.0, 4.0)); // equal
}

#[test]
fn point_polygon() {
    let a = sq(0.0, 0.0, 4.0, 4.0);
    check2(&a, &pt(2.0, 2.0)); // inside
    check2(&a, &pt(4.0, 2.0)); // on edge
    check2(&a, &pt(0.0, 0.0)); // on vertex
    check2(&a, &pt(9.0, 9.0)); // outside
}

#[test]
fn point_point() {
    check2(&pt(1.0, 1.0), &pt(1.0, 1.0)); // equal
    check2(&pt(1.0, 1.0), &pt(2.0, 2.0)); // different
}

#[test]
fn line_line() {
    check2(&line(0.0, 0.0, 4.0, 4.0), &line(0.0, 4.0, 4.0, 0.0)); // cross
    check2(&line(0.0, 0.0, 4.0, 0.0), &line(4.0, 0.0, 8.0, 0.0)); // touch end (collinear)
    check2(&line(0.0, 0.0, 4.0, 0.0), &line(2.0, 0.0, 6.0, 0.0)); // collinear overlap
    check2(&line(0.0, 0.0, 4.0, 0.0), &line(0.0, 1.0, 4.0, 1.0)); // parallel disjoint
    check2(&line(0.0, 0.0, 4.0, 4.0), &line(0.0, 0.0, 4.0, 4.0)); // equal
    check2(&line(0.0, 0.0, 4.0, 0.0), &line(2.0, 0.0, 2.0, 4.0)); // T-touch
}

#[test]
fn line_polygon() {
    let a = sq(0.0, 0.0, 4.0, 4.0);
    check2(&a, &line(2.0, 2.0, 3.0, 3.0)); // line inside
    check2(&a, &line(-1.0, 2.0, 5.0, 2.0)); // line crossing through
    check2(&a, &line(4.0, 0.0, 4.0, 4.0)); // line on edge
    check2(&a, &line(5.0, 0.0, 9.0, 4.0)); // line outside
    check2(&a, &line(2.0, 2.0, 6.0, 2.0)); // interior to exterior
}

#[test]
fn point_line() {
    let l = line(0.0, 0.0, 4.0, 0.0);
    check2(&l, &pt(2.0, 0.0)); // on interior
    check2(&l, &pt(0.0, 0.0)); // on boundary endpoint
    check2(&l, &pt(2.0, 1.0)); // off
}

#[test]
fn polygon_with_hole() {
    let outer = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ];
    let hole = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)];
    let donut = Geometry::Polygon(Polygon::new(
        LineString::from(outer.to_vec()),
        vec![LineString::from(hole.to_vec())],
    ));
    check2(&donut, &pt(5.0, 5.0)); // in the hole → outside
    check2(&donut, &pt(1.0, 1.0)); // in the ring body
    check2(&donut, &sq(4.0, 4.0, 6.0, 6.0)); // square inside the hole
}

#[test]
fn adversarial_areal_relate_matches_geo() {
    let holed = |shell: &[(f64, f64)], hole: &[(f64, f64)]| {
        Geometry::Polygon(Polygon::new(
            LineString::from(shell.to_vec()),
            vec![LineString::from(hole.to_vec())],
        ))
    };
    let shell = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ];
    let hole = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)];
    let donut = holed(&shell, &hole);

    // outer contains inner (interiors overlap; inner is within outer)
    check2(&sq(0.0, 0.0, 10.0, 10.0), &sq(3.0, 3.0, 7.0, 7.0));
    // a thin horizontal bar crossed by a thin vertical bar: interiors overlap in
    // a tiny region and no vertex of one lies inside the other — this is the case
    // that pure interior-sampling could miss but a ring-crossing check catches.
    check2(&sq(0.0, 0.0, 10.0, 1.0), &sq(4.0, -5.0, 5.0, 5.0));
    // square exactly filling the hole touches the donut along the hole ring only
    check2(&donut, &sq(3.0, 3.0, 7.0, 7.0));
    // square overlapping both the hole and the solid ring (proper crossing)
    check2(&donut, &sq(5.0, 5.0, 12.0, 12.0));
    // a small square deep in the solid ring is within the donut
    check2(&donut, &sq(0.5, 0.5, 2.5, 2.5));
    // identical donuts are topologically equal
    check2(&donut, &holed(&shell, &hole));
    // two donuts whose holes overlap but bodies also overlap
    let donut2 = holed(
        &[
            (5.0, 5.0),
            (15.0, 5.0),
            (15.0, 15.0),
            (5.0, 15.0),
            (5.0, 5.0),
        ],
        &[
            (8.0, 8.0),
            (12.0, 8.0),
            (12.0, 12.0),
            (8.0, 12.0),
            (8.0, 8.0),
        ],
    );
    check2(&donut, &donut2);
}

#[test]
fn multipolygon_cases() {
    let mp = Geometry::MultiPolygon(MultiPolygon::new(vec![
        Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (2.0, 0.0),
                (2.0, 2.0),
                (0.0, 2.0),
                (0.0, 0.0),
            ]),
            vec![],
        ),
        Polygon::new(
            LineString::from(vec![
                (5.0, 5.0),
                (7.0, 5.0),
                (7.0, 7.0),
                (5.0, 7.0),
                (5.0, 5.0),
            ]),
            vec![],
        ),
    ]));
    check2(&mp, &pt(1.0, 1.0)); // inside first
    check2(&mp, &pt(6.0, 6.0)); // inside second
    check2(&mp, &pt(3.0, 3.0)); // between → outside
    check2(&mp, &sq(1.0, 1.0, 6.0, 6.0)); // spans both
}

// ---- randomized stress vs geo ----

struct Lcg(u64);
impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
    fn n(&mut self, m: i64) -> f64 {
        (self.next_u32() as i64 % m) as f64
    }
}

fn rrect(r: &mut Lcg, g: i64) -> Geometry {
    let (x, y) = (r.n(g), r.n(g));
    let (w, h) = (1.0 + r.n(g / 2), 1.0 + r.n(g / 2));
    sq(x, y, x + w, y + h)
}
fn rline(r: &mut Lcg, g: i64) -> Geometry {
    // avoid degenerate zero-length lines (invalid input; geo's handling of them
    // is undefined and not something we conform to)
    let (x0, y0) = (r.n(g), r.n(g));
    let (mut x1, mut y1) = (r.n(g), r.n(g));
    if x0 == x1 && y0 == y1 {
        x1 += 1.0;
        y1 += 1.0;
    }
    line(x0, y0, x1, y1)
}
fn rpt(r: &mut Lcg, g: i64) -> Geometry {
    pt(r.n(g), r.n(g))
}

#[test]
fn stress_polygon_polygon() {
    let mut r = Lcg(0xA11CE);
    for _ in 0..500 {
        check2(&rrect(&mut r, 8), &rrect(&mut r, 8));
    }
}

#[test]
fn stress_point_and_line_vs_polygon() {
    let mut r = Lcg(0xB0B);
    for _ in 0..400 {
        let poly = rrect(&mut r, 8);
        check2(&poly, &rpt(&mut r, 9));
        check2(&poly, &rline(&mut r, 9));
    }
}

#[test]
fn stress_line_line_and_misc() {
    let mut r = Lcg(0xCAFE);
    for _ in 0..400 {
        check2(&rline(&mut r, 6), &rline(&mut r, 6));
    }
    let mut r = Lcg(0xD00D);
    for _ in 0..200 {
        let a = ls(&[(r.n(8), r.n(8)), (r.n(8), r.n(8)), (r.n(8), r.n(8))]);
        check2(&a, &rrect(&mut r, 8));
    }
}

fn rtri(r: &mut Lcg, g: i64) -> Option<Geometry> {
    let p = [(r.n(g), r.n(g)), (r.n(g), r.n(g)), (r.n(g), r.n(g))];
    let area2 = (p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1);
    if area2.abs() < 1.0 {
        return None; // skip degenerate/near-collinear triangles
    }
    Some(poly(&[p[0], p[1], p[2], p[0]]))
}

#[test]
fn stress_triangles_vs_areas() {
    let mut r = Lcg(0x7711);
    let mut n = 0;
    while n < 400 {
        let Some(t) = rtri(&mut r, 9) else { continue };
        check2(&t, &rrect(&mut r, 9));
        if let Some(t2) = rtri(&mut r, 9) {
            check2(&t, &t2);
        }
        n += 1;
    }
}
