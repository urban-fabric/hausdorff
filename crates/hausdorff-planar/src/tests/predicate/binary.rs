//! The binary predicates, from three directions at once.
//!
//! These were three separate files — `predicates.rs`, `intersects.rs` and
//! `r2/contains.rs` — holding two, one and two tests respectively, all answering
//! the one question the suite index already bundled them under: do the predicates
//! agree with each other and with the matrix they dispatch to? Three files of one
//! or two tests is how `Contains` ended up covered in two places that neither
//! referenced.
//!
//! `SpatialPredicate::eval` must agree with the DE-9IM `relate` matrix it
//! dispatches to — in particular the `Intersects`/`Disjoint` direct path and the
//! point-in-areal `Contains`/`Within` fast path must equal the full matrix — and
//! `transpose` must hold (`a P b` ⟺ `b P.transpose() a`).

use crate::space::r2::SpatialPredicate;
use crate::space::r2::relate;
use crate::space::r2::{
    Contains, Coord, Geometry, Intersects, LineString, MultiPolygon, Point, Polygon, Rect,
};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn f(&mut self, half: f64) -> f64 {
        (f64::from((self.next() & 0x00FF_FFFF) as u32) / f64::from(0x00FF_FFFF_u32) * 2.0 - 1.0)
            * half
    }
}

fn make(rng: &mut Rng, kind: u64) -> Geometry {
    match kind % 4 {
        0 => Geometry::Point(Point::new(rng.f(10.0), rng.f(10.0))),
        1 => Geometry::LineString(LineString::from(vec![
            (rng.f(10.0), rng.f(10.0)),
            (rng.f(10.0), rng.f(10.0)),
            (rng.f(10.0), rng.f(10.0)),
        ])),
        2 => {
            let (x, y) = (rng.f(8.0), rng.f(8.0));
            let (w, h) = (1.0 + rng.f(4.0).abs(), 1.0 + rng.f(4.0).abs());
            Geometry::Polygon(Polygon::new(
                LineString::from(vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h), (x, y)]),
                vec![],
            ))
        }
        _ => Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (rng.f(10.0), rng.f(10.0)),
                (rng.f(10.0), rng.f(10.0)),
                (rng.f(10.0), rng.f(10.0)),
            ]),
            vec![],
        )),
    }
}

const ALL: [SpatialPredicate; 10] = [
    SpatialPredicate::Intersects,
    SpatialPredicate::Contains,
    SpatialPredicate::Within,
    SpatialPredicate::Covers,
    SpatialPredicate::CoveredBy,
    SpatialPredicate::Disjoint,
    SpatialPredicate::Touches,
    SpatialPredicate::Crosses,
    SpatialPredicate::Equals,
    SpatialPredicate::Overlaps,
];

#[test]
fn eval_agrees_with_relate_and_transpose_holds() {
    let mut rng = Rng(0x5DEE_CE66_D1B5_4A32);
    for _ in 0..1500 {
        let ka = rng.next();
        let a = make(&mut rng, ka);
        let kb = rng.next();
        let b = make(&mut rng, kb);
        let m = relate(&a, &b);

        assert_eq!(
            SpatialPredicate::Intersects.eval(&a, &b),
            m.is_intersects(),
            "intersects"
        );
        assert_eq!(
            SpatialPredicate::Disjoint.eval(&a, &b),
            m.is_disjoint(),
            "disjoint"
        );
        assert_eq!(
            SpatialPredicate::Contains.eval(&a, &b),
            m.is_contains(),
            "contains"
        );
        assert_eq!(
            SpatialPredicate::Within.eval(&a, &b),
            m.is_within(),
            "within"
        );
        assert_eq!(
            SpatialPredicate::Covers.eval(&a, &b),
            m.is_covers(),
            "covers"
        );
        assert_eq!(
            SpatialPredicate::CoveredBy.eval(&a, &b),
            m.is_coveredby(),
            "coveredby"
        );
        assert_eq!(
            SpatialPredicate::Touches.eval(&a, &b),
            m.is_touches(),
            "touches"
        );
        assert_eq!(
            SpatialPredicate::Equals.eval(&a, &b),
            m.is_equal_topo(),
            "equals"
        );

        for p in ALL {
            assert_eq!(
                p.eval(&a, &b),
                p.transpose().eval(&b, &a),
                "transpose {p:?}"
            );
        }
    }
}

#[test]
fn geometry_contains_point_is_areal_only() {
    let sq = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    assert!(sq.contains(&Coord::new(2.0, 2.0)));
    assert!(!sq.contains(&Coord::new(9.0, 9.0)));
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (4.0, 0.0)]));
    assert!(!line.contains(&Coord::new(2.0, 0.0)));
}

// ---- folded in from the former `r2/contains.rs` ----

fn holed_square() -> Polygon {
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
    let s = holed_square();
    assert!(s.contains(&Point::new(0.5, 0.5)));
    assert!(!s.contains(&Point::new(2.0, 2.0)));
    assert!(!s.contains(&Point::new(9.0, 9.0)));
}

#[test]
fn multipolygon_and_rect() {
    let mp = MultiPolygon::new(vec![holed_square()]);
    assert!(mp.contains(&Coord::new(0.5, 0.5)));
    let r = Rect::new((0.0, 0.0), (2.0, 2.0));
    assert!(r.contains(&Point::new(1.0, 1.0)));
    assert!(!r.contains(&Point::new(2.0, 2.0)));
}

// ---- folded in from the former `intersects.rs` ----

fn hpoly(coords: &[(f64, f64)]) -> Geometry {
    Geometry::Polygon(Polygon::new(LineString::from(coords.to_vec()), vec![]))
}
fn hpt(x: f64, y: f64) -> Geometry {
    Geometry::Point(Point::new(x, y))
}
fn hline(c: &[(f64, f64)]) -> Geometry {
    Geometry::LineString(LineString::from(c.to_vec()))
}

#[test]
fn intersects_matches_golden() {
    let sq_a = hpoly(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]);
    // (geometry, whether it intersects the 4×4 square at the origin)
    let cases: Vec<(Geometry, bool)> = vec![
        (
            hpoly(&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]),
            true,
        ), // overlapping
        (
            hpoly(&[(9.0, 9.0), (10.0, 9.0), (10.0, 10.0), (9.0, 10.0)]),
            false,
        ), // disjoint
        (
            hpoly(&[(4.0, 0.0), (8.0, 0.0), (8.0, 4.0), (4.0, 4.0)]),
            true,
        ), // shared edge
        (hpt(2.0, 2.0), true),                         // interior point
        (hpt(9.0, 9.0), false),                        // exterior point
        (hline(&[(2.0, -1.0), (2.0, 5.0)]), true),     // crossing line
        (hline(&[(-5.0, -5.0), (-1.0, -1.0)]), false), // disjoint line
    ];
    for (hb, want) in cases {
        let ours = sq_a.intersects(&hb);
        assert_eq!(ours, want, "intersects mismatch for {hb:?}");
    }
}
