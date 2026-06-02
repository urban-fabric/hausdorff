//! `SpatialPredicate::eval` must agree with the DE-9IM `relate` matrix it
//! dispatches to — in particular the `Intersects`/`Disjoint` direct path and the
//! point-in-areal `Contains`/`Within` fast path must equal the full matrix — and
//! `transpose` must hold (`a P b` ⟺ `b P.transpose() a`).

use hausdorff::space::r2::SpatialPredicate;
use hausdorff::space::r2::relate::relate;
use hausdorff::{Contains, Coord, Geometry, LineString, Point, Polygon};

// Continuous coordinates (not an integer grid), so points never land exactly on
// an edge — that boundary case is the one documented divergence between the
// ray-cast fast path and the matrix, and it is measure-zero here.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
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
                LineString::from(vec![
                    (x, y),
                    (x + w, y),
                    (x + w, y + h),
                    (x, y + h),
                    (x, y),
                ]),
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

        // The predicates whose `is_*` take no dimension argument: eval must equal
        // the matrix exactly, including the Intersects/Disjoint and point-in-areal
        // Contains/Within short-circuits.
        assert_eq!(SpatialPredicate::Intersects.eval(&a, &b), m.is_intersects(), "intersects");
        assert_eq!(SpatialPredicate::Disjoint.eval(&a, &b), m.is_disjoint(), "disjoint");
        assert_eq!(SpatialPredicate::Contains.eval(&a, &b), m.is_contains(), "contains");
        assert_eq!(SpatialPredicate::Within.eval(&a, &b), m.is_within(), "within");
        assert_eq!(SpatialPredicate::Covers.eval(&a, &b), m.is_covers(), "covers");
        assert_eq!(SpatialPredicate::CoveredBy.eval(&a, &b), m.is_coveredby(), "coveredby");
        assert_eq!(SpatialPredicate::Touches.eval(&a, &b), m.is_touches(), "touches");
        assert_eq!(SpatialPredicate::Equals.eval(&a, &b), m.is_equal_topo(), "equals");

        // `a P b` ⟺ `b P.transpose() a` for every predicate (covers Crosses /
        // Overlaps symmetry too).
        for p in ALL {
            assert_eq!(p.eval(&a, &b), p.transpose().eval(&b, &a), "transpose {p:?}");
        }
    }
}

#[test]
fn geometry_contains_point_is_areal_only() {
    let sq = Geometry::Polygon(Polygon::new(
        LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]),
        vec![],
    ));
    assert!(sq.contains(&Coord::new(2.0, 2.0)));
    assert!(!sq.contains(&Coord::new(9.0, 9.0)));
    // A non-areal geometry contains no point under this (areal) predicate.
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (4.0, 0.0)]));
    assert!(!line.contains(&Coord::new(2.0, 0.0)));
}
