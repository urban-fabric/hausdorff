//! Native Intersects predicate conformance against golden expected results
//! (a 4×4 square at the origin against overlapping / disjoint / edge-sharing /
//! interior / exterior / crossing configurations). No third-party geometry crate.

use hausdorff::space::r2::Intersects;
use hausdorff::{Geometry, LineString, Point, Polygon};

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
        (hpoly(&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]), true), // overlapping
        (hpoly(&[(9.0, 9.0), (10.0, 9.0), (10.0, 10.0), (9.0, 10.0)]), false), // disjoint
        (hpoly(&[(4.0, 0.0), (8.0, 0.0), (8.0, 4.0), (4.0, 4.0)]), true), // shared edge
        (hpt(2.0, 2.0), true),                                            // interior point
        (hpt(9.0, 9.0), false),                                           // exterior point
        (hline(&[(2.0, -1.0), (2.0, 5.0)]), true),                        // crossing line
        (hline(&[(-5.0, -5.0), (-1.0, -1.0)]), false),                    // disjoint line
    ];
    for (hb, want) in cases {
        let ours = sq_a.intersects(&hb);
        assert_eq!(ours, want, "intersects mismatch for {hb:?}");
    }
}
