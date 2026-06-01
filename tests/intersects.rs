//! Native Intersects predicate conformance vs geo's Intersects.

use geo::Intersects as _;
use hausdorff::r2::Intersects;
use hausdorff::{Geometry, LineString, Point, Polygon};

fn hpoly(coords: &[(f64, f64)]) -> Geometry {
    Geometry::Polygon(Polygon::new(LineString::from(coords.to_vec()), vec![]))
}
fn gpoly(coords: &[(f64, f64)]) -> geo::Geometry<f64> {
    geo::Geometry::Polygon(geo::Polygon::new(
        geo::LineString(coords.iter().map(|&(x, y)| geo::Coord { x, y }).collect()),
        vec![],
    ))
}
fn hpt(x: f64, y: f64) -> Geometry {
    Geometry::Point(Point::new(x, y))
}
fn gpt(x: f64, y: f64) -> geo::Geometry<f64> {
    geo::Geometry::Point(geo::Point::new(x, y))
}
fn hline(c: &[(f64, f64)]) -> Geometry {
    Geometry::LineString(LineString::from(c.to_vec()))
}
fn gline(c: &[(f64, f64)]) -> geo::Geometry<f64> {
    geo::Geometry::LineString(geo::LineString(
        c.iter().map(|&(x, y)| geo::Coord { x, y }).collect(),
    ))
}

#[test]
fn intersects_matches_geo() {
    let sq_a = (
        hpoly(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
        gpoly(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
    );
    let cases: Vec<(Geometry, geo::Geometry<f64>)> = vec![
        // overlapping square
        (hpoly(&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]),
         gpoly(&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)])),
        // disjoint square
        (hpoly(&[(9.0, 9.0), (10.0, 9.0), (10.0, 10.0), (9.0, 10.0)]),
         gpoly(&[(9.0, 9.0), (10.0, 9.0), (10.0, 10.0), (9.0, 10.0)])),
        // touching at an edge
        (hpoly(&[(4.0, 0.0), (8.0, 0.0), (8.0, 4.0), (4.0, 4.0)]),
         gpoly(&[(4.0, 0.0), (8.0, 0.0), (8.0, 4.0), (4.0, 4.0)])),
        // point inside / outside
        (hpt(2.0, 2.0), gpt(2.0, 2.0)),
        (hpt(9.0, 9.0), gpt(9.0, 9.0)),
        // line crossing / clear
        (hline(&[(2.0, -1.0), (2.0, 5.0)]), gline(&[(2.0, -1.0), (2.0, 5.0)])),
        (hline(&[(-5.0, -5.0), (-1.0, -1.0)]), gline(&[(-5.0, -5.0), (-1.0, -1.0)])),
    ];
    for (hb, gb) in cases {
        let ours = sq_a.0.intersects(&hb);
        let theirs = sq_a.1.intersects(&gb);
        assert_eq!(ours, theirs, "intersects mismatch: ours={ours} geo={theirs} for {gb:?}");
    }
}
