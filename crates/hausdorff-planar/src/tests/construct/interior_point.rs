//! Conformance for the native InteriorPoint: the deterministic cases (the
//! widest-scanline midpoint, vertex picks, centers) must match golden values
//! captured from the reference; for the rest the contract must hold — the
//! returned point intersects the geometry, and lies in the interior of areal
//! geometries with non-empty interiors.

use crate::space::r2::InteriorPoint;
use crate::space::r2::{
    Contains, Geometry, Intersects, Line, LineString, MultiLineString, MultiPoint, MultiPolygon,
    Point, Polygon, Rect, Triangle,
};

fn ls(coords: &[(f64, f64)]) -> LineString {
    LineString::from(coords.to_vec())
}

fn poly(rings: &[&[(f64, f64)]]) -> Polygon {
    Polygon::new(
        ls(rings[0]),
        rings[1..].iter().map(|r| ls(r)).collect::<Vec<_>>(),
    )
}

#[test]
fn deterministic_cases_match_golden() {
    // Golden interior points captured from the reference (geo) at conversion time.
    const GOLDEN: [Option<(f64, f64)>; 12] = [
        Some((1.0, 2.0)),
        Some((2.0, 2.0)),
        Some((2.0, 5.0)),
        Some((8.5, 3.5)),
        Some((2.0, 2.0)),
        Some((0.0, 0.0)),
        Some((1.0, 0.0)),
        Some((3.0, 2.0)),
        Some((7.0, -3.0)),
        Some((0.0, 0.0)),
        Some((2.0, 1.0)),
        Some((1.0, 1.0)),
    ];
    let fixtures: Vec<(&str, Geometry)> = vec![
        (
            "rhombus",
            Geometry::Polygon(poly(&[&[(-2.0, 1.0), (1.0, 3.0), (4.0, 1.0), (1.0, -1.0)]])),
        ),
        (
            "square",
            Geometry::Polygon(poly(&[&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]])),
        ),
        (
            "rectangle with centered hole",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
                &[(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0)],
            ])),
        ),
        (
            "multi polygon picks the widest member",
            Geometry::MultiPolygon(MultiPolygon::new(vec![
                poly(&[&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]]),
                poly(&[&[(5.0, 0.0), (12.0, 0.0), (12.0, 7.0), (5.0, 7.0)]]),
            ])),
        ),
        (
            "line string interior vertex",
            Geometry::LineString(ls(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (9.0, 9.0)])),
        ),
        (
            "two point line string",
            Geometry::LineString(ls(&[(0.0, 0.0), (1.0, 1.0)])),
        ),
        (
            "multi line string",
            Geometry::MultiLineString(MultiLineString::new(vec![
                ls(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]),
                ls(&[(0.0, 5.0), (1.0, 5.0), (2.0, 5.0)]),
            ])),
        ),
        (
            "multi point",
            Geometry::MultiPoint(MultiPoint::from(vec![(5.0, 1.0), (1.0, 3.0), (3.0, 2.0)])),
        ),
        ("point", Geometry::Point(Point::new(7.0, -3.0))),
        ("line", Geometry::Line(Line::new((0.0, 0.0), (2.0, 4.0)))),
        ("rect", Geometry::Rect(Rect::new((0.0, 0.0), (4.0, 2.0)))),
        (
            "triangle",
            Geometry::Triangle(Triangle::new(
                crate::space::r2::Coord::new(0.0, 0.0),
                crate::space::r2::Coord::new(3.0, 0.0),
                crate::space::r2::Coord::new(0.0, 3.0),
            )),
        ),
    ];
    for ((name, geom), golden) in fixtures.iter().zip(GOLDEN) {
        assert_eq!(geom.interior_point().map(|p| p.x_y()), golden, "{name}");
    }
}

#[test]
fn rhombus_matches_geo_doc_value() {
    let rhombus = poly(&[&[(-2.0, 1.0), (1.0, 3.0), (4.0, 1.0), (1.0, -1.0)]]);
    assert_eq!(rhombus.interior_point(), Some(Point::new(1.0, 2.0)));
}

#[test]
fn contract_point_is_on_the_geometry() {
    let shapes: Vec<Geometry> = vec![
        Geometry::Polygon(poly(&[&[
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 6.0),
            (3.0, 3.0),
            (0.0, 6.0),
        ]])),
        Geometry::Polygon(poly(&[
            &[(0.0, 0.0), (12.0, 0.0), (12.0, 10.0), (0.0, 10.0)],
            &[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)],
            &[(5.0, 1.0), (8.0, 2.0), (7.0, 4.0), (4.0, 3.0)],
            &[(2.0, 5.0), (6.0, 6.0), (5.0, 9.0), (1.0, 8.0)],
        ])),
        Geometry::MultiPolygon(MultiPolygon::new(vec![
            poly(&[&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]]),
            poly(&[&[(8.0, 8.0), (9.0, 8.0), (9.0, 9.0), (8.0, 9.0)]]),
        ])),
        Geometry::LineString(ls(&[(0.0, 0.0), (3.0, 7.0), (-2.0, 5.0), (4.0, 4.0)])),
        Geometry::MultiPoint(MultiPoint::from(vec![(0.0, 0.0), (10.0, 10.0)])),
        Geometry::GeometryCollection(vec![
            Geometry::Point(Point::new(0.0, 0.0)),
            Geometry::Polygon(poly(&[&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]])),
        ]),
    ];

    for geom in &shapes {
        let pt = geom.interior_point().expect("non-empty geometry");
        assert!(
            geom.intersects(&Geometry::Point(pt)),
            "{geom:?} does not intersect its interior point {pt:?}"
        );
    }
}

#[test]
fn contract_areal_interior_points_are_contained() {
    let polygons = vec![
        poly(&[&[(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (3.0, 3.0), (0.0, 6.0)]]),
        poly(&[
            &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
            &[(1.0, 1.0), (9.0, 1.0), (9.0, 9.0), (1.0, 9.0)],
        ]),
        poly(&[&[(0.0, 0.0), (100.0, 0.0), (100.0, 0.5), (0.0, 0.5)]]),
    ];
    for polygon in &polygons {
        let pt = polygon.interior_point().expect("non-empty polygon");
        assert!(
            polygon.contains(&pt),
            "{polygon:?} interior point {pt:?} not strictly inside"
        );
    }
}

#[test]
fn empty_geometries_yield_none() {
    assert_eq!(LineString::new(vec![]).interior_point(), None);
    assert_eq!(MultiPoint::new(vec![]).interior_point(), None);
    assert_eq!(MultiPolygon::new(vec![]).interior_point(), None);
    assert_eq!(MultiLineString::new(vec![]).interior_point(), None);
}
