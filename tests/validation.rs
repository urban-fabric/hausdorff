//! Conformance for the native OGC validity checks against geo's Validation.
//! The `is_valid` verdict must agree with geo across a battery of valid and
//! invalid fixtures, and the typed error variants must match expectations on
//! the canonical invalid shapes.

mod common;

use geo::algorithm::Validation as _;
use hausdorff::space::r2::validation::{
    CoordIndex, GeometryIndex, InvalidLineString, InvalidMultiPolygon, InvalidPolygon, RingRole,
    Validation,
};
use hausdorff::{Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Polygon};

fn ls(coords: &[(f64, f64)]) -> LineString {
    LineString::from(coords.to_vec())
}

fn poly(rings: &[&[(f64, f64)]]) -> Polygon {
    Polygon::new(
        ls(rings[0]),
        rings[1..].iter().map(|r| ls(r)).collect::<Vec<_>>(),
    )
}

/// Both libraries must agree on the validity verdict for `geom`.
fn check_verdict(name: &str, geom: &Geometry) {
    let ours = geom.is_valid();
    let geo_verdict = common::to_geo_geometry(geom).is_valid();
    assert_eq!(
        ours, geo_verdict,
        "{name}: hausdorff said {ours}, geo said {geo_verdict}"
    );
}

fn fixtures() -> Vec<(&'static str, Geometry, bool)> {
    vec![
        (
            "valid square",
            Geometry::Polygon(poly(&[&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]])),
            true,
        ),
        (
            "valid concave polygon",
            Geometry::Polygon(poly(&[&[
                (0.0, 0.0),
                (6.0, 0.0),
                (6.0, 6.0),
                (3.0, 3.0),
                (0.0, 6.0),
            ]])),
            true,
        ),
        (
            "valid complex polygon with multiple holes",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (12.0, 0.0), (12.0, 10.0), (0.0, 10.0)],
                &[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)],
                &[(5.0, 1.0), (8.0, 2.0), (7.0, 4.0), (4.0, 3.0)],
                &[(2.0, 5.0), (6.0, 6.0), (5.0, 9.0), (1.0, 8.0)],
                &[(9.0, 5.0), (11.0, 5.0), (11.0, 9.0), (9.0, 9.0)],
            ])),
            true,
        ),
        (
            "bowtie polygon",
            Geometry::Polygon(poly(&[&[(0.0, 0.0), (4.0, 0.0), (0.0, 2.0), (4.0, 2.0)]])),
            false,
        ),
        (
            "spike polygon",
            Geometry::Polygon(poly(&[&[
                (0.0, 0.0),
                (4.0, 0.0),
                (4.0, 4.0),
                (2.0, 4.0),
                (2.0, 6.0),
                (2.0, 4.0),
                (0.0, 4.0),
            ]])),
            false,
        ),
        (
            "hole outside shell",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
                &[(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0)],
            ])),
            false,
        ),
        (
            "hole partially outside shell",
            Geometry::Polygon(poly(&[
                &[(0.5, 0.5), (3.0, 0.5), (3.0, 2.5), (0.5, 2.5)],
                &[(1.0, 1.0), (1.0, 2.0), (2.5, 2.0), (3.5, 1.0)],
            ])),
            false,
        ),
        (
            "holes touching at a point",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
                &[(1.0, 2.0), (2.0, 1.0), (3.0, 2.0), (2.0, 3.0)],
                &[(3.0, 2.0), (4.0, 1.0), (5.0, 2.0), (4.0, 3.0)],
            ])),
            true,
        ),
        (
            "holes touching on a line",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
                &[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)],
                &[(3.0, 1.0), (5.0, 1.0), (5.0, 3.0), (3.0, 3.0)],
            ])),
            false,
        ),
        (
            "overlapping holes",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
                &[(1.0, 1.0), (5.0, 1.0), (5.0, 5.0), (1.0, 5.0)],
                &[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)],
            ])),
            false,
        ),
        (
            "nested holes",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
                &[(1.0, 1.0), (1.0, 9.0), (9.0, 9.0), (9.0, 1.0)],
                &[(2.0, 2.0), (2.0, 8.0), (8.0, 8.0), (8.0, 2.0)],
            ])),
            false,
        ),
        (
            "hole touching shell at a point",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (4.0, 1.0), (4.0, 4.0), (0.0, 4.0)],
                &[(0.0, 2.0), (2.0, 1.0), (3.0, 2.0), (2.0, 3.0)],
            ])),
            true,
        ),
        (
            "hole touching shell on a line",
            Geometry::Polygon(poly(&[
                &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
                &[(0.0, 2.0), (0.0, 1.0), (2.0, 1.0), (3.0, 2.0), (2.0, 3.0)],
            ])),
            false,
        ),
        (
            "nan coordinate in polygon",
            Geometry::Polygon(poly(&[&[
                (0.0, 0.0),
                (4.0, 0.0),
                (f64::NAN, 4.0),
                (0.0, 4.0),
            ]])),
            false,
        ),
        (
            "infinite coordinate in line string",
            Geometry::LineString(ls(&[(0.0, 0.0), (f64::INFINITY, 1.0)])),
            false,
        ),
        (
            "unclosed ring is auto-closed and valid",
            Geometry::Polygon(poly(&[&[(0.0, 0.0), (1.0, 1.0), (0.0, 1.0)]])),
            true,
        ),
        (
            "ring with too few distinct points",
            Geometry::Polygon(poly(&[&[(0.0, 0.0), (1.0, 1.0)]])),
            false,
        ),
        (
            "valid line string",
            Geometry::LineString(ls(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])),
            true,
        ),
        (
            "line string with one distinct point",
            Geometry::LineString(ls(&[(0.0, 0.0), (0.0, 0.0)])),
            false,
        ),
        (
            "multi line string with a degenerate member",
            Geometry::MultiLineString(MultiLineString::new(vec![
                ls(&[(0.0, 0.0), (1.0, 1.0)]),
                ls(&[(2.0, 2.0), (2.0, 2.0)]),
            ])),
            false,
        ),
        (
            "valid multi point",
            Geometry::MultiPoint(MultiPoint::from(vec![(0.0, 0.0), (1.0, 1.0)])),
            true,
        ),
        (
            "multi point with nan",
            Geometry::MultiPoint(MultiPoint::from(vec![(0.0, 0.0), (f64::NAN, 1.0)])),
            false,
        ),
        (
            "degenerate line",
            Geometry::Line(Line::new((1.0, 1.0), (1.0, 1.0))),
            false,
        ),
        (
            "valid multi polygon",
            Geometry::MultiPolygon(MultiPolygon::new(vec![
                poly(&[&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]]),
                poly(&[&[(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0)]]),
            ])),
            true,
        ),
        (
            "multi polygon members touching at a point",
            Geometry::MultiPolygon(MultiPolygon::new(vec![
                poly(&[&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]]),
                poly(&[&[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0)]]),
            ])),
            true,
        ),
        (
            "multi polygon members sharing an edge",
            Geometry::MultiPolygon(MultiPolygon::new(vec![
                poly(&[&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]]),
                poly(&[&[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0)]]),
            ])),
            false,
        ),
        (
            "multi polygon members overlapping",
            Geometry::MultiPolygon(MultiPolygon::new(vec![
                poly(&[&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]]),
                poly(&[&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]]),
            ])),
            false,
        ),
        (
            "geometry collection with an invalid member",
            Geometry::GeometryCollection(vec![
                Geometry::Point(hausdorff::Point::new(0.0, 0.0)),
                Geometry::LineString(ls(&[(0.0, 0.0), (0.0, 0.0)])),
            ]),
            false,
        ),
    ]
}

#[test]
fn verdicts_match_geo() {
    for (name, geom, _) in fixtures() {
        check_verdict(name, &geom);
    }
}

#[test]
fn verdicts_match_expectations() {
    for (name, geom, expected) in fixtures() {
        assert_eq!(geom.is_valid(), expected, "{name}");
    }
}

#[test]
fn bowtie_reports_self_intersection() {
    let bowtie = poly(&[&[(0.0, 0.0), (4.0, 0.0), (0.0, 2.0), (4.0, 2.0)]]);
    assert_eq!(
        bowtie.validation_errors(),
        vec![InvalidPolygon::SelfIntersection(RingRole::Exterior)]
    );
}

#[test]
fn hole_outside_shell_reports_containment_failure() {
    let polygon = poly(&[
        &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        &[(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0)],
    ]);
    assert_eq!(
        polygon.validation_errors(),
        vec![InvalidPolygon::InteriorRingNotContainedInExteriorRing(
            RingRole::Interior(0)
        )]
    );
}

#[test]
fn nan_coordinates_report_each_offending_index() {
    let line = ls(&[(f64::NAN, 0.0), (1.0, 1.0), (2.0, f64::NAN)]);
    assert_eq!(
        line.validation_errors(),
        vec![
            InvalidLineString::NonFiniteCoord(CoordIndex(0)),
            InvalidLineString::NonFiniteCoord(CoordIndex(2)),
        ]
    );
}

#[test]
fn multi_polygon_errors_carry_member_indices() {
    let mp = MultiPolygon::new(vec![
        poly(&[&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]]),
        poly(&[&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]]),
    ]);
    assert_eq!(
        mp.check_validation().unwrap_err(),
        InvalidMultiPolygon::ElementsOverlaps(GeometryIndex(0), GeometryIndex(1))
    );
}

#[test]
fn first_error_and_error_list_are_consistent() {
    for (name, geom, _) in fixtures() {
        let errors = geom.validation_errors();
        match geom.check_validation() {
            Ok(()) => assert!(errors.is_empty(), "{name}"),
            Err(first) => assert_eq!(Some(&first), errors.first(), "{name}"),
        }
    }
}
