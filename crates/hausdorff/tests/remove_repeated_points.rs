//! Conformance for the native RemoveRepeatedPoints: the result must equal the
//! input with consecutive duplicates collapsed (full first-occurrence dedup for
//! multi-points), the `_mut` variant must match the owning form, and rings must
//! stay closed. Expected outputs are written out independently (the operation is
//! a trivial dedup); no third-party geometry crate.

use hausdorff::space::r2::RemoveRepeatedPoints;
use hausdorff::{Geometry, LineString, MultiLineString, MultiPoint, MultiPolygon, Polygon};

fn ls(coords: &[(f64, f64)]) -> LineString {
    LineString::from(coords.to_vec())
}

/// Ours must equal the independently-specified `expected` cleanup, and the
/// in-place variant must match the owning one.
fn check(name: &str, geom: &Geometry, expected: &Geometry) {
    let ours = geom.remove_repeated_points();
    assert_eq!(&ours, expected, "{name}: owning variant");

    let mut in_place = geom.clone();
    in_place.remove_repeated_points_mut();
    assert_eq!(in_place, ours, "{name}: _mut variant");

    // The operation is idempotent: a cleaned geometry is a fixed point.
    assert_eq!(ours.remove_repeated_points(), ours, "{name}: idempotent");
}

#[test]
fn matches_expected_across_fixtures() {
    let dup_line = ls(&[
        (0.0, 0.0),
        (1.0, 1.0),
        (1.0, 1.0),
        (1.0, 1.0),
        (2.0, 2.0),
        (2.0, 2.0),
        (0.0, 0.0),
    ]);
    let dedup_line = ls(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (0.0, 0.0)]);
    let dup_ring = ls(&[
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (0.0, 4.0),
        (0.0, 0.0),
    ]);
    let dedup_ring = ls(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]);

    let cases: Vec<(&str, Geometry, Geometry)> = vec![
        (
            "line string",
            Geometry::LineString(dup_line.clone()),
            Geometry::LineString(dedup_line.clone()),
        ),
        (
            "line string without duplicates",
            Geometry::LineString(ls(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)])),
            Geometry::LineString(ls(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)])),
        ),
        (
            "multi point full dedup",
            Geometry::MultiPoint(MultiPoint::from(vec![
                (0.0, 0.0),
                (1.0, 1.0),
                (1.0, 1.0),
                (2.0, 2.0),
                (0.0, 0.0),
            ])),
            // Full first-occurrence dedup drops the repeated (1,1) and trailing (0,0).
            Geometry::MultiPoint(MultiPoint::from(vec![(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)])),
        ),
        (
            "polygon with duplicated ring points",
            Geometry::Polygon(Polygon::new(
                dup_ring.clone(),
                vec![ls(&[
                    (1.0, 1.0),
                    (2.0, 1.0),
                    (2.0, 1.0),
                    (2.0, 2.0),
                    (1.0, 1.0),
                ])],
            )),
            Geometry::Polygon(Polygon::new(
                dedup_ring.clone(),
                vec![ls(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 1.0)])],
            )),
        ),
        (
            "multi line string",
            Geometry::MultiLineString(MultiLineString::new(vec![
                dup_line.clone(),
                ls(&[(5.0, 5.0), (5.0, 5.0), (6.0, 6.0)]),
            ])),
            Geometry::MultiLineString(MultiLineString::new(vec![
                dedup_line.clone(),
                ls(&[(5.0, 5.0), (6.0, 6.0)]),
            ])),
        ),
        (
            "multi polygon",
            Geometry::MultiPolygon(MultiPolygon::new(vec![
                Polygon::new(dup_ring, vec![]),
                Polygon::new(
                    ls(&[
                        (10.0, 10.0),
                        (11.0, 11.0),
                        (11.0, 11.0),
                        (10.0, 12.0),
                        (10.0, 10.0),
                    ]),
                    vec![],
                ),
            ])),
            Geometry::MultiPolygon(MultiPolygon::new(vec![
                Polygon::new(dedup_ring, vec![]),
                Polygon::new(
                    ls(&[(10.0, 10.0), (11.0, 11.0), (10.0, 12.0), (10.0, 10.0)]),
                    vec![],
                ),
            ])),
        ),
        (
            "point is unchanged",
            Geometry::Point(hausdorff::Point::new(3.0, 4.0)),
            Geometry::Point(hausdorff::Point::new(3.0, 4.0)),
        ),
        (
            "geometry collection recurses",
            Geometry::GeometryCollection(vec![
                Geometry::LineString(dup_line),
                Geometry::MultiPoint(MultiPoint::from(vec![(1.0, 1.0), (1.0, 1.0)])),
            ]),
            Geometry::GeometryCollection(vec![
                Geometry::LineString(dedup_line),
                Geometry::MultiPoint(MultiPoint::from(vec![(1.0, 1.0)])),
            ]),
        ),
    ];

    for (name, geom, expected) in &cases {
        check(name, geom, expected);
    }
}

#[test]
fn deduped_polygon_rings_stay_closed() {
    let polygon = Polygon::new(
        ls(&[
            (0.0, 0.0),
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (4.0, 4.0),
            (0.0, 4.0),
        ]),
        vec![],
    );
    let cleaned = polygon.remove_repeated_points();
    assert!(cleaned.exterior().is_closed());
    assert_eq!(cleaned.exterior().0.len(), 5);
}
