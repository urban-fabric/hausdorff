//! Conformance for the native RemoveRepeatedPoints against geo's: the same
//! input must produce coordinate-identical output (consecutive dedup for lines
//! and rings, full first-occurrence dedup for multi-points), and the `_mut`
//! variants must agree with the owning forms.

mod common;

use geo::RemoveRepeatedPoints as _;
use hausdorff::space::r2::RemoveRepeatedPoints;
use hausdorff::{Geometry, LineString, MultiLineString, MultiPoint, MultiPolygon, Polygon};

fn ls(coords: &[(f64, f64)]) -> LineString {
    LineString::from(coords.to_vec())
}

/// Ours and geo's cleanup of the same geometry must be identical, and the
/// in-place variant must match the owning one.
fn check(name: &str, geom: &Geometry) {
    let ours = geom.remove_repeated_points();
    let geo_result = common::to_geo_geometry(geom).remove_repeated_points();
    assert_eq!(
        common::to_geo_geometry(&ours),
        geo_result,
        "{name}: owning variant"
    );

    let mut in_place = geom.clone();
    in_place.remove_repeated_points_mut();
    assert_eq!(in_place, ours, "{name}: _mut variant");
}

#[test]
fn matches_geo_across_fixtures() {
    let dup_line = ls(&[
        (0.0, 0.0),
        (1.0, 1.0),
        (1.0, 1.0),
        (1.0, 1.0),
        (2.0, 2.0),
        (2.0, 2.0),
        (0.0, 0.0),
    ]);
    let dup_ring = ls(&[
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (0.0, 4.0),
        (0.0, 0.0),
    ]);
    let fixtures: Vec<(&str, Geometry)> = vec![
        ("line string", Geometry::LineString(dup_line.clone())),
        (
            "line string without duplicates",
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
        ),
        (
            "multi line string",
            Geometry::MultiLineString(MultiLineString::new(vec![
                dup_line.clone(),
                ls(&[(5.0, 5.0), (5.0, 5.0), (6.0, 6.0)]),
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
        ),
        (
            "point is unchanged",
            Geometry::Point(hausdorff::Point::new(3.0, 4.0)),
        ),
        (
            "geometry collection recurses",
            Geometry::GeometryCollection(vec![
                Geometry::LineString(dup_line),
                Geometry::MultiPoint(MultiPoint::from(vec![(1.0, 1.0), (1.0, 1.0)])),
            ]),
        ),
    ];

    for (name, geom) in &fixtures {
        check(name, geom);
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
