use super::*;

fn ls(coords: &[(f64, f64)]) -> LineString {
    LineString::from(coords.to_vec())
}

#[test]
fn simple_polygon_is_valid() {
    let polygon = Polygon::new(ls(&[(0.0, 0.0), (1.0, 1.0), (0.0, 1.0)]), vec![]);
    assert!(polygon.is_valid());
    assert!(polygon.validation_errors().is_empty());
}

#[test]
fn bowtie_polygon_is_self_intersecting() {
    let polygon = Polygon::new(
        ls(&[(0.0, 0.0), (4.0, 0.0), (0.0, 2.0), (4.0, 2.0)]),
        vec![],
    );
    assert_eq!(
        polygon.check_validation().unwrap_err(),
        InvalidPolygon::SelfIntersection(RingRole::Exterior)
    );
}

#[test]
fn spike_is_self_intersecting() {
    let polygon = Polygon::new(
        ls(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (2.0, 6.0),
            (2.0, 4.0),
            (0.0, 4.0),
        ]),
        vec![],
    );
    assert_eq!(
        polygon.validation_errors(),
        vec![InvalidPolygon::SelfIntersection(RingRole::Exterior)]
    );
}

#[test]
fn too_few_points_in_ring() {
    let polygon = Polygon::new(ls(&[(0.0, 0.0), (1.0, 1.0)]), vec![]);
    assert_eq!(
        polygon.validation_errors(),
        vec![InvalidPolygon::TooFewPointsInRing(RingRole::Exterior)]
    );
}

#[test]
fn hole_outside_shell_not_contained() {
    let polygon = Polygon::new(
        ls(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
        vec![ls(&[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0)])],
    );
    assert_eq!(
        polygon.check_validation().unwrap_err(),
        InvalidPolygon::InteriorRingNotContainedInExteriorRing(RingRole::Interior(0))
    );
}

#[test]
fn holes_touching_at_point_are_valid_on_a_line_invalid() {
    let shell = ls(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
    let point_touch = Polygon::new(
        shell.clone(),
        vec![
            ls(&[(1.0, 2.0), (2.0, 1.0), (3.0, 2.0), (2.0, 3.0)]),
            ls(&[(3.0, 2.0), (4.0, 1.0), (5.0, 2.0), (4.0, 3.0)]),
        ],
    );
    assert!(point_touch.is_valid());

    let line_touch = Polygon::new(
        shell,
        vec![
            ls(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]),
            ls(&[(3.0, 1.0), (5.0, 1.0), (5.0, 3.0), (3.0, 3.0)]),
        ],
    );
    assert_eq!(
        line_touch.check_validation().unwrap_err(),
        InvalidPolygon::IntersectingRingsOnALine(RingRole::Interior(0), RingRole::Interior(1))
    );
}

#[test]
fn nested_holes_intersect_on_an_area() {
    let polygon = Polygon::new(
        ls(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        vec![
            ls(&[(1.0, 1.0), (1.0, 9.0), (9.0, 9.0), (9.0, 1.0)]),
            ls(&[(2.0, 2.0), (2.0, 8.0), (8.0, 8.0), (8.0, 2.0)]),
        ],
    );
    assert_eq!(
        polygon.check_validation().unwrap_err(),
        InvalidPolygon::IntersectingRingsOnAnArea(RingRole::Interior(0), RingRole::Interior(1))
    );
}

#[test]
fn nan_coordinate_is_reported_with_its_index() {
    let line = ls(&[(0.0, 0.0), (f64::NAN, 1.0), (2.0, 2.0)]);
    assert_eq!(
        line.validation_errors(),
        vec![InvalidLineString::NonFiniteCoord(CoordIndex(1))]
    );
}

#[test]
fn multipolygon_overlap_and_touch() {
    let a = Polygon::new(
        ls(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
        vec![],
    );
    let overlapping = Polygon::new(
        ls(&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]),
        vec![],
    );
    let mp = MultiPolygon::new(vec![a.clone(), overlapping]);
    assert_eq!(
        mp.check_validation().unwrap_err(),
        InvalidMultiPolygon::ElementsOverlaps(GeometryIndex(0), GeometryIndex(1))
    );

    let edge_sharing = Polygon::new(
        ls(&[(4.0, 0.0), (8.0, 0.0), (8.0, 4.0), (4.0, 4.0)]),
        vec![],
    );
    let mp = MultiPolygon::new(vec![a.clone(), edge_sharing]);
    assert_eq!(
        mp.check_validation().unwrap_err(),
        InvalidMultiPolygon::ElementsTouchOnALine(GeometryIndex(0), GeometryIndex(1))
    );

    let disjoint = Polygon::new(
        ls(&[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0)]),
        vec![],
    );
    assert!(MultiPolygon::new(vec![a, disjoint]).is_valid());
}

#[test]
fn geometry_collection_reports_member_index() {
    let g = Geometry::GeometryCollection(vec![
        Geometry::Point(Point::new(0.0, 0.0)),
        Geometry::LineString(ls(&[(0.0, 0.0), (0.0, 0.0)])),
    ]);
    assert_eq!(
        g.validation_errors(),
        vec![InvalidGeometry::InvalidGeometryCollection(
            GeometryIndex(1),
            Box::new(InvalidGeometry::InvalidLineString(
                InvalidLineString::TooFewPoints
            )),
        )]
    );
}
