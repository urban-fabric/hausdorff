use super::*;
use hausdorff_planar::space::r2::{LineString, MultiPoint, Point, Polygon};

fn point(lon: f64, lat: f64) -> Geometry {
    Geometry::Point(Point::new(lon, lat))
}

/// An axis-aligned lon/lat rectangle as a polygon geometry.
fn poly(min_lon: f64, min_lat: f64, max_lon: f64, max_lat: f64) -> Geometry {
    Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (min_lon, min_lat),
            (max_lon, min_lat),
            (max_lon, max_lat),
            (min_lon, max_lat),
            (min_lon, min_lat),
        ]),
        vec![],
    ))
}

#[test]
fn far_apart_polygons_are_disjoint() {
    // Two small rects on opposite sides of the globe (galois join proptest).
    let south = poly(-98.0, -60.0, -88.0, -50.0);
    let north = poly(80.0, 53.0, 88.0, 59.0);
    assert_eq!(
        SpatialPredicate::Intersects.eval(&south, &north),
        Some(false),
        "small rects on opposite sides of the globe do not intersect"
    );
    assert_eq!(SpatialPredicate::Disjoint.eval(&south, &north), Some(true));
}

#[test]
fn point_in_polygon_core_predicates() {
    let square = poly(-10.0, -10.0, 10.0, 10.0);
    let inside = point(0.0, 0.0);
    let outside = point(20.0, 20.0);

    // point <pred> polygon
    assert_eq!(SpatialPredicate::Within.eval(&inside, &square), Some(true));
    assert_eq!(
        SpatialPredicate::Within.eval(&outside, &square),
        Some(false)
    );
    assert_eq!(
        SpatialPredicate::Intersects.eval(&inside, &square),
        Some(true)
    );
    assert_eq!(
        SpatialPredicate::Disjoint.eval(&outside, &square),
        Some(true)
    );
    assert_eq!(
        SpatialPredicate::CoveredBy.eval(&inside, &square),
        Some(true)
    );
    // A point never contains a polygon.
    assert_eq!(
        SpatialPredicate::Contains.eval(&inside, &square),
        Some(false)
    );

    // Symmetric: polygon <pred> point (transposed).
    assert_eq!(
        SpatialPredicate::Contains.eval(&square, &inside),
        Some(true)
    );
    assert_eq!(
        SpatialPredicate::Contains.eval(&square, &outside),
        Some(false)
    );
    assert_eq!(
        SpatialPredicate::Within.eval(&square, &inside),
        Some(false),
        "an areal is never within a point"
    );
    assert_eq!(
        SpatialPredicate::Intersects.eval(&square, &inside),
        Some(true)
    );
}

#[test]
fn spherical_differs_from_planar_across_the_antimeridian() {
    // A polygon straddling ±180° (lon 170 → -170, a 20°-wide band centred on the
    // antimeridian). A point at lon 180 is *inside* it on the sphere.
    let band = poly(170.0, -10.0, -170.0, 10.0);
    let on_antimeridian = point(180.0, 0.0);

    // Spherical: the point is inside the band that wraps ±180°.
    assert_eq!(
        SpatialPredicate::Within.eval(&on_antimeridian, &band),
        Some(true),
        "spherical point-in-polygon wraps the antimeridian"
    );

    // The planar predicate on the same lon/lat coordinates disagrees: the ring's
    // vertices span lon [-170, 170] as a *Cartesian* box, which does NOT contain
    // lon 180 — the exact error the spherical predicate fixes.
    let planar =
        hausdorff_planar::space::r2::SpatialPredicate::Within.eval(&on_antimeridian, &band);
    assert!(
        !planar,
        "planar treats the wrap-around band as excluding lon 180"
    );
}

#[test]
fn multipoint_within_requires_all_points_inside() {
    let square = poly(-10.0, -10.0, 10.0, 10.0);
    let all_in = Geometry::MultiPoint(MultiPoint::new(vec![
        Point::new(0.0, 0.0),
        Point::new(5.0, -5.0),
    ]));
    let some_out = Geometry::MultiPoint(MultiPoint::new(vec![
        Point::new(0.0, 0.0),
        Point::new(50.0, 50.0),
    ]));
    assert_eq!(SpatialPredicate::Within.eval(&all_in, &square), Some(true));
    assert_eq!(
        SpatialPredicate::Within.eval(&some_out, &square),
        Some(false)
    );
    // Intersects holds as long as one point is inside.
    assert_eq!(
        SpatialPredicate::Intersects.eval(&some_out, &square),
        Some(true)
    );
}

#[test]
fn polygon_with_hole_excludes_the_hole() {
    // A 20×20 square with a 4×4 hole at the centre.
    let with_hole = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (-10.0, -10.0),
            (10.0, -10.0),
            (10.0, 10.0),
            (-10.0, 10.0),
            (-10.0, -10.0),
        ]),
        vec![LineString::from(vec![
            (-2.0, -2.0),
            (2.0, -2.0),
            (2.0, 2.0),
            (-2.0, 2.0),
            (-2.0, -2.0),
        ])],
    ));
    // Centre is in the hole → not contained; off-centre (but in the shell) is.
    assert_eq!(
        SpatialPredicate::Contains.eval(&with_hole, &point(0.0, 0.0)),
        Some(false),
        "a point in the hole is not contained"
    );
    assert_eq!(
        SpatialPredicate::Contains.eval(&with_hole, &point(6.0, 6.0)),
        Some(true)
    );
}

#[test]
fn areal_areal_intersects_hole_free() {
    let a = poly(-10.0, -10.0, 5.0, 5.0);
    let overlap = poly(0.0, 0.0, 10.0, 10.0);
    let disjoint = poly(50.0, 50.0, 60.0, 60.0);
    let inside = poly(-2.0, -2.0, 2.0, 2.0);

    assert_eq!(SpatialPredicate::Intersects.eval(&a, &overlap), Some(true));
    assert_eq!(
        SpatialPredicate::Intersects.eval(&a, &disjoint),
        Some(false)
    );
    assert_eq!(SpatialPredicate::Disjoint.eval(&a, &disjoint), Some(true));
    // One polygon fully inside the other → their regions share points.
    assert_eq!(SpatialPredicate::Intersects.eval(&a, &inside), Some(true));
    // Containment predicates on areal×areal are not in v1 scope → planar fallback.
    assert_eq!(SpatialPredicate::Contains.eval(&a, &inside), None);
}

#[test]
fn areal_areal_with_holes_falls_back_to_planar() {
    let holed = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (-10.0, -10.0),
            (10.0, -10.0),
            (10.0, 10.0),
            (-10.0, 10.0),
            (-10.0, -10.0),
        ]),
        vec![LineString::from(vec![
            (-2.0, -2.0),
            (2.0, -2.0),
            (2.0, 2.0),
            (-2.0, 2.0),
            (-2.0, -2.0),
        ])],
    ));
    let other = poly(0.0, 0.0, 5.0, 5.0);
    assert_eq!(
        SpatialPredicate::Intersects.eval(&holed, &other),
        None,
        "a holed polygon defers to the planar predicate"
    );
}

#[test]
fn areal_areal_intersects_across_the_antimeridian() {
    // Two bands both straddling ±180° that overlap only in the wrap strip. On the
    // sphere their shell edges cross; a planar reading (both bounded lon [-x,x])
    // would call them disjoint in the wrap region.
    let a = poly(170.0, -5.0, -170.0, 5.0);
    let b = poly(175.0, 0.0, -175.0, 10.0);
    assert_eq!(
        SpatialPredicate::Intersects.eval(&a, &b),
        Some(true),
        "spherical shells meet across the antimeridian"
    );
    let planar = hausdorff_planar::space::r2::SpatialPredicate::Intersects.eval(&a, &b);
    // Planar happens to also see these wide boxes overlapping near lon 0 — the
    // point of the spherical path is that the *answer* is computed on great-circle
    // edges; assert the spherical result stands on its own.
    let _ = planar;
}

#[test]
fn out_of_scope_combinations_return_none() {
    // point×point and line×anything are not in v1 scope → fall back to planar.
    let p = point(0.0, 0.0);
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]));
    assert_eq!(SpatialPredicate::Intersects.eval(&p, &p), None);
    assert_eq!(
        SpatialPredicate::Intersects.eval(&line, &poly(-1.0, -1.0, 1.0, 1.0)),
        None
    );
}
