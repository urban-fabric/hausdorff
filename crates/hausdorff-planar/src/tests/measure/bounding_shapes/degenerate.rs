//! Degenerate input across all three descriptors at once: empty, single-point,
//! collinear and repeated-coordinate geometries must return `None` or a sane
//! answer rather than panicking.

use super::*;

#[test]
fn degenerate_input_is_handled_rather_than_panicking() {
    let empty = MultiPoint::new(vec![]);
    assert!(
        empty.minimum_bounding_circle().is_none(),
        "no points, no circle"
    );
    assert!(empty.oriented_envelope().is_none());

    // A single point: a circle of radius zero, and a degenerate envelope.
    let one = MultiPoint::new(vec![Point(Coord::new(3.0, 4.0))]);
    let c = one
        .minimum_bounding_circle()
        .expect("one point still has a circle");
    close(c.radius, 0.0, 1e-12, "a single point has radius zero");
    close(c.centre.x, 3.0, 1e-12, "centred on the point");
    let env = one.oriented_envelope().expect("degenerate envelope");
    close(oriented_envelope_area(env), 0.0, 1e-12, "zero area");

    // Collinear points: the circumcircle branch has no solution, so the two-point
    // branch has to carry it.
    let line = MultiPoint::new(vec![
        Point(Coord::new(0.0, 0.0)),
        Point(Coord::new(1.0, 1.0)),
        Point(Coord::new(2.0, 2.0)),
        Point(Coord::new(3.0, 3.0)),
    ]);
    let c = line
        .minimum_bounding_circle()
        .expect("collinear still bounds");
    close(c.radius, (18.0_f64).sqrt() / 2.0, 1e-9, "half the span");
    for p in line.to_coords() {
        assert!(c.contains(p), "{p:?} outside the collinear bounding circle");
    }
    let env = line.oriented_envelope().expect("collinear envelope");
    close(env.width, 0.0, 1e-9, "a collinear set has no width");

    // A geometry with no interior has no inscribed circle.
    assert!(
        LineString::from(vec![(0.0, 0.0), (1.0, 1.0)])
            .maximum_inscribed_circle(1e-3)
            .is_none(),
        "a linestring has no interior to inscribe into"
    );

    // Non-finite coordinates are dropped rather than poisoning the result.
    let nasty = MultiPoint::new(vec![
        Point(Coord::new(0.0, 0.0)),
        Point(Coord::new(f64::NAN, 1.0)),
        Point(Coord::new(2.0, 0.0)),
        Point(Coord::new(f64::INFINITY, 0.0)),
    ]);
    let c = nasty
        .minimum_bounding_circle()
        .expect("the finite points still bound");
    assert!(
        c.radius.is_finite() && c.centre.is_finite(),
        "non-finite input must not poison the circle, got {c:?}"
    );
    close(c.radius, 1.0, 1e-9, "the two finite points are 2 apart");
}
