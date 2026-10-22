//! Buffer properties, stated as areas and containments rather than as
//! coordinate lists.
//!
//! A buffer is an approximation — the arcs are polygonal — so asserting exact
//! coordinates would pin the arc refinement rather than the operation. What is
//! worth asserting is what the operation means: the area a known shape's
//! buffer must have, that refining the arcs converges on it from below, that
//! dilation only ever adds and erosion only ever removes, and that a point's
//! buffer contains exactly the points within the distance.

use super::*;
use crate::space::r2::{Area, Contains, LineString, MultiPoint, Point, Polygon};

const PI: f64 = std::f64::consts::PI;

fn unit_square() -> Geometry {
    Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.0, 0.0),
        ]),
        Vec::new(),
    ))
}

fn area(polygons: &MultiPolygon) -> f64 {
    polygons.unsigned_area()
}

#[test]
fn a_points_buffer_approaches_the_disc_area_from_below() {
    let point = Geometry::Point(Point::new(0.0, 0.0));
    let exact = PI * 4.0;
    let mut previous = 0.0;
    for segments in [2u32, 4, 8, 16, 64, 256] {
        let buffered = buffer(&point, 2.0, BufferStyle::with_quadrant_segments(segments));
        let got = area(&buffered);
        // An inscribed polygon is always smaller than its circle, and refining
        // only adds area.
        assert!(
            got < exact,
            "{segments} segments gave {got}, which exceeds the disc area {exact}"
        );
        assert!(
            got > previous,
            "{segments} segments gave {got}, not more than {previous}"
        );
        previous = got;
    }
    // At the default refinement the error should already be well under a
    // percent; this is the claim the module documentation makes.
    let default = area(&buffer(&point, 2.0, BufferStyle::default()));
    assert!(
        (exact - default) / exact < 0.01,
        "default refinement is {default} against {exact}"
    );
}

#[test]
fn a_squares_buffer_has_the_area_of_its_parts() {
    // Dilating a w x h rectangle by d gives the rectangle, four d-wide
    // collars, and four quarter-discs that together make one disc:
    //   w*h + 2d(w + h) + pi d^2
    let (w, h, d) = (1.0, 1.0, 0.5);
    let exact = w * h + 2.0 * d * (w + h) + PI * d * d;
    let got = area(&buffer(
        &unit_square(),
        d,
        BufferStyle::with_quadrant_segments(256),
    ));
    assert!(
        (got - exact).abs() / exact < 1e-3,
        "buffered square area {got} against {exact}"
    );
}

#[test]
fn a_segments_buffer_has_the_area_of_a_stadium() {
    // A length-L segment dilated by d is a rectangle plus two half-discs:
    //   2dL + pi d^2
    let (length, d) = (3.0, 0.75);
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (length, 0.0)]));
    let exact = 2.0 * d * length + PI * d * d;
    let got = area(&buffer(&line, d, BufferStyle::with_quadrant_segments(256)));
    assert!(
        (got - exact).abs() / exact < 1e-3,
        "buffered segment area {got} against {exact}"
    );
}

#[test]
fn dilation_contains_the_original_and_erosion_is_contained_by_it() {
    let square = unit_square();
    let style = BufferStyle::default();
    let grown = buffer(&square, 0.25, style);
    let shrunk = buffer(&square, -0.25, style);
    let original = match &square {
        Geometry::Polygon(polygon) => polygon.clone(),
        _ => unreachable!("fixture is a polygon"),
    };

    // Area ordering is the cheap statement of the containment chain.
    assert!(area(&shrunk) < original.unsigned_area());
    assert!(original.unsigned_area() < area(&grown));

    // And the eroded square is the concentric one, so its area is known
    // exactly: side 1 - 2*0.25 = 0.5.
    assert!(
        (area(&shrunk) - 0.25).abs() < 1e-6,
        "eroded square area {}",
        area(&shrunk)
    );

    // Every original corner survives dilation; the centre survives erosion.
    for corner in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        assert!(
            grown.contains(&Point::new(corner.0, corner.1)),
            "dilation dropped the corner {corner:?}"
        );
    }
    assert!(shrunk.contains(&Point::new(0.5, 0.5)));
}

#[test]
fn eroding_past_the_inradius_empties_the_region() {
    // The unit square's inradius is 0.5, so retreating further than that
    // leaves nothing. An implementation that returned the input unchanged, or
    // a mirrored region, would fail here.
    let shrunk = buffer(&unit_square(), -0.6, BufferStyle::default());
    assert!(
        area(&shrunk) < 1e-9,
        "eroding past the inradius left area {}",
        area(&shrunk)
    );
}

#[test]
fn erosion_of_a_curve_or_a_point_is_empty() {
    // Neither has an interior to retreat into. Empty is the answer; returning
    // the input would claim a region that is not there.
    let style = BufferStyle::default();
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]));
    let point = Geometry::Point(Point::new(0.0, 0.0));
    assert!(buffer(&line, -1.0, style).0.is_empty());
    assert!(buffer(&point, -1.0, style).0.is_empty());
}

#[test]
fn zero_distance_returns_the_areal_part() {
    let style = BufferStyle::default();
    // Areal input keeps its area; a curve has none to keep.
    let square = buffer(&unit_square(), 0.0, style);
    assert!((area(&square) - 1.0).abs() < 1e-12);
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]));
    assert!(buffer(&line, 0.0, style).0.is_empty());
}

#[test]
fn overlapping_parts_are_dissolved_rather_than_double_counted() {
    // Two points closer together than twice the distance have overlapping
    // discs. The union must charge the overlap once, so the area is strictly
    // less than two discs.
    let points = Geometry::MultiPoint(MultiPoint::new(vec![
        Point::new(0.0, 0.0),
        Point::new(1.0, 0.0),
    ]));
    let style = BufferStyle::with_quadrant_segments(256);
    let joined = area(&buffer(&points, 1.0, style));
    let two_discs = 2.0 * PI;
    assert!(
        joined < two_discs,
        "overlap was double counted: {joined} against {two_discs}"
    );
    // And the result is one connected region rather than two.
    assert_eq!(buffer(&points, 1.0, style).0.len(), 1);

    // Far apart, the discs stay separate and the area is additive.
    let apart = Geometry::MultiPoint(MultiPoint::new(vec![
        Point::new(0.0, 0.0),
        Point::new(100.0, 0.0),
    ]));
    let separate = buffer(&apart, 1.0, style);
    assert_eq!(separate.0.len(), 2);
    assert!((area(&separate) - two_discs).abs() / two_discs < 1e-3);
}

#[test]
fn a_hole_is_preserved_and_closes_under_dilation() {
    // A square annulus: dilating from both boundaries narrows the hole, and a
    // large enough distance closes it entirely. An implementation that only
    // buffered the exterior ring would never close the hole.
    let outer = LineString::from(vec![
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ]);
    let hole = LineString::from(vec![
        (4.0, 4.0),
        (6.0, 4.0),
        (6.0, 6.0),
        (4.0, 6.0),
        (4.0, 4.0),
    ]);
    let annulus = Geometry::Polygon(Polygon::new(outer, vec![hole]));
    let style = BufferStyle::with_quadrant_segments(32);

    // A small dilation keeps the hole, smaller than it was.
    let small = buffer(&annulus, 0.25, style);
    assert_eq!(small.0.len(), 1);
    assert_eq!(
        small.0[0].interiors().len(),
        1,
        "a small dilation should keep the hole"
    );

    // A dilation past the hole's inradius closes it.
    let large = buffer(&annulus, 1.5, style);
    assert_eq!(large.0.len(), 1);
    assert!(
        large.0[0].interiors().is_empty(),
        "a dilation past the hole's inradius should close it"
    );
}

#[test]
fn a_non_finite_distance_describes_no_region() {
    let style = BufferStyle::default();
    for distance in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            buffer(&unit_square(), distance, style).0.is_empty(),
            "{distance} should produce no region"
        );
    }
}

#[test]
fn the_arc_refinement_is_clamped_rather_than_trusted() {
    // Accuracy, not meaning, so every value is interpretable — but one
    // segment per quadrant would turn a round join into a mitred one, and an
    // unbounded value is unbounded work driven by an argument.
    assert_eq!(
        BufferStyle::with_quadrant_segments(0).quadrant_segments(),
        2
    );
    assert_eq!(
        BufferStyle::with_quadrant_segments(1).quadrant_segments(),
        2
    );
    assert_eq!(
        BufferStyle::with_quadrant_segments(u32::MAX).quadrant_segments(),
        256
    );
    assert_eq!(BufferStyle::default().quadrant_segments(), 8);
}

#[test]
fn a_degenerate_segment_still_buffers_to_its_endpoint_disc() {
    // A zero-length segment contributes no rectangle; the endpoint discs are
    // what keep the result non-empty.
    let repeated = Geometry::LineString(LineString::from(vec![(1.0, 1.0), (1.0, 1.0)]));
    let style = BufferStyle::with_quadrant_segments(256);
    let buffered = buffer(&repeated, 1.0, style);
    assert_eq!(buffered.0.len(), 1);
    assert!((area(&buffered) - PI).abs() / PI < 1e-3);
}

#[test]
fn every_refinement_level_unions_without_a_degenerate_meeting() {
    // The discs are rotated by half a step precisely so no vertex lands on an
    // offset rectangle's edge. Without that, a 3-unit segment buffered by 0.75
    // aborts the union sweep at 4 and at 256 segments per quadrant — and
    // succeeds at 8 through 128 — because the sweep cannot label a boundary
    // meeting that touches instead of crossing.
    //
    // Swept across refinements rather than asserted at the default, since the
    // failure is a property of which coordinates coincide and the default was
    // one of the values that happened to work.
    let line = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (3.0, 0.0)]));
    let exact = 2.0 * 0.75 * 3.0 + PI * 0.75 * 0.75;
    let mut previous = 0.0;
    for segments in [2u32, 4, 8, 16, 32, 64, 128, 256] {
        let buffered = buffer(&line, 0.75, BufferStyle::with_quadrant_segments(segments));
        assert_eq!(
            buffered.0.len(),
            1,
            "{segments} segments gave {} parts; a segment's buffer is one region",
            buffered.0.len()
        );
        let got = area(&buffered);
        assert!(
            got < exact && got > previous,
            "{segments} segments gave {got}, which is not between {previous} and {exact}"
        );
        previous = got;
    }
}
