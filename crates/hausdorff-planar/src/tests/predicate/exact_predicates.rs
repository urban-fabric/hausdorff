//! The exact integer predicates, pushed to the two places they can fail:
//! **near-degeneracy**, where the whole point is to beat `f64`, and the
//! **`i128` ceiling**, where the answer must be a reported overflow rather than a
//! wrapped sign.
//!
//! These predicates exist because a floating-point determinant can return the
//! *wrong sign* — not an imprecise magnitude, a wrong sign — when three points are
//! nearly collinear. A wrong sign is categorically worse than an imprecise number
//! because the sweep and the convex-hull walk branch on it: one bad sign and the
//! algorithm takes a structurally impossible path and produces nonsense or loops.
//! So the tests below are built around a sharp question: *is there an input where
//! `f64` gets the sign wrong and the exact predicate gets it right?* If the answer
//! is ever no, the exact path has stopped earning its cost.
//!
//! **The oracle is construction, not comparison.** Asserting that the exact
//! predicate agrees with an independently-computed `i128` determinant would be
//! circular — it is the same arithmetic. Instead every case is *built* from a known
//! answer: take a point exactly on the line `ab`, then move it one unit in `y`. The
//! cross product is then exactly `±(b.x − a.x)`, whose sign is known in closed form
//! before any predicate runs. One unit on an integer grid is the smallest possible
//! perturbation, so these are the hardest inputs that still have a defined answer.
//!
//! The `i128` ceiling is not one number, because the predicates have different
//! degrees. `orient` is a 2×2 determinant of differences, so it needs roughly
//! `2·log2(coordinate)` bits; `in_circle` is a 4×4 determinant involving squared
//! norms and needs roughly `4·log2(coordinate)`. A coordinate magnitude near `2^42`
//! is therefore comfortable for `orient` and already past the ceiling for
//! `in_circle`, and the tests record that asymmetry rather than picking one bound
//! and pretending it covers both.

use crate::space::r2::{
    ExactCoord, ExactPredicateError, InCirclePosition, OrientationSign, in_circle_position, orient,
    orientation_determinant, point_on_segment, segment_intersection, segments_properly_cross,
};

/// The floating-point determinant the exact path exists to correct. Written out
/// here rather than imported so the test does not need a dependency on
/// `hausdorff-kernels` just to demonstrate a rounding failure — this *is* the
/// formula, and the point is what `f64` does to it.
fn orient_f64(a: ExactCoord, b: ExactCoord, c: ExactCoord) -> f64 {
    let (ax, ay) = (a.x as f64, a.y as f64);
    let (bx, by) = (b.x as f64, b.y as f64);
    let (cx, cy) = (c.x as f64, c.y as f64);
    (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
}

#[test]
fn exact_orientation_beats_f64_on_near_collinear_triples() {
    // `c` sits exactly on the line through `a` and `b`, then is nudged one unit up.
    // The determinant is then exactly `b.x - a.x`, so the true orientation is known
    // without computing anything: positive `b.x - a.x` means counter-clockwise.
    //
    // The magnitudes are chosen above 2^53 so that `c.y + 1` is *not representable*
    // in f64 — the nudge vanishes on conversion and the float sees three collinear
    // points. That is the failure this predicate exists to prevent.
    let mut float_got_it_wrong = 0;
    for shift in [54u32, 56, 58, 60] {
        let m = 1i128 << shift;
        let a = ExactCoord::new(0, 0);
        let b = ExactCoord::new(m, m);
        // On the line y = x, nudged up by one.
        let c = ExactCoord::new(m / 2, m / 2 + 1);

        let exact = orient(a, b, c).expect("within the i128 range");
        assert_eq!(
            exact,
            OrientationSign::CounterClockwise,
            "2^{shift}: a one-unit nudge above the line is counter-clockwise, \
             and the determinant is exactly b.x - a.x = {m}"
        );
        assert_eq!(
            orientation_determinant(a, b, c).expect("in range"),
            m,
            "2^{shift}: the determinant has a known closed form"
        );

        // The same points below the line give the opposite sign, exactly.
        let c_below = ExactCoord::new(m / 2, m / 2 - 1);
        assert_eq!(
            orient(a, b, c_below).expect("in range"),
            OrientationSign::Clockwise,
            "2^{shift}: one unit below the line is clockwise"
        );

        if orient_f64(a, b, c) == 0.0 {
            float_got_it_wrong += 1;
        }
    }
    assert!(
        float_got_it_wrong > 0,
        "these cases are supposed to be ones f64 cannot resolve; if f64 now gets \
         them all right the test has stopped testing anything"
    );
}

#[test]
fn exactly_collinear_points_are_collinear_at_every_scale() {
    // Collinearity is the knife edge: it must be reported as collinear and not as a
    // tiny nonzero sign. Points on y = 2x, at scales spanning the usable range.
    for shift in 0..42u32 {
        let m = 1i128 << shift;
        let a = ExactCoord::new(0, 0);
        let b = ExactCoord::new(m, 2 * m);
        let c = ExactCoord::new(3 * m, 6 * m);
        assert_eq!(
            orient(a, b, c).expect("in range"),
            OrientationSign::Collinear,
            "2^{shift}: (0,0), (m,2m), (3m,6m) are exactly collinear"
        );
        assert_eq!(orientation_determinant(a, b, c).expect("in range"), 0);
    }
    // Degenerate triples: repeated and coincident points are collinear, not an error.
    let p = ExactCoord::new(7, -3);
    let q = ExactCoord::new(11, 5);
    for (a, b, c) in [(p, p, p), (p, p, q), (p, q, p), (p, q, q)] {
        assert_eq!(
            orient(a, b, c).expect("in range"),
            OrientationSign::Collinear,
            "a triple with a repeated point spans no area"
        );
    }
}

#[test]
fn orientation_is_antisymmetric_under_a_swap() {
    // Swapping two vertices must flip the sign, and a cyclic rotation must preserve
    // it. These are properties of the determinant, so they hold for every input,
    // which makes them a sweep rather than a case.
    let pts = [
        ExactCoord::new(0, 0),
        ExactCoord::new(5, 1),
        ExactCoord::new(2, 9),
        ExactCoord::new(-7, 3),
        ExactCoord::new(1 << 30, -(1 << 29)),
        ExactCoord::new(-(1 << 31), 1 << 20),
    ];
    for a in pts {
        for b in pts {
            for c in pts {
                let abc = orient(a, b, c).expect("in range");
                let bac = orient(b, a, c).expect("in range");
                let bca = orient(b, c, a).expect("in range");
                let flipped = match bac {
                    OrientationSign::Clockwise => OrientationSign::CounterClockwise,
                    OrientationSign::CounterClockwise => OrientationSign::Clockwise,
                    OrientationSign::Collinear => OrientationSign::Collinear,
                };
                assert_eq!(
                    abc, flipped,
                    "swapping the first two vertices must flip the sign"
                );
                assert_eq!(abc, bca, "a cyclic rotation must preserve the sign");
            }
        }
    }
}

#[test]
fn overflow_is_reported_not_wrapped() {
    // The guarantee the module's docs make is "overflow reported instead of
    // silently wrapping". A wrapped sign would be indistinguishable from a real
    // answer, which is the worst possible failure for a predicate whose entire
    // value is sign correctness.
    let big = i128::MAX;
    let extremes = [
        (
            ExactCoord::new(big, big),
            ExactCoord::new(big, -big),
            ExactCoord::new(-big, big),
        ),
        (
            ExactCoord::new(i128::MIN, 0),
            ExactCoord::new(big, 0),
            ExactCoord::new(0, big),
        ),
        (
            ExactCoord::new(1i128 << 64, 0),
            ExactCoord::new(0, 1i128 << 64),
            ExactCoord::new(-(1i128 << 64), 0),
        ),
        // Collinear but far out of range: still overflow, not a lucky zero.
        (
            ExactCoord::new(0, 0),
            ExactCoord::new(big, big),
            ExactCoord::new(-big, -big),
        ),
    ];
    for (a, b, c) in extremes {
        assert_eq!(
            orient(a, b, c),
            Err(ExactPredicateError::Overflow),
            "{a:?} {b:?} {c:?} exceeds i128 and must say so"
        );
        assert_eq!(
            orientation_determinant(a, b, c),
            Err(ExactPredicateError::Overflow)
        );
        assert_eq!(
            in_circle_position(a, b, c, ExactCoord::new(1, 1)),
            Err(ExactPredicateError::Overflow)
        );
    }
    // `i128::MIN` has no positive counterpart, so negating it is its own overflow
    // path; worth a case of its own because `checked_neg` is easy to omit.
    let m = ExactCoord::new(i128::MIN, i128::MIN);
    assert!(
        orient(m, ExactCoord::new(0, 0), ExactCoord::new(1, 1)).is_err(),
        "i128::MIN must not wrap through negation"
    );
}

#[test]
fn the_two_predicates_have_different_usable_ranges() {
    // `orient` is degree 2 and `in_circle` degree 4, so there is a band of
    // coordinate magnitudes where one works and the other cannot. Pinned so the
    // asymmetry is a documented property rather than a surprise at a call site.
    let m = 1i128 << 42;
    let (a, b, c) = (
        ExactCoord::new(m, 0),
        ExactCoord::new(0, m),
        ExactCoord::new(-m, 0),
    );
    assert!(
        orient(a, b, c).is_ok(),
        "a degree-2 determinant is comfortable at 2^42"
    );
    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(1, 1)),
        Err(ExactPredicateError::Overflow),
        "a degree-4 determinant is already past the ceiling at 2^42"
    );
    // At a modest scale both work.
    let s = 1i128 << 20;
    let (a, b, c) = (
        ExactCoord::new(s, 0),
        ExactCoord::new(0, s),
        ExactCoord::new(-s, 0),
    );
    assert!(orient(a, b, c).is_ok());
    assert!(in_circle_position(a, b, c, ExactCoord::new(0, 0)).is_ok());
}

#[test]
fn in_circle_agrees_with_a_known_circumcircle() {
    // The unit circle scaled by r, with three points on it and a fourth placed
    // inside, on, and outside by construction — so the expected answer is known
    // without computing a determinant.
    let r = 1_000_000i128;
    let a = ExactCoord::new(r, 0);
    let b = ExactCoord::new(0, r);
    let c = ExactCoord::new(-r, 0);
    // (0, -r) is the fourth point of the same circle.
    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(0, -r)).expect("in range"),
        InCirclePosition::OnCircle,
        "the antipode of b is on the same circle"
    );
    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(0, 0)).expect("in range"),
        InCirclePosition::Inside,
        "the centre is inside"
    );
    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(2 * r, 0)).expect("in range"),
        InCirclePosition::Outside,
        "twice the radius is outside"
    );
    // One unit inside and one unit outside — the near-degenerate case, where the
    // floating-point version would be at its least trustworthy.
    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(r - 1, 0)).expect("in range"),
        InCirclePosition::Inside,
        "one unit inside the rim is inside"
    );
    assert_eq!(
        in_circle_position(a, b, c, ExactCoord::new(r + 1, 0)).expect("in range"),
        InCirclePosition::Outside,
        "one unit outside the rim is outside"
    );
    // A collinear defining triple has no circumcircle, and that is a distinct error
    // from overflow.
    assert_eq!(
        in_circle_position(
            ExactCoord::new(0, 0),
            ExactCoord::new(1, 1),
            ExactCoord::new(2, 2),
            ExactCoord::new(5, 0)
        ),
        Err(ExactPredicateError::DegenerateCircumcircle),
        "collinear points define no circle; this must not be reported as overflow"
    );
}

#[test]
fn segment_predicates_handle_touching_and_collinear_overlap() {
    let o = ExactCoord::new(0, 0);
    let e = ExactCoord::new(10, 10);
    // A proper crossing.
    assert_eq!(
        segments_properly_cross(o, e, ExactCoord::new(0, 10), ExactCoord::new(10, 0)),
        Ok(true)
    );
    // Touching at an endpoint is *not* a proper crossing — the distinction the
    // sweep depends on.
    assert_eq!(
        segments_properly_cross(o, e, e, ExactCoord::new(20, 0)),
        Ok(false),
        "sharing an endpoint is a touch, not a proper crossing"
    );
    // T-junction: an endpoint lying in the other segment's interior.
    assert_eq!(
        segments_properly_cross(o, e, ExactCoord::new(5, 5), ExactCoord::new(10, 0)),
        Ok(false),
        "an endpoint on the interior of the other segment is not a proper crossing"
    );
    // Collinear overlap, containment, and disjointness along one line.
    for (c, d, label) in [
        (
            ExactCoord::new(5, 5),
            ExactCoord::new(15, 15),
            "partial overlap",
        ),
        (ExactCoord::new(2, 2), ExactCoord::new(8, 8), "contained"),
        (ExactCoord::new(20, 20), ExactCoord::new(30, 30), "disjoint"),
        (
            ExactCoord::new(10, 10),
            ExactCoord::new(20, 20),
            "endpoint touch",
        ),
    ] {
        let got = segment_intersection(o, e, c, d);
        assert!(
            got.is_ok(),
            "{label}: should classify, not error, got {got:?}"
        );
        assert_eq!(
            segments_properly_cross(o, e, c, d),
            Ok(false),
            "{label}: collinear segments never *properly* cross"
        );
    }
    // Zero-length segments: a degenerate "segment" is a point.
    assert_eq!(
        point_on_segment(o, o, o),
        Ok(true),
        "a point lies on the degenerate segment at itself"
    );
    assert_eq!(
        point_on_segment(o, o, ExactCoord::new(1, 0)),
        Ok(false),
        "a different point does not"
    );
    assert_eq!(
        segments_properly_cross(o, o, o, o),
        Ok(false),
        "two degenerate segments do not properly cross"
    );
    // Endpoint containment, inclusive.
    assert_eq!(point_on_segment(o, e, o), Ok(true));
    assert_eq!(point_on_segment(o, e, e), Ok(true));
    assert_eq!(point_on_segment(o, e, ExactCoord::new(5, 5)), Ok(true));
    assert_eq!(point_on_segment(o, e, ExactCoord::new(5, 6)), Ok(false));
    assert_eq!(
        point_on_segment(o, e, ExactCoord::new(11, 11)),
        Ok(false),
        "collinear but beyond the endpoint is not on the segment"
    );
}
