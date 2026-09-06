//! Pinch geometry: shapes whose boundary touches itself at a point.
//!
//! The suite's other fixtures are convex or grid-aligned, so nothing here was
//! covered before — which is how a union could return less area than one of its
//! operands without a test noticing (#44).

use super::*;

/// A hexagon and a triangle captured from an isochrone band, translated to the
/// origin. Verbatim, because the defect they provoke is a topological one and
/// rounding these coordinates loses it.
fn hexagon_and_triangle() -> (MultiPolygon, MultiPolygon) {
    let a = hmp(&[&[&[
        (0.07580900000002089, -0.1411330000000035),
        (-0.055163999999976454, 0.08880199999998695),
        (-0.05946699999998373, 0.08293199999999601),
        (-0.060502999999982876, 0.08360000000000056),
        (-0.0852938073957148, 0.11397130958639679),
        (-0.10349099999997691, 0.11004299999998679),
    ]]]);
    let b = hmp(&[&[&[
        (-0.060116323943645966, 0.08335067605633739),
        (-0.04725199999998608, 0.08127499999998733),
        (-0.046615999999986, 0.07875199999998728),
    ]]]);
    (a, b)
}

// Both sliver cases fail today. They are the reproducer on #44, kept runnable
// and named so the gap is visible rather than absent: `cargo test -- --ignored`
// runs them, and the ignore comes off with the fix.
#[test]
#[ignore = "reproduces #44: contour tracing closes rings with an invented chord"]
fn union_of_a_sliver_and_a_triangle_covers_both() {
    let (a, b) = hexagon_and_triangle();
    let (aa, ba) = (ha(&a), ha(&b));
    let union = ha(&a.union(&b));
    assert!(
        union >= aa.max(ba) - EPS,
        "A∪B must cover A: union {union} vs A {aa}"
    );
    let inter = ha(&a.intersection(&b));
    assert!(
        (union - (aa + ba - inter)).abs() < EPS,
        "A∪B = A+B-A∩B: {union} vs {}",
        aa + ba - inter
    );
}

#[test]
#[ignore = "reproduces #44: contour tracing closes rings with an invented chord"]
fn difference_of_a_sliver_and_a_triangle_is_consistent() {
    let (a, b) = hexagon_and_triangle();
    let inter = ha(&a.intersection(&b));
    let diff = ha(&a.difference(&b));
    assert!(
        (diff - (ha(&a) - inter)).abs() < EPS,
        "A\\B = A-A∩B: {diff} vs {}",
        ha(&a) - inter
    );
}

/// Two unit squares meeting at exactly one vertex. The union has to keep both.
#[test]
fn squares_touching_at_one_vertex_keep_both_lobes() {
    let a = hmp(&[&[&sq(0.0, 0.0, 1.0, 1.0)]]);
    let b = hmp(&[&[&sq(-1.0, -1.0, 0.0, 0.0)]]);
    check_ops(&a, &b, 0.0);
    assert!((ha(&a.union(&b)) - 2.0).abs() < EPS);
}

/// A square with a second square hanging off one corner, unioned with a box
/// that overlaps only the upper lobe — the pinch survives into the result.
#[test]
fn pinch_survives_into_the_result() {
    let bowtie_ish = hmp(&[&[&sq(0.0, 0.0, 2.0, 2.0)], &[&sq(-2.0, -2.0, 0.0, 0.0)]]);
    let overlap = hmp(&[&[&sq(1.0, 1.0, 3.0, 3.0)]]);
    // Upper lobe 4 + lower lobe 4 + the part of `overlap` outside them (4 - 1).
    check_ops(&bowtie_ish, &overlap, 1.0);
    assert!((ha(&bowtie_ish.union(&overlap)) - 11.0).abs() < EPS);
}
