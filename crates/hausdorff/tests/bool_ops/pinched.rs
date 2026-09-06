//! Pinch geometry: shapes whose boundary touches itself at a point.
//!
//! The suite's other fixtures are convex or grid-aligned, so nothing here was
//! covered before — which is how a union could return less area than one of its
//! operands without a test noticing (#44).

use super::*;

/// Load a WKT fixture from `tests/fixtures`. The large pinch cases live on disk
/// rather than inline because their coordinates must survive verbatim.
fn fixture(name: &str) -> MultiPolygon {
    let path = format!("{}/tests/fixtures/{name}.wkt", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    match hausdorff::format::wkt::from_wkt(text.trim()).expect("fixture parses") {
        hausdorff::Geometry::MultiPolygon(mp) => mp,
        other => panic!("expected a multipolygon, got {other:?}"),
    }
}

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

// #44 was two defects, and the reproducers in this file cover both.
//
// *Contour tracing* closed a ring with a chord that was never in the result,
// because at a point where four or more result edges met it continued along
// whichever one happened to sit next in the sorted event list. Tracing now
// orients every result edge so the interior is on its left and pairs them at a
// shared point by angle, so the continuation follows from the labelling.
//
// *Region labelling* then left the tracer edge sets it could not close at all —
// a point with an **odd** number of result edges meeting, which no closed
// boundary has. The union below is the cheapest reproducer of that in the
// suite: twelve result segments, two of them degree-1. The triangle's first
// vertex lies on a hexagon edge to within two ULPs, so the sweep splits the edge
// there while ordering the segments as though it had not, and labels a triangle
// edge inside the hexagon when it is outside. Dropping it opens the corner.
#[test]
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

/// A second #44 reproducer, from the same isochrone-band run as the fixtures in
/// [`super::isochrone_bands`] but a different pair. It failed *differently* from
/// the hexagon/triangle above — that one is three vertices against six, this is
/// 574 against 393 with many self-touching pinches — so it is kept as the check
/// that a repair works at scale and not just on the small case. Stored as WKT
/// because the coordinates must survive verbatim; the conversion from the
/// captured WKB was checked to re-encode byte for byte.
fn pinch_band() -> (MultiPolygon, MultiPolygon) {
    (fixture("pinch_band_a"), fixture("pinch_band_b"))
}

#[test]
fn union_of_a_pinched_band_pair_covers_both() {
    let (a, b) = pinch_band();
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
fn difference_of_a_pinched_band_pair_is_consistent() {
    let (a, b) = pinch_band();
    let inter = ha(&a.intersection(&b));
    let diff = ha(&a.difference(&b));
    assert!(diff <= ha(&a) + EPS, "A\\B ≤ A: {diff} vs {}", ha(&a));
    assert!(
        (diff - (ha(&a) - inter)).abs() < EPS,
        "A\\B = A-A∩B: {diff} vs {}",
        ha(&a) - inter
    );
}

/// The band pair that held out longest. Its operands share exact vertices in two
/// places and cross again a millionth of a degree further along, which left the
/// crossing pair separated by a third segment in the sweep's status — never
/// adjacent, so never tested, so never split. An edge of one polygon then ran
/// from inside the other to outside it with nothing marking where it changed,
/// and the tracer was handed a vertex with a single edge at it.
///
/// The whole set-algebra battery, because the failure showed in every operation.
/// The intersection area is a golden constant in the sense of this file's header
/// — captured once the invariants below all held, and kept so a regression that
/// preserved the identities while moving the answer would still be caught.
#[test]
fn shared_vertex_band_pair_satisfies_the_set_algebra() {
    let (a, b) = (
        fixture("shared_vertex_band_a"),
        fixture("shared_vertex_band_b"),
    );
    check_ops(&a, &b, 0.013_051_614_871_983_47);
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
