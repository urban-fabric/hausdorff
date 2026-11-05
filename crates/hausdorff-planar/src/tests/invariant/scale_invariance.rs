//! Two symmetries planar geometry must respect, and the constants that broke them.
//!
//! Translating a shape cannot change its area, and scaling it by `s` must scale
//! every area by `s²`. Both are exact statements about the plane, so any code that
//! violates them has an absolute number in it that should have been relative to
//! something. Each bug this suite covers was exactly that, and they came in two
//! flavours.
//!
//! **Position-scaled where it should be scale-free.** The shoelace sum
//! `Σ xᵢyᵢ₊₁ − xᵢ₊₁yᵢ` is translation-invariant in exact arithmetic and badly
//! conditioned in `f64`: it builds an answer the size of the *shape* out of terms
//! the size of the shape's *position squared*. Once `|coord|² > 2⁵³ ≈ 9.0e15` the
//! products stop being exact integers and a unit contribution falls below their
//! ULP, putting the cliff at a coordinate magnitude near `1e8`. A unit square
//! reported area `0` at an offset of `1e8`, `16384` at `1e10` and `134217728` at
//! `1e12`; `is_ccw` called it clockwise at all three; its centroid came back 3.3e9
//! away from it; and `polygonize` returned *no polygons at all* for a 3×3 grid of
//! unit cells at `1e9`, reporting all 27 edges as cut edges, because face
//! selection reads the sign of that sum. Separately, `relate` grew a segment's
//! bounding box by `1e-9 · |coord|` — 40 units at `1e10` — so two triangles
//! meeting at a point were reported as touching along a line, which made a valid
//! `MultiPolygon` invalid and `make_valid` fail on geometry it had just repaired
//! correctly.
//!
//! **Absolute where it should be relative.** `relate`'s collinearity test floored
//! its tolerance at an area of `1.0`, so for features below roughly `3e-5` units
//! the floor dominated every `orient` and all three points of any triangle read as
//! collinear. A metre is about `1e-5` of a degree, so that is ordinary data too.
//!
//! None of this needs a reference implementation to test: the area of a unit
//! square is 1 wherever it sits, and the two symmetries supply the rest.
//!
//! That the same mistake appeared in the shared shoelace kernel *and* four private
//! copies of it is why these tests are grouped by the hazard rather than filed
//! next to each function — a per-function test would have been written beside only
//! one of them.

use crate::space::r2::validation::Validation;
use crate::space::r2::{
    Area, BooleanOps, Centroid, Coord, LineString, MakeValid, MultiPolygon, Polygon, Polygonize,
    SimplifyTopologyPreserving, Winding,
};

/// Offsets spanning the conditioning cliff: `1e6` was always fine, `1e8` is where
/// the products leave exact-integer range, and `1e15` is past the point where the
/// *coordinates themselves* can still resolve a unit feature.
const OFFSETS: [f64; 6] = [0.0, 1e6, 1e8, 1e10, 1e12, 1e15];

fn unit_square(ox: f64, oy: f64) -> Polygon {
    Polygon::new(
        LineString::from(vec![
            (ox, oy),
            (ox + 1.0, oy),
            (ox + 1.0, oy + 1.0),
            (ox, oy + 1.0),
            (ox, oy),
        ]),
        vec![],
    )
}

#[test]
fn area_is_translation_invariant() {
    for off in OFFSETS {
        let square = unit_square(off, off);
        assert!(
            (square.unsigned_area() - 1.0).abs() < 1e-9,
            "a unit square at offset {off:e} has area {}, not 1",
            square.unsigned_area()
        );

        // A triangle too, because it is the smallest ring that has an area at all
        // and the one the exact-predicate layer is built on.
        let triangle = Polygon::new(
            LineString::from(vec![
                (off, off),
                (off + 1.0, off),
                (off, off + 1.0),
                (off, off),
            ]),
            vec![],
        );
        assert!(
            (triangle.unsigned_area() - 0.5).abs() < 1e-9,
            "a unit triangle at offset {off:e} has area {}, not 0.5",
            triangle.unsigned_area()
        );

        // A hole must still subtract: the shell and the hole are summed
        // separately, so a conditioning failure in one and not the other shows up
        // here and nowhere else.
        let holed = Polygon::new(
            LineString::from(vec![
                (off, off),
                (off + 4.0, off),
                (off + 4.0, off + 4.0),
                (off, off + 4.0),
                (off, off),
            ]),
            vec![LineString::from(vec![
                (off + 1.0, off + 1.0),
                (off + 2.0, off + 1.0),
                (off + 2.0, off + 2.0),
                (off + 1.0, off + 2.0),
                (off + 1.0, off + 1.0),
            ])],
        );
        assert!(
            (holed.unsigned_area() - 15.0).abs() < 1e-9,
            "a 4x4 square with a unit hole at offset {off:e} has area {}, not 15",
            holed.unsigned_area()
        );
    }
}

#[test]
fn winding_is_translation_invariant() {
    // The sign of the same sum, which makes this the sharper test: an area that is
    // merely inaccurate is still usually the right sign, but past the cliff the
    // sum is noise and the sign is a coin toss. Winding decides shell versus hole
    // during normalisation, so getting it wrong inverts geometry rather than
    // nudging a number.
    for off in OFFSETS {
        let ccw = unit_square(off, off);
        assert!(
            ccw.exterior().is_ccw(),
            "the counter-clockwise unit square at offset {off:e} is not reported ccw"
        );

        let mut reversed: Vec<Coord> = ccw.exterior().0.clone();
        reversed.reverse();
        assert!(
            !LineString::new(reversed).is_ccw(),
            "the clockwise unit square at offset {off:e} is reported ccw"
        );
    }
}

#[test]
fn the_centroid_is_translation_covariant() {
    // The centroid accumulates first moments, so its terms go as `|coord|³` against
    // an answer of size `|coord|·area` — one power worse conditioned than the area.
    // Translating the input must translate the answer by exactly the same amount.
    for off in OFFSETS {
        let c = unit_square(off, off)
            .centroid()
            .expect("a square has a centroid");
        let (want_x, want_y) = (off + 0.5, off + 0.5);
        // Relative to the offset: at 1e15 the coordinate itself cannot resolve 0.5,
        // so an absolute bound would be testing `f64`, not the algorithm.
        let scale = off.abs().max(1.0);
        assert!(
            (c.0.x - want_x).abs() <= 1e-9 * scale && (c.0.y - want_y).abs() <= 1e-9 * scale,
            "the centroid of a unit square at offset {off:e} is {:?}, not ({want_x}, {want_y})",
            c.0
        );
    }
}

#[test]
fn polygonize_recovers_every_face_far_from_the_origin() {
    // The failure that led to all of the above. `trace` keeps a traced face when
    // the signed area of its ring is positive, so once that sum is noise the faces
    // are discarded and every edge is reported as bounding nothing. A 3x3 grid of
    // unit cells plus a diagonal gives 12 faces totalling area 9 — the diagonal
    // halves three of the cells — and that must hold wherever the grid sits.
    let grid = |off: f64| -> Vec<LineString> {
        let mut out = Vec::new();
        for i in 0..4 {
            let t = f64::from(i);
            out.push(LineString::from(vec![(off, off + t), (off + 3.0, off + t)]));
            out.push(LineString::from(vec![(off + t, off), (off + t, off + 3.0)]));
        }
        out.push(LineString::from(vec![(off, off), (off + 3.0, off + 3.0)]));
        out
    };

    for off in [0.0_f64, 1e6, 1e8, 1e9, 1e12] {
        let got = grid(off).as_slice().polygonize();
        let total: f64 = got.polygons.iter().map(Area::unsigned_area).sum();
        assert_eq!(
            got.polygons.iter().count(),
            12,
            "offset {off:e}: expected 12 faces, got {} (with {} dangles and {} cut edges) \
             — zero faces is the signature of the shoelace sign going to noise",
            got.polygons.iter().count(),
            got.dangles.len(),
            got.cut_edges.len()
        );
        assert!(
            (total - 9.0).abs() < 1e-6,
            "offset {off:e}: the faces total {total}, not 9"
        );
        assert!(
            got.is_fully_closed(),
            "offset {off:e}: a closed grid left {} dangles and {} cut edges",
            got.dangles.len(),
            got.cut_edges.len()
        );
    }
}

#[test]
fn repair_and_simplification_survive_a_distant_origin() {
    // `make_valid` normalises ring winding from the same sign, and
    // `simplify_topology` compares the area before and against after to decide
    // whether a removal is acceptable. Both read the sum this module is about.
    for off in [0.0_f64, 1e8, 1e10, 1e12] {
        // A bowtie: self-intersecting, so repair has real work to do, and its two
        // lobes have a known total area of 1.
        let bowtie = Polygon::new(
            LineString::from(vec![
                (off, off),
                (off + 2.0, off + 2.0),
                (off + 2.0, off),
                (off, off + 2.0),
                (off, off),
            ]),
            vec![],
        );
        let repaired = bowtie.make_valid().expect("a bowtie is repairable");
        assert!(
            (repaired.unsigned_area() - 2.0).abs() < 1e-6,
            "offset {off:e}: the repaired bowtie has area {}, not 2",
            repaired.unsigned_area()
        );

        // Simplification must not change an area it cannot measure.
        let square = MultiPolygon::new(vec![unit_square(off, off)]);
        let simplified = square.simplify_topology_preserving(1e-9);
        assert!(
            (simplified.unsigned_area() - 1.0).abs() < 1e-6,
            "offset {off:e}: simplifying a unit square gave area {}, not 1",
            simplified.unsigned_area()
        );
    }
}

#[test]
fn set_operations_are_translation_invariant() {
    // The sweep has its own copy of the ring-area sum, used when it assembles
    // output rings into polygons. Two unit squares overlapping in a quarter give
    // areas known in advance at every offset.
    for off in [0.0_f64, 1e8, 1e10, 1e12] {
        let a = MultiPolygon::new(vec![unit_square(off, off)]);
        let b = MultiPolygon::new(vec![unit_square(off + 0.5, off + 0.5)]);
        for (name, got, want) in [
            ("union", a.union(&b).unsigned_area(), 1.75),
            ("intersection", a.intersection(&b).unsigned_area(), 0.25),
            ("difference", a.difference(&b).unsigned_area(), 0.75),
        ] {
            assert!(
                (got - want).abs() < 1e-6,
                "offset {off:e}: {name} of two overlapping unit squares is {got}, not {want}"
            );
        }
    }
}

/// Two triangles meeting at exactly one point, at offset `off` and feature size
/// `s`. A point touch is legal for a `MultiPolygon`; a *line* touch is not, so
/// this is the smallest geometry that can tell the two apart.
///
/// Both triangles put an edge on the line `y = x` on opposite sides of the shared
/// vertex, which is what makes it adversarial: each one's far vertex is exactly
/// collinear with the other's edge, so the only thing separating "point touch"
/// from "line touch" is whether that vertex is correctly judged to lie *beyond*
/// the end of the edge.
fn touching_triangles(off: f64, s: f64) -> MultiPolygon {
    MultiPolygon::new(vec![
        Polygon::new(
            LineString::from(vec![
                (off + s, off + s),
                (off, off + 2.0 * s),
                (off, off),
                (off + s, off + s),
            ]),
            vec![],
        ),
        Polygon::new(
            LineString::from(vec![
                (off + s, off + s),
                (off + 2.0 * s, off),
                (off + 2.0 * s, off + 2.0 * s),
                (off + s, off + s),
            ]),
            vec![],
        ),
    ])
}

#[test]
fn a_point_touch_is_not_a_line_touch_far_from_the_origin() {
    for off in [0.0_f64, 1e6, 1e8, 1e10, 1e12] {
        let pair = touching_triangles(off, 1.0);
        assert!(
            pair.is_valid(),
            "offset {off:e}: two triangles meeting at one point were rejected as {:?}",
            pair.validation_errors()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
        );
        assert!(
            (pair.unsigned_area() - 2.0).abs() < 1e-6,
            "offset {off:e}: the pair has area {}, not 2",
            pair.unsigned_area()
        );
    }
}

#[test]
fn predicates_survive_very_small_features() {
    // The other direction: an absolute floor on an area tolerance swallows a
    // geometry smaller than the floor. `1e-9` of a degree is about 0.1 mm, so this
    // is a precision a caller can legitimately be working at.
    for s in [1.0_f64, 1e-3, 1e-6, 1e-9, 1e-12] {
        let pair = touching_triangles(0.0, s);
        assert!(
            pair.is_valid(),
            "feature size {s:e}: two triangles meeting at one point were rejected as {:?}",
            pair.validation_errors()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
        );
        // And the areas must scale as s², which is the scaling half of the symmetry.
        let want = 2.0 * s * s;
        assert!(
            (pair.unsigned_area() - want).abs() <= 1e-9 * want.max(f64::MIN_POSITIVE),
            "feature size {s:e}: area {} is not {want}",
            pair.unsigned_area()
        );
    }
}

#[test]
fn areas_scale_with_the_square_of_the_input() {
    // Pure scaling, with no offset, so a failure here is a floor rather than a
    // conditioning problem.
    for s in [1e-9_f64, 1e-6, 1e-3, 1.0, 1e3, 1e6, 1e9] {
        let square = Polygon::new(
            LineString::from(vec![(0.0, 0.0), (s, 0.0), (s, s), (0.0, s), (0.0, 0.0)]),
            vec![],
        );
        let got = square.unsigned_area();
        assert!(
            (got - s * s).abs() <= 1e-12 * (s * s),
            "a square of side {s:e} has area {got}, not {}",
            s * s
        );
        assert!(
            square.exterior().is_ccw(),
            "a square of side {s:e} is not reported counter-clockwise"
        );
    }
}
