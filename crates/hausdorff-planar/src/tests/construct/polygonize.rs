//! Polygonize: faces recovered from a line soup, and the leftovers named.
//!
//! The oracle is **closed-form area**, because a grid of lines encloses cells
//! whose areas are known before any code runs. A 2x2 grid of unit cells is four
//! polygons of area 1; nothing about that needs a reference implementation.
//!
//! The second oracle is the **round trip against the boundary operation**:
//! polygonizing a polygon's own ring edges must give that polygon back. That ties
//! this module to the geometry model rather than only to itself.
//!
//! The cases that matter most are the ones where the input does *not* close, since
//! that is what real data does. A dangle and a cut edge are both "this edge bounds
//! no face", and the suite checks they are distinguished — the difference tells a
//! user what to fix, so collapsing them into one bucket would make the diagnostic
//! useless.

use crate::space::r2::validation::Validation;
use crate::space::r2::{Area, Contains, Coord, LineString, Point, Polygon, Polygonize};

fn line(pts: &[(f64, f64)]) -> LineString {
    LineString::from(pts.to_vec())
}

fn close(got: f64, want: f64, what: &str) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1.0),
        "{what}: got {got}, want {want}"
    );
}

#[test]
fn a_single_closed_ring_becomes_one_polygon() {
    let ring = [line(&[
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (0.0, 0.0),
    ])];
    let out = ring.polygonize();
    assert_eq!(out.polygons.0.len(), 1, "one ring, one face");
    close(out.polygons.unsigned_area(), 16.0, "a 4x4 square");
    assert!(
        out.is_fully_closed(),
        "a closed ring leaves nothing over, got {:?} / {:?}",
        out.dangles,
        out.cut_edges
    );
    assert!(out.polygons.is_valid(), "the output must be valid geometry");
}

#[test]
fn four_lines_forming_a_square_close_into_it() {
    // The actual use case: the boundary arrives as separate lines, drawn in
    // whatever direction the surveyor walked, and must still close.
    let edges = [
        line(&[(0.0, 0.0), (2.0, 0.0)]),
        line(&[(2.0, 2.0), (2.0, 0.0)]), // drawn backwards on purpose
        line(&[(2.0, 2.0), (0.0, 2.0)]),
        line(&[(0.0, 0.0), (0.0, 2.0)]),
    ];
    let out = edges.polygonize();
    assert_eq!(out.polygons.0.len(), 1, "four edges enclose one face");
    close(out.polygons.unsigned_area(), 4.0, "a 2x2 square");
    assert!(out.is_fully_closed(), "nothing should be left over");
}

#[test]
fn a_grid_of_lines_yields_every_cell() {
    // Closed-form oracle: a 3x3 lattice of lines over [0,2]^2 bounds four unit
    // cells. This is the shape of a street grid becoming a block layer.
    let mut edges = Vec::new();
    for i in 0..3 {
        let t = f64::from(i);
        edges.push(line(&[(0.0, t), (2.0, t)]));
        edges.push(line(&[(t, 0.0), (t, 2.0)]));
    }
    let out = edges.polygonize();
    assert_eq!(out.polygons.0.len(), 4, "a 2x2 grid has four cells");
    close(out.polygons.unsigned_area(), 4.0, "total area");
    for poly in out.polygons.iter() {
        close(poly.unsigned_area(), 1.0, "each cell is a unit square");
        assert!(poly.interiors().is_empty(), "a grid cell has no hole");
    }
    assert!(out.is_fully_closed());
}

#[test]
fn lines_crossing_without_a_shared_vertex_are_noded_first() {
    // The reason the arrangement is needed rather than plain ring-finding: these
    // two lines cross in their interiors and share no vertex, so a naive
    // graph walk would see four loose ends and find nothing.
    let edges = [
        line(&[(0.0, 0.0), (2.0, 2.0)]),
        line(&[(0.0, 2.0), (2.0, 0.0)]),
        line(&[(0.0, 0.0), (0.0, 2.0)]),
        line(&[(2.0, 0.0), (2.0, 2.0)]),
        line(&[(0.0, 0.0), (2.0, 0.0)]),
        line(&[(0.0, 2.0), (2.0, 2.0)]),
    ];
    let out = edges.polygonize();
    // The two diagonals cut the square into four triangles.
    assert_eq!(
        out.polygons.0.len(),
        4,
        "the crossing diagonals must be noded, giving four triangles, got {}",
        out.polygons.0.len()
    );
    close(out.polygons.unsigned_area(), 4.0, "total area is preserved");
    for poly in out.polygons.iter() {
        close(poly.unsigned_area(), 1.0, "each triangle is a quarter");
    }
}

#[test]
fn a_dangle_is_reported_and_not_silently_dropped() {
    // A line stopping short of the boundary it was meant to meet. The face still
    // closes; the spur has a loose end and bounds nothing.
    let edges = [
        line(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]),
        line(&[(2.0, 2.0), (3.0, 3.0)]), // free-floating spur inside the square
    ];
    let out = edges.polygonize();
    assert_eq!(out.polygons.0.len(), 1, "the square still closes");
    close(out.polygons.unsigned_area(), 16.0, "and keeps its area");
    assert!(
        !out.dangles.is_empty(),
        "the spur must be reported as a dangle, not discarded"
    );
    assert!(
        !out.is_fully_closed(),
        "a dangle means the input did not fully close"
    );
}

#[test]
fn a_chain_of_spurs_peels_off_entirely() {
    // Removing a dangle can expose another, so the peel has to iterate. A
    // single-pass implementation would leave the inner links of this chain behind
    // and then report them as cut edges, which is the wrong diagnosis.
    let edges = [
        line(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]),
        line(&[(2.0, 2.0), (2.5, 2.5)]),
        line(&[(2.5, 2.5), (3.0, 3.0)]),
        line(&[(3.0, 3.0), (3.5, 3.5)]),
    ];
    let out = edges.polygonize();
    assert_eq!(out.polygons.0.len(), 1);
    close(out.polygons.unsigned_area(), 16.0, "the face is unaffected");
    assert!(
        out.cut_edges.is_empty(),
        "a chain of spurs is all dangles; none should be misreported as cut edges: {:?}",
        out.cut_edges
    );
    assert!(
        out.dangles.len() >= 3,
        "every link of the chain is a dangle"
    );
}

#[test]
fn a_bridge_between_two_loops_is_reported_as_a_cut_edge() {
    // Both ends connected, yet it bounds no face: removing it would disconnect the
    // network. A different fix from a dangle, hence a different bucket.
    let edges = [
        // Left square.
        line(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]),
        // Right square, detached.
        line(&[(3.0, 0.0), (4.0, 0.0), (4.0, 1.0), (3.0, 1.0), (3.0, 0.0)]),
        // The connector.
        line(&[(1.0, 0.0), (3.0, 0.0)]),
    ];
    let out = edges.polygonize();
    assert_eq!(out.polygons.0.len(), 2, "two faces close");
    close(out.polygons.unsigned_area(), 2.0, "two unit squares");
    assert!(
        !out.cut_edges.is_empty(),
        "the connector bounds nothing and both its ends are joined, so it is a cut edge"
    );
    assert!(
        out.dangles.is_empty(),
        "the connector has no loose end, so it must not be called a dangle: {:?}",
        out.dangles
    );
}

#[test]
fn a_ring_inside_a_ring_becomes_a_shell_with_a_hole() {
    // Depth classification, which is what distinguishes this from naive
    // ring-finding: the inner ring is a hole, not a second polygon.
    let edges = [
        line(&[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]),
        line(&[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)]),
    ];
    let out = edges.polygonize();
    assert_eq!(
        out.polygons.0.len(),
        1,
        "a nested pair is one polygon with a hole, not two polygons"
    );
    assert_eq!(out.polygons.0[0].interiors().len(), 1, "exactly one hole");
    close(
        out.polygons.unsigned_area(),
        100.0 - 16.0,
        "shell minus hole",
    );
    assert!(
        !out.polygons.0[0].contains(&Point(Coord::new(5.0, 5.0))),
        "the hole interior is outside the polygon"
    );
    assert!(
        out.polygons.0[0].contains(&Point(Coord::new(1.0, 1.0))),
        "the ring between them is inside"
    );
    assert!(out.is_fully_closed());
}

#[test]
fn three_levels_of_nesting_give_an_island_in_a_lake() {
    // Depth parity, not just "inside something": the innermost ring is at depth 2,
    // so it is a shell again — an island in a lake — rather than a hole of a hole.
    let edges = [
        line(&[
            (0.0, 0.0),
            (12.0, 0.0),
            (12.0, 12.0),
            (0.0, 12.0),
            (0.0, 0.0),
        ]),
        line(&[
            (2.0, 2.0),
            (10.0, 2.0),
            (10.0, 10.0),
            (2.0, 10.0),
            (2.0, 2.0),
        ]),
        line(&[(4.0, 4.0), (8.0, 4.0), (8.0, 8.0), (4.0, 8.0), (4.0, 4.0)]),
    ];
    let out = edges.polygonize();
    assert_eq!(
        out.polygons.0.len(),
        2,
        "the outer ring and the island are both shells, got {}",
        out.polygons.0.len()
    );
    // Outer ring area 144 - 64 = 80; island 16.
    close(
        out.polygons.unsigned_area(),
        80.0 + 16.0,
        "ring plus island",
    );
    assert!(
        out.polygons.iter().any(|p| p.interiors().len() == 1),
        "the outer shell carries the lake as a hole"
    );
    assert!(out.polygons.is_valid());
}

#[test]
fn polygonizing_a_polygons_own_boundary_returns_it() {
    // Ties the module to the geometry model: polygonize is the inverse of taking a
    // boundary, so this round trip must hold.
    let originals = [
        Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (5.0, 0.0),
                (5.0, 5.0),
                (0.0, 5.0),
                (0.0, 0.0),
            ]),
            vec![],
        ),
        Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 10.0),
                (0.0, 10.0),
                (0.0, 0.0),
            ]),
            vec![LineString::from(vec![
                (2.0, 2.0),
                (4.0, 2.0),
                (4.0, 4.0),
                (2.0, 4.0),
                (2.0, 2.0),
            ])],
        ),
    ];
    for original in &originals {
        let out = crate::space::r2::Geometry::Polygon(original.clone()).polygonize();
        close(
            out.polygons.unsigned_area(),
            original.unsigned_area(),
            "round-tripped area",
        );
        assert_eq!(out.polygons.0.len(), 1, "one polygon in, one polygon out");
        assert_eq!(
            out.polygons.0[0].interiors().len(),
            original.interiors().len(),
            "holes must survive the round trip"
        );
        assert!(
            out.is_fully_closed(),
            "a closed boundary leaves nothing over"
        );
    }
}

#[test]
fn a_duplicated_boundary_does_not_become_a_sliver() {
    // A shared parcel edge digitized from both sides arrives twice. Without
    // deduplication the arrangement would carry it as two coincident edges and the
    // trace would find a zero-area face between them.
    let edges = [
        line(&[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]),
        line(&[(0.0, 0.0), (3.0, 0.0)]), // the bottom edge again
        line(&[(3.0, 0.0), (0.0, 0.0)]), // and once more, reversed
    ];
    let out = edges.polygonize();
    assert_eq!(out.polygons.0.len(), 1, "still one face");
    close(out.polygons.unsigned_area(), 9.0, "and the right area");
    for poly in out.polygons.iter() {
        assert!(
            poly.unsigned_area() > 1e-9,
            "no zero-area sliver should appear"
        );
    }
}

#[test]
fn degenerate_and_hostile_input_is_handled() {
    use crate::space::r2::{Geometry, MultiPoint};

    // Nothing to polygonize.
    for empty in [
        Vec::<LineString>::new(),
        vec![line(&[])],
        vec![line(&[(1.0, 1.0)])],
        vec![line(&[(1.0, 1.0), (1.0, 1.0)])],
    ] {
        let out = empty.polygonize();
        assert!(
            out.polygons.0.is_empty(),
            "nothing enclosed, so no polygons: {out:?}"
        );
    }

    // An open chain closes nothing but must not panic, and every edge is a dangle.
    let open = [line(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])];
    let out = open.polygonize();
    assert!(out.polygons.0.is_empty(), "an open chain encloses nothing");
    assert!(!out.dangles.is_empty(), "its edges are dangles");

    // Points contribute no edges.
    let pts = [Geometry::MultiPoint(MultiPoint::new(vec![
        Point(Coord::new(0.0, 0.0)),
        Point(Coord::new(1.0, 1.0)),
    ]))];
    assert!(pts.polygonize().polygons.0.is_empty());

    // Non-finite coordinates cannot be noded; an empty result rather than a panic.
    let nasty = [line(&[(0.0, 0.0), (f64::NAN, 1.0), (1.0, 1.0), (0.0, 0.0)])];
    let out = nasty.polygonize();
    assert!(
        out.polygons.0.is_empty(),
        "a non-finite coordinate makes the input unusable, got {out:?}"
    );
}
