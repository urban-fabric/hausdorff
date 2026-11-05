//! The Delaunay properties themselves: counts, tiling, the empty-circumcircle
//! guarantee, orientation, and the nearest-neighbour subgraph.

use super::*;

#[test]
fn a_square_triangulates_into_two_triangles() {
    // Four cocircular points: the Delaunay triangulation is *not unique*, both
    // diagonals are valid, and the count is what is determined. Euler: 2*4-2-4 = 2.
    let t = delaunay(&pts(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]));
    assert_eq!(t.triangles.len(), 2, "a square is two triangles");
    assert_eq!(t.points.len(), 4);
    let area: f64 = t
        .to_triangles()
        .iter()
        .map(|x| x.to_polygon().unsigned_area())
        .sum();
    assert!(
        (area - 1.0).abs() < 1e-12,
        "the triangles tile the square, got {area}"
    );
    assert!(
        t.is_delaunay(),
        "cocircular points are *on* the circle, not inside"
    );
}

#[test]
fn triangle_count_matches_eulers_formula() {
    // 2n - 2 - h, exactly. An integer identity, so there is no tolerance to hide
    // behind: a lost or duplicated triangle fails it.
    let cases: Vec<Vec<Coord>> = vec![
        pts(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]),
        pts(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (2.0, 2.0)]),
        pts(&[
            (0.0, 0.0),
            (10.0, 1.0),
            (5.0, 9.0),
            (3.0, 3.0),
            (7.0, 4.0),
            (5.0, 6.0),
        ]),
        scatter(30, 100.0, 7),
        scatter(80, 1.0, 11),
        // Geographic magnitudes, where the untranslated determinant loses bits.
        scatter(40, 0.5, 13)
            .into_iter()
            .map(|c| Coord::new(c.x - 122.4, c.y + 37.8))
            .collect(),
        // Projected magnitudes, worse still.
        scatter(40, 5000.0, 17)
            .into_iter()
            .map(|c| Coord::new(c.x + 551_000.0, c.y + 4_182_000.0))
            .collect(),
    ];
    for (i, points) in cases.iter().enumerate() {
        let t = delaunay(points);
        let n = t.points.len();
        let h = hull_count(&t.points);
        let expected = 2 * n - 2 - h;
        assert_eq!(
            t.triangles.len(),
            expected,
            "case {i}: {n} points with {h} on the hull should give {expected} triangles, \
             got {}",
            t.triangles.len()
        );
    }
}

#[test]
fn the_triangles_tile_the_convex_hull_exactly() {
    // Area conservation, which catches overlaps and gaps that the triangle count
    // can miss when two errors cancel.
    for (i, points) in [
        scatter(40, 100.0, 3),
        scatter(60, 2.0, 5),
        scatter(25, 1e4, 19),
    ]
    .iter()
    .enumerate()
    {
        let t = delaunay(points);
        let (sum, hull) = areas_translated(&t);
        assert!(
            (sum - hull).abs() <= 1e-9 * hull.max(1.0),
            "case {i}: triangles total {sum} but the hull is {hull}"
        );
    }
}

#[test]
fn the_delaunay_property_holds_on_scattered_input() {
    // The definition. Checked on enough scatters, at enough magnitudes, that the
    // coordinate translation in the implementation is actually being exercised.
    for (i, points) in [
        scatter(50, 100.0, 23),
        scatter(70, 1.0, 29),
        scatter(40, 1e6, 31),
        scatter(45, 0.001, 37),
    ]
    .iter()
    .enumerate()
    {
        let t = delaunay(points);
        assert!(
            !t.triangles.is_empty(),
            "case {i}: a scatter of {} points must triangulate",
            points.len()
        );
        assert!(
            t.is_delaunay(),
            "case {i}: the empty-circumcircle property must hold"
        );
    }
}

#[test]
fn every_triangle_is_counter_clockwise_and_non_degenerate() {
    // Orientation is part of the contract, and a zero-area triangle in the output
    // would break anything that computes a normal or interpolates over it.
    let t = delaunay(&scatter(50, 10.0, 41));
    for (i, &[a, b, c]) in t.triangles.iter().enumerate() {
        let (pa, pb, pc) = (t.points[a], t.points[b], t.points[c]);
        let cross = (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x);
        assert!(
            cross > 0.0,
            "triangle {i} is clockwise or degenerate (cross = {cross})"
        );
        assert!(a != b && b != c && a != c, "triangle {i} repeats a vertex");
    }
}

#[test]
fn the_delaunay_graph_contains_the_nearest_neighbour_of_every_point() {
    // A known property worth leaning on: the Delaunay graph is a supergraph of
    // the nearest-neighbour graph. It is an independent check on the *edges*
    // rather than the triangles, and it would fail if the triangulation were
    // merely some valid triangulation rather than the Delaunay one.
    let points = scatter(40, 50.0, 43);
    let t = delaunay(&points);
    let edges: std::collections::HashSet<(usize, usize)> = t.edges().into_iter().collect();
    for (i, &p) in t.points.iter().enumerate() {
        let nearest = t
            .points
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .min_by(|(_, a), (_, b)| {
                let da = (a.x - p.x).powi(2) + (a.y - p.y).powi(2);
                let db = (b.x - p.x).powi(2) + (b.y - p.y).powi(2);
                da.total_cmp(&db)
            })
            .map(|(j, _)| j)
            .expect("more than one point");
        let key = (i.min(nearest), i.max(nearest));
        assert!(
            edges.contains(&key),
            "the Delaunay graph must contain the edge from {i} to its nearest \
             neighbour {nearest}"
        );
    }
}

#[test]
fn degenerate_input_is_handled_rather_than_panicking() {
    // Too few points, duplicates, collinear, and non-finite — none of which has a
    // triangulation, and none of which may panic.
    for (name, input) in [
        ("empty", Vec::new()),
        ("one point", pts(&[(0.0, 0.0)])),
        ("two points", pts(&[(0.0, 0.0), (1.0, 1.0)])),
        ("all duplicates", pts(&[(1.0, 1.0), (1.0, 1.0), (1.0, 1.0)])),
        (
            "collinear",
            pts(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]),
        ),
        (
            "non-finite",
            pts(&[(0.0, 0.0), (f64::NAN, 1.0), (f64::INFINITY, 2.0)]),
        ),
    ] {
        let t = delaunay(&input);
        assert!(
            t.triangles.is_empty(),
            "{name} has no triangulation, got {} triangles",
            t.triangles.len()
        );
        assert!(
            t.is_delaunay(),
            "{name}: an empty triangulation is vacuously Delaunay"
        );
        let v = voronoi(&input);
        assert!(
            v.cells.iter().all(Option::is_none),
            "{name}: no site can have a bounded cell"
        );
    }

    // Duplicates alongside real points are deduplicated, not triangulated twice.
    let withdups = pts(&[(0.0, 0.0), (0.0, 0.0), (4.0, 0.0), (4.0, 0.0), (2.0, 3.0)]);
    let t = delaunay(&withdups);
    assert_eq!(
        t.points.len(),
        3,
        "duplicates collapse to three distinct sites"
    );
    assert_eq!(t.triangles.len(), 1, "three points are one triangle");
}

#[test]
fn the_trait_form_works_on_geometry() {
    let g = Geometry::MultiPoint(MultiPoint::new(
        pts(&[(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (0.0, 5.0), (2.0, 3.0)])
            .into_iter()
            .map(Point)
            .collect(),
    ));
    let t = g.delaunay_triangulation();
    assert!(!t.triangles.is_empty(), "the trait form must triangulate");
    assert!(t.is_delaunay());
    assert_eq!(
        t.to_multi_polygon().0.len(),
        t.triangles.len(),
        "every triangle becomes a polygon"
    );
    let v = g.voronoi_diagram();
    assert_eq!(v.sites.len(), g.to_coords().len());
    // A polygon's own vertices triangulate too, which is the TIN use case.
    let poly = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    assert!(poly.delaunay_triangulation().is_delaunay());
}

// ---- adversarial inputs ----
//
// These are the cases that broke the first implementation, kept as assertions so
// they cannot break it again. Two of them found a real defect: a global scan for
// the Bowyer-Watson cavity selected a disconnected region when many points were
// cocircular, and the re-fan then produced overlapping triangles -- 16 points on a
// circle gave 22 where only 14 exist. The cavity is now grown by breadth-first
// search across triangle adjacency, which cannot select a disconnected set.
