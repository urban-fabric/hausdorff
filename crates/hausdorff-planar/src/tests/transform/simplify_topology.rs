//! Topology-preserving simplification, against the thing plain RDP gets wrong.
//!
//! The oracle here is **validity**, which is unusually strong for a simplifier:
//! the crate already has an independent validator, so "the output is valid
//! geometry" is checkable without any reference implementation and without a
//! tolerance. If a chord crosses the ring it came from, `is_valid` says so.
//!
//! The sharpest case is the **comparison against plain `simplify`**. A test that
//! only showed the topology-preserving variant producing valid output would not
//! show it was needed — plain RDP produces valid output most of the time too. So
//! the suite constructs a ring where plain RDP demonstrably self-intersects and
//! then requires the guarded variant not to. If plain RDP ever stops failing that
//! case, the test says so rather than passing quietly.

use crate::space::r2::validation::Validation;
use crate::space::r2::{
    Area, Coord, Geometry, LineString, MultiPolygon, Polygon, Simplify, SimplifyTopologyPreserving,
};

fn ring(pts: &[(f64, f64)]) -> LineString {
    let mut v: Vec<Coord> = pts.iter().map(|&(x, y)| Coord::new(x, y)).collect();
    if v.first() != v.last() {
        v.push(v[0]);
    }
    LineString::new(v)
}

/// A tight spiral ribbon: three turns, densely sampled, with the arms 0.4 apart.
///
/// This is the fixture that actually defeats plain RDP, and finding one took a
/// search rather than a guess — a deep narrow bay, the intuitive candidate, does
/// *not* break it, because RDP keeps the vertex farthest from the chord and that
/// vertex is the bottom of the bay. A spiral is different: every arm is nearly
/// straight, so RDP straightens each one, and a straightened arm cuts across its
/// neighbour. The 0.4 gap is what makes the arms close enough for the chords to
/// cross at a tolerance of 1.
///
/// Verified below to be valid as built and invalid after plain `simplify`, so the
/// comparison it anchors is real rather than assumed.
fn tight_spiral() -> Polygon {
    let (turns, samples, gap) = (3.0_f64, 400usize, 0.4_f64);
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for i in 0..=samples {
        let t = i as f64 / samples as f64;
        let th = t * turns * std::f64::consts::TAU;
        let r = 1.0 + 9.0 * t;
        pts.push((r * th.cos(), r * th.sin()));
    }
    // Back along a slightly smaller radius, closing the ribbon.
    for i in (0..=samples).rev() {
        let t = i as f64 / samples as f64;
        let th = t * turns * std::f64::consts::TAU;
        let r = 1.0 + 9.0 * t - gap;
        pts.push((r * th.cos(), r * th.sin()));
    }
    Polygon::new(ring(&pts), vec![])
}

#[test]
fn a_simplified_ring_stays_valid_where_plain_rdp_does_not() {
    let original = tight_spiral();
    assert!(original.is_valid(), "the fixture itself must be valid");

    // Large enough to tempt RDP into collapsing the spike.
    let epsilon = 1.0;
    let plain = original.simplify(epsilon);
    let guarded = original.simplify_topology_preserving(epsilon);

    assert!(
        guarded.is_valid(),
        "the topology-preserving result must be valid, got {guarded:?}"
    );
    // The point of the exercise: plain RDP is expected to break here. If it does
    // not, this fixture has stopped demonstrating anything and should be replaced
    // rather than silently passing.
    assert!(
        !plain.is_valid(),
        "this fixture is supposed to be one plain RDP mangles; if plain RDP now \
         keeps it valid the test is no longer testing the difference"
    );
}

#[test]
fn output_is_valid_across_tolerances_and_fixtures() {
    // A sweep rather than a case: the veto has to hold at every tolerance, not
    // just the one the fixture was tuned for. A tolerance larger than the geometry
    // is the interesting end — that is where everything wants to collapse.
    let fixtures: Vec<Polygon> = vec![
        tight_spiral(),
        // A simple concave ring with a narrow inlet — valid as built, and the
        // intuitive candidate for breaking RDP even though it does not.
        Polygon::new(
            ring(&[
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 10.0),
                (5.2, 10.0),
                (5.2, 0.6),
                (4.9, 0.5),
                (4.6, 0.6),
                (4.6, 10.0),
                (0.0, 10.0),
            ]),
            vec![],
        ),
        // A polygon with a hole, so the hole rings are exercised too.
        Polygon::new(
            ring(&[(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)]),
            vec![ring(&[
                (5.0, 5.0),
                (15.0, 5.0),
                (15.0, 6.0),
                (6.0, 6.0),
                (6.0, 15.0),
                (5.0, 15.0),
            ])],
        ),
        // A near-degenerate sliver: the case that collapses or inverts.
        Polygon::new(
            ring(&[
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 0.001),
                (5.0, 0.0005),
                (0.0, 0.001),
            ]),
            vec![],
        ),
    ];
    for (i, poly) in fixtures.iter().enumerate() {
        assert!(poly.is_valid(), "fixture {i} must start valid");
        let before = poly.unsigned_area();
        for epsilon in [0.0, 1e-6, 0.01, 0.5, 1.0, 5.0, 100.0] {
            let out = poly.simplify_topology_preserving(epsilon);
            assert!(
                out.is_valid(),
                "fixture {i} at epsilon {epsilon} became invalid: {out:?}"
            );
            let after = out.unsigned_area();
            assert!(
                after > 0.0,
                "fixture {i} at epsilon {epsilon} collapsed to zero area"
            );
            // Orientation is preserved, so the signed areas keep their sign. Checked
            // via unsigned area staying in a sane band rather than exactly: the
            // whole point of simplifying is to change the area a little.
            assert!(
                after <= before * 4.0 + 1.0,
                "fixture {i} at epsilon {epsilon} ballooned from {before} to {after}, \
                 which means a ring inverted"
            );
        }
    }
}

#[test]
fn a_ring_is_never_reduced_below_a_triangle() {
    // A "simplified" ring with no area is a deleted feature, not a smaller one. At
    // an absurd tolerance everything wants to collapse, and the floor is what keeps
    // the output a polygon.
    let poly = Polygon::new(
        ring(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (0.0, 3.0),
        ]),
        vec![],
    );
    for epsilon in [10.0, 1e3, 1e9] {
        let out = poly.simplify_topology_preserving(epsilon);
        assert!(
            out.exterior().0.len() >= 4,
            "a ring must keep at least three distinct vertices at epsilon {epsilon}, \
             got {}",
            out.exterior().0.len()
        );
        assert!(out.unsigned_area() > 0.0, "and must still bound area");
        assert!(out.is_valid());
    }
}

#[test]
fn simplification_actually_removes_vertices_when_it_safely_can() {
    // The veto must not be so conservative that the operation does nothing. A
    // convex ring with many nearly-collinear vertices has no way to
    // self-intersect, so every removal within tolerance should be taken.
    let mut pts: Vec<(f64, f64)> = Vec::new();
    let mut i = 0;
    while i <= 100 {
        let t = f64::from(i) / 100.0;
        // A gentle arc: dense vertices, tiny deviation from the chord.
        pts.push((t * 10.0, (t * std::f64::consts::PI).sin() * 0.05));
        i += 1;
    }
    pts.push((10.0, 5.0));
    pts.push((0.0, 5.0));
    let poly = Polygon::new(ring(&pts), vec![]);
    let before = poly.exterior().0.len();
    let out = poly.simplify_topology_preserving(0.5);
    assert!(
        out.exterior().0.len() < before / 4,
        "a smooth convex arc should simplify hard: {} vertices from {before}",
        out.exterior().0.len()
    );
    assert!(out.is_valid());
}

#[test]
fn an_open_line_is_plain_douglas_peucker() {
    // Documented behaviour: an open line has no ring to self-intersect against, so
    // the guarded variant gives the same answer as the plain one. Pinned so a
    // caller is not surprised to find the two differ.
    let line = LineString::from(vec![
        (0.0, 0.0),
        (1.0, 0.1),
        (2.0, -0.1),
        (3.0, 0.05),
        (4.0, 0.0),
        (5.0, 3.0),
    ]);
    for epsilon in [0.0, 0.2, 1.0] {
        assert_eq!(
            line.simplify_topology_preserving(epsilon).0,
            line.simplify(epsilon).0,
            "an open line gets plain RDP at epsilon {epsilon}"
        );
    }
}

#[test]
fn multipolygons_and_collections_are_simplified_per_part() {
    let mp = MultiPolygon::new(vec![tight_spiral(), {
        let mut p = tight_spiral();
        // Translate the second part well clear of the first.
        p = Polygon::new(
            LineString::new(
                p.exterior()
                    .0
                    .iter()
                    .map(|c| Coord::new(c.x + 100.0, c.y))
                    .collect(),
            ),
            vec![],
        );
        p
    }]);
    assert!(mp.is_valid(), "the fixture must start valid");
    let out = mp.simplify_topology_preserving(1.0);
    assert_eq!(out.0.len(), 2, "both parts survive");
    assert!(out.is_valid(), "and the result is valid: {out:?}");

    let gc = Geometry::GeometryCollection(vec![
        Geometry::Polygon(tight_spiral()),
        Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 0.2), (2.0, 0.0)])),
        Geometry::Point(crate::space::r2::Point(Coord::new(3.0, 3.0))),
    ]);
    let simplified = gc.simplify_topology_preserving(0.5);
    match simplified {
        Geometry::GeometryCollection(parts) => {
            assert_eq!(parts.len(), 3, "every part is carried through");
            assert!(
                matches!(parts[2], Geometry::Point(_)),
                "a point has nothing to simplify and passes through unchanged"
            );
        }
        other => panic!("a collection must stay a collection, got {other:?}"),
    }
}

#[test]
fn degenerate_input_is_returned_rather_than_mangled() {
    // Too few vertices to simplify, and non-finite coordinates. Neither should
    // panic, and neither should produce something worse than it started.
    for pts in [
        vec![(0.0, 0.0)],
        vec![(0.0, 0.0), (1.0, 0.0)],
        vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
    ] {
        let p = Polygon::new(ring(&pts), vec![]);
        let out = p.simplify_topology_preserving(1.0);
        assert_eq!(
            out.exterior().0.len(),
            p.exterior().0.len(),
            "a minimal ring has nothing to remove"
        );
    }
    let nasty = Polygon::new(
        ring(&[
            (0.0, 0.0),
            (f64::NAN, 1.0),
            (2.0, 2.0),
            (0.0, 3.0),
            (1.0, 4.0),
        ]),
        vec![],
    );
    // Only requirement: it returns. A NaN cannot be compared, so no guarantee is
    // made about which vertices survive.
    let _ = nasty.simplify_topology_preserving(0.5);
}
