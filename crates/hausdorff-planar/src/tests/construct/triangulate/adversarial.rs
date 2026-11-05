//! Inputs chosen to break the in-circle predicate: cocircular points, lattices,
//! near-duplicates, and extreme magnitudes.
//!
//! This is the suite that found the cavity-selection defect — 16 points on a
//! circle produced 22 triangles where 14 exist.

use super::*;

#[test]
fn cocircular_input_does_not_produce_overlapping_triangles() {
    // The defect this suite found. Points on a common circle make every in-circle
    // determinant near zero, so they round inconsistently; a cavity chosen by
    // global scan is then disconnected and the re-fan overlaps. Growing the cavity
    // outward from a seed removes the failure mode rather than tightening a
    // tolerance.
    let circle: Vec<Coord> = (0..16)
        .map(|i| {
            let a = std::f64::consts::TAU * f64::from(i) / 16.0;
            Coord::new(a.cos() * 10.0, a.sin() * 10.0)
        })
        .collect();
    // All 16 are on the hull, so Euler gives 2*16 - 2 - 16 = 14 exactly.
    let t = assert_sound("16 points on a circle", &circle);
    assert_eq!(t.triangles.len(), 14, "a 16-gon fans into 14 triangles");

    let mut with_centre = circle.clone();
    with_centre.push(Coord::new(0.0, 0.0));
    let t = assert_sound("circle plus centre", &with_centre);
    assert_eq!(
        t.triangles.len(),
        16,
        "adding the centre gives one triangle per circle edge"
    );

    // Other cocircular families, at several sizes, since the failure depended on
    // how many near-zero determinants had to agree.
    for sides in [3usize, 4, 5, 6, 7, 8, 12, 24, 31] {
        let ring: Vec<Coord> = (0..sides)
            .map(|i| {
                let a = std::f64::consts::TAU * (i as f64) / (sides as f64);
                Coord::new(a.cos(), a.sin())
            })
            .collect();
        let t = assert_sound(&format!("regular {sides}-gon"), &ring);
        assert_eq!(
            t.triangles.len(),
            sides - 2,
            "a convex {sides}-gon fans into {} triangles",
            sides - 2
        );
    }
}

#[test]
fn lattices_triangulate_soundly() {
    // A lattice is the other massively-cocircular case: every 2x2 block shares a
    // circumcircle. It is also the case where Euler's `h` counts points *along*
    // the hull edges, not just the four corners.
    for side in [2usize, 3, 4, 5, 7] {
        let mut lattice = Vec::new();
        for i in 0..side {
            for j in 0..side {
                lattice.push(Coord::new(i as f64, j as f64));
            }
        }
        assert_sound(&format!("{side}x{side} lattice"), &lattice);
    }
    // A lattice with one point nudged off it, which breaks the symmetry the
    // cocircular case relies on.
    let mut jittered = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            jittered.push(Coord::new(f64::from(i), f64::from(j)));
        }
    }
    jittered.push(Coord::new(2.0 + 1e-9, 2.0));
    assert_sound("lattice plus a nudged point", &jittered);
}

#[test]
fn near_duplicate_points_are_merged_rather_than_confusing_the_predicate() {
    // Two points a few ulp apart are distinct to `==` but indistinguishable to a
    // degree-4 determinant, so keeping both asks the predicate a question it
    // cannot answer. They are merged within a relative tolerance instead.
    let mut points = scatter(20, 10.0, 101);
    let base = points[0];
    for k in 1..5 {
        points.push(Coord::new(base.x + 1e-12 * f64::from(k), base.y));
    }
    let t = assert_sound("scatter plus near-duplicates", &points);
    assert_eq!(
        t.points.len(),
        20,
        "the near-duplicates must collapse into the point they shadow"
    );

    // Exact duplicates, en masse.
    let repeated: Vec<Coord> = (0..30)
        .map(|i| {
            [
                Coord::new(0.0, 0.0),
                Coord::new(1.0, 0.0),
                Coord::new(0.0, 1.0),
            ][i % 3]
        })
        .collect();
    let t = assert_sound("three points, thirty copies", &repeated);
    assert_eq!(t.points.len(), 3);
    assert_eq!(t.triangles.len(), 1);
}

#[test]
fn extreme_shapes_and_magnitudes_triangulate_soundly() {
    // A long thin strip: every point is on the hull, so Euler gives 2n - 2 - n.
    let mut strip = Vec::new();
    for i in 0..40 {
        strip.push(Coord::new(f64::from(i) * 100.0, 0.0));
        strip.push(Coord::new(f64::from(i) * 100.0, 1e-4));
    }
    assert_sound("long thin strip", &strip);

    // A collinear run embedded in a scatter: the run's points are interior, and a
    // zero-area fan triangle would show up as an area mismatch.
    let mut mixed = scatter(20, 50.0, 103);
    for i in 0..10 {
        mixed.push(Coord::new(f64::from(i) * 5.0, 25.0));
    }
    assert_sound("scatter plus a collinear run", &mixed);

    // Magnitudes spanning eighteen orders, which is what the centroid translation
    // in both the builder and the checker exists for.
    for (name, shift, scale, seed) in [
        ("at the origin", 0.0, 1.0, 107),
        ("geographic", -122.0, 0.5, 109),
        ("projected", 551_000.0, 5000.0, 113),
        ("at 1e9", 1e9, 1.0, 127),
        ("tiny", 0.0, 1e-9, 131),
    ] {
        let points: Vec<Coord> = scatter(25, scale, seed)
            .into_iter()
            .map(|c| Coord::new(c.x + shift, c.y + shift))
            .collect();
        assert_sound(&format!("scatter {name}"), &points);
    }
}
