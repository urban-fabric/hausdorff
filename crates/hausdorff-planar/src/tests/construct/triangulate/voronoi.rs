//! The Voronoi dual: cell ownership, unbounded hull cells, and nearest-site
//! lookup, plus the same adversarial inputs run through the dual.

use super::*;

#[test]
fn voronoi_cells_contain_only_points_nearer_their_own_site() {
    // The defining property of a Voronoi cell, sampled over each closed cell: a
    // point inside cell `i` must be nearer site `i` than any other site.
    let points = scatter(30, 50.0, 47);
    let v = voronoi(&points);
    assert!(
        v.cells.iter().any(Option::is_some),
        "an interior site should get a closed cell"
    );
    for (i, cell) in v.cells.iter().enumerate() {
        let Some(poly) = cell else { continue };
        let site = v.sites[i];
        // Sample the cell: its own vertices pulled toward the site, which are
        // inside by convexity.
        for vert in poly.exterior().0.iter() {
            let probe = Coord::new(
                site.x + (vert.x - site.x) * 0.5,
                site.y + (vert.y - site.y) * 0.5,
            );
            let own = (probe.x - site.x).powi(2) + (probe.y - site.y).powi(2);
            for (j, other) in v.sites.iter().enumerate() {
                if j == i {
                    continue;
                }
                let d = (probe.x - other.x).powi(2) + (probe.y - other.y).powi(2);
                assert!(
                    d >= own - 1e-6 * own.max(1.0),
                    "a point in cell {i} is nearer site {j}: {d} < {own}"
                );
            }
        }
    }
}

#[test]
fn hull_sites_have_no_bounded_cell_and_say_so() {
    // The honest `None`: a hull site's true cell is unbounded, so there is no
    // polygon. Pinned because inventing a clipped rectangle instead would be
    // geometry the caller never asked for.
    let points = pts(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (5.0, 5.0),
    ]);
    let v = voronoi(&points);
    assert_eq!(v.sites.len(), 5);
    // The centre is the only interior site.
    let closed: Vec<usize> = v
        .cells
        .iter()
        .enumerate()
        .filter(|(_, c)| c.is_some())
        .map(|(i, _)| i)
        .collect();
    for &i in &closed {
        assert!(
            (v.sites[i].x - 5.0).abs() < 1e-9 && (v.sites[i].y - 5.0).abs() < 1e-9,
            "only the interior site should have a bounded cell, got site {:?}",
            v.sites[i]
        );
    }
    // And the clip gives every site a cell.
    let bounds = crate::space::r2::Rect::new(Coord::new(-5.0, -5.0), Coord::new(15.0, 15.0));
    let clipped = v.clipped_to(bounds);
    assert_eq!(clipped.len(), 5);
    assert!(
        clipped.iter().filter(|c| c.is_some()).count() >= closed.len(),
        "clipping must not lose the cells that already closed"
    );
}

#[test]
fn nearest_site_agrees_with_a_brute_force_scan() {
    // The question a Voronoi diagram is built to answer, checked against the
    // scan it is supposed to replace.
    let points = scatter(25, 100.0, 53);
    let v = voronoi(&points);
    for probe in scatter(40, 100.0, 59) {
        let found = nearest_site(&v, probe).expect("non-empty");
        let brute = v
            .sites
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                let da = (a.x - probe.x).powi(2) + (a.y - probe.y).powi(2);
                let db = (b.x - probe.x).powi(2) + (b.y - probe.y).powi(2);
                da.total_cmp(&db)
            })
            .map(|(i, _)| i)
            .expect("non-empty");
        let d_found = (v.sites[found].x - probe.x).powi(2) + (v.sites[found].y - probe.y).powi(2);
        let d_brute = (v.sites[brute].x - probe.x).powi(2) + (v.sites[brute].y - probe.y).powi(2);
        // Compare distances rather than indices: equidistant sites are a genuine
        // tie and either answer is right.
        assert!(
            (d_found - d_brute).abs() <= 1e-6 * d_brute.max(1.0),
            "nearest_site picked {found} at {d_found} but {brute} is at {d_brute}"
        );
    }
}

#[test]
fn voronoi_survives_the_same_adversarial_inputs() {
    // The dual inherits whatever the triangulation got wrong, so the cases above
    // are run through it too. The oracle is the defining property: a point in a
    // cell is nearer that cell's site than any other.
    let mut lattice = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            lattice.push(Coord::new(f64::from(i), f64::from(j)));
        }
    }
    let circle: Vec<Coord> = (0..12)
        .map(|i| {
            let a = std::f64::consts::TAU * f64::from(i) / 12.0;
            Coord::new(a.cos() * 5.0, a.sin() * 5.0)
        })
        .collect();
    let mut circle_centre = circle.clone();
    circle_centre.push(Coord::new(0.0, 0.0));

    for (name, points) in [
        ("lattice", lattice),
        ("circle", circle),
        ("circle plus centre", circle_centre),
        ("scatter", scatter(30, 50.0, 137)),
    ] {
        let v = voronoi(&points);
        assert_eq!(
            v.sites.len(),
            v.cells.len(),
            "{name}: one cell slot per site"
        );
        for (i, cell) in v.cells.iter().enumerate() {
            let Some(poly) = cell else { continue };
            let site = v.sites[i];
            for vert in &poly.exterior().0 {
                // Halfway from the site to a cell vertex is inside the cell, since
                // a Voronoi cell is convex.
                let probe = Coord::new(
                    site.x + (vert.x - site.x) * 0.5,
                    site.y + (vert.y - site.y) * 0.5,
                );
                let own = (probe.x - site.x).powi(2) + (probe.y - site.y).powi(2);
                for (j, other) in v.sites.iter().enumerate() {
                    if j == i {
                        continue;
                    }
                    let d = (probe.x - other.x).powi(2) + (probe.y - other.y).powi(2);
                    assert!(
                        d >= own - 1e-6 * own.max(1.0),
                        "{name}: a point in cell {i} is nearer site {j}"
                    );
                }
            }
        }
    }
}
