//! The twelve pentagons, and the grid identities that hold around them.
//!
//! H3 tiles the sphere by projecting an icosahedron, and an icosahedron has
//! twelve vertices. A hexagonal grid cannot close around a vertex of degree five,
//! so at each of those twelve places one cell is a **pentagon** and the grid is
//! locally wrong-by-construction. Every traversal, hierarchy, and area operation
//! has to special-case them, which makes pentagons the single most likely place
//! for a port of H3 to be subtly wrong.
//!
//! The oracles here are **counting identities**, which is what makes them worth
//! having: each one is a closed form in `k`, derived rather than transcribed, so a
//! wrong number cannot be mistaken for a tolerance problem.
//!
//! | shape | ring at distance `k` | disk out to `k` |
//! |---|---|---|
//! | hexagon | `6k` | `1 + 3k(k+1)` |
//! | pentagon | `5k` | `1 + 5k(k+1)/2` |
//!
//! **The ideal counts only hold where the disk contains no pentagon.** This is the
//! trap, and it is why the hexagon cases below are measured at fine resolutions in
//! open ocean rather than on a base cell. At resolution 0 there are 122 cells and
//! 12 of them are pentagons, so almost every base cell has a pentagon within two
//! steps and its disk is legitimately *short* of the ideal — 18 instead of 19 at
//! `k = 2`. A test that asserted the ideal count on a base cell would be asserting
//! that pentagons do not exist, and "fixing" the implementation to satisfy it would
//! break the grid. So each case states which regime it is in, and the
//! pentagon-adjacency is asserted rather than assumed.

use crate::index::h3::H3Cell;

/// A point in the South Atlantic, far from every icosahedron vertex. At fine
/// resolutions its neighbourhood is pure hexagon, which is the regime where the
/// ideal counts apply.
const OPEN_OCEAN: (f64, f64) = (0.0, -30.0);

#[test]
fn there_are_exactly_twelve_pentagons_among_122_base_cells() {
    let base: Vec<H3Cell> = H3Cell::base_cells().collect();
    assert_eq!(base.len(), 122, "H3 has 122 resolution-0 base cells");
    let pentagons = base.iter().filter(|c| c.is_pentagon()).count();
    assert_eq!(
        pentagons, 12,
        "one pentagon per icosahedron vertex, so exactly twelve"
    );
}

#[test]
fn a_pentagon_stays_a_pentagon_down_the_hierarchy() {
    // The pentagon is at a vertex of the icosahedron, so the centre child of a
    // pentagon is itself a pentagon at every finer resolution. Its other children
    // are ordinary hexagons.
    for p in H3Cell::base_cells().filter(|c| c.is_pentagon()) {
        let mut cell = p;
        for res in 1..=8u8 {
            let kids = cell.children_at_resolution(res);
            let pent_kids = kids.iter().filter(|c| c.is_pentagon()).count();
            assert_eq!(
                pent_kids, 1,
                "exactly one descendant of a pentagon at resolution {res} is a pentagon"
            );
            cell = *kids
                .iter()
                .find(|c| c.is_pentagon())
                .expect("the centre child");
        }
    }
}

#[test]
fn pentagon_descendant_counts_follow_the_closed_form() {
    // A hexagon has 7 children, a pentagon 6 (one of the seven subsequences is
    // deleted), and only the pentagon's centre child is itself a pentagon. So the
    // count at resolution r is 1 + 5·(7^r − 1)/6 rather than 7^r.
    for p in H3Cell::base_cells().filter(|c| c.is_pentagon()) {
        for res in 1..=6u32 {
            let expected = 1 + 5 * (7u64.pow(res) - 1) / 6;
            let got = p.children_at_resolution(res as u8).len() as u64;
            assert_eq!(
                got, expected,
                "pentagon descendants at resolution {res} should be {expected}, got {got}"
            );
            assert_eq!(
                p.children_count(res as u8),
                expected,
                "children_count must agree with the enumeration at resolution {res}"
            );
        }
    }
    // And the hexagon case, for contrast: a plain 7^r.
    let hex = H3Cell::base_cells()
        .find(|c| !c.is_pentagon())
        .expect("a hexagonal base cell");
    for res in 1..=6u32 {
        assert_eq!(
            hex.children_at_resolution(res as u8).len() as u64,
            7u64.pow(res),
            "a hexagon has 7^r descendants"
        );
    }
}

#[test]
fn pentagon_rings_have_five_k_cells_and_hexagon_rings_six_k() {
    let pent = H3Cell::base_cells()
        .find(|c| c.is_pentagon())
        .expect("a pentagon");
    for k in 1..5u32 {
        let ring = pent
            .grid_disk_distances(k)
            .into_iter()
            .filter(|(_, d)| *d == k)
            .count();
        assert_eq!(
            ring as u32,
            5 * k,
            "the ring at distance {k} from a pentagon holds 5k cells"
        );
    }

    let hex = H3Cell::for_point(OPEN_OCEAN.0, OPEN_OCEAN.1, 8).expect("a fine-resolution cell");
    for k in 1..5u32 {
        let disk = hex.grid_disk_distances(k);
        assert!(
            !disk.iter().any(|(c, _)| c.is_pentagon()),
            "this neighbourhood is supposed to be pentagon-free; the ideal count \
             does not hold otherwise"
        );
        let ring = disk.into_iter().filter(|(_, d)| *d == k).count();
        assert_eq!(
            ring as u32,
            6 * k,
            "the ring at distance {k} from an interior hexagon holds 6k cells"
        );
    }
}

#[test]
fn disk_counts_match_the_ideal_only_where_no_pentagon_intrudes() {
    // The pentagon-free regime: the ideal hexagonal disk, exactly, at several
    // resolutions so this is not a one-cell coincidence.
    for res in [5u8, 8, 10, 13] {
        let c = H3Cell::for_point(OPEN_OCEAN.0, OPEN_OCEAN.1, res).expect("valid cell");
        for k in 0..5u32 {
            let disk = c.grid_disk(k);
            assert!(
                !disk.iter().any(|c| c.is_pentagon()),
                "resolution {res} k={k} was expected to be pentagon-free"
            );
            assert_eq!(
                disk.len() as u32,
                1 + 3 * k * (k + 1),
                "resolution {res}, k={k}: an interior hexagon's disk is 1 + 3k(k+1)"
            );
        }
    }

    // The pentagon's own disk: 1 + 5k(k+1)/2, five-fold rather than six-fold.
    let pent = H3Cell::base_cells()
        .find(|c| c.is_pentagon())
        .expect("a pentagon");
    for k in 0..5u32 {
        assert_eq!(
            pent.grid_disk(k).len() as u32,
            1 + 5 * k * (k + 1) / 2,
            "a pentagon's disk out to k={k} is 1 + 5k(k+1)/2"
        );
    }

    // The distorted regime, asserted as such: a resolution-0 hexagon next to a
    // pentagon has *fewer* than the ideal, and that is correct behaviour.
    let distorted = H3Cell::base_cells()
        .find(|c| !c.is_pentagon() && c.grid_disk(2).iter().any(|c| c.is_pentagon()))
        .expect("a base-cell hexagon with a pentagon in its 2-disk");
    assert!(
        (distorted.grid_disk(2).len() as u32) < 1 + 3 * 2 * 3,
        "a disk reaching a pentagon must come up short of the ideal count, not meet it"
    );
}

#[test]
fn disks_never_repeat_a_cell_even_around_a_pentagon() {
    // A traversal that walks into the pentagon's missing sector and back out is the
    // way a k-ring implementation starts emitting duplicates.
    let mut cells: Vec<H3Cell> = H3Cell::base_cells().filter(|c| c.is_pentagon()).collect();
    cells.push(H3Cell::for_point(OPEN_OCEAN.0, OPEN_OCEAN.1, 9).expect("valid"));
    for c in cells {
        for k in 0..6u32 {
            let disk = c.grid_disk(k);
            let unique: std::collections::HashSet<H3Cell> = disk.iter().copied().collect();
            assert_eq!(
                unique.len(),
                disk.len(),
                "grid_disk({k}) around {c:?} emitted a duplicate"
            );
            // Every cell in the disk must be at the right resolution and valid.
            for d in &disk {
                assert_eq!(d.resolution(), c.resolution());
            }
        }
    }
}

#[test]
fn grid_disk_is_symmetric_as_a_neighbour_relation() {
    // `b` within k of `a` must imply `a` within k of `b`. Adjacency is symmetric
    // even where the grid is not uniform, and a pentagon is exactly where an
    // asymmetric implementation would show up.
    let pent = H3Cell::base_cells()
        .find(|c| c.is_pentagon())
        .expect("a pentagon");
    // Walk down to a finer resolution so the neighbourhood is bigger than the
    // whole base-cell grid.
    let mut cell = pent;
    for res in 1..=4u8 {
        cell = *cell
            .children_at_resolution(res)
            .iter()
            .find(|c| c.is_pentagon())
            .expect("centre child");
    }
    for k in 1..4u32 {
        for neighbour in cell.grid_disk(k) {
            assert!(
                neighbour.grid_disk(k).contains(&cell),
                "{cell:?} reaches {neighbour:?} within {k} but not the reverse"
            );
        }
    }
}
