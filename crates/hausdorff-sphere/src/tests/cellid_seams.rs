//! `S2CellId` where the cube meets the sphere: face corners and edges, the level
//! extremes, and the bit-level invariants the whole index rests on.
//!
//! An `S2CellId` is not an opaque handle — the hierarchy is *encoded in the bits*,
//! and that is the only reason the index is useful. `parent` is a mask,
//! containment is an integer range comparison, and numeric order is Hilbert order,
//! so a database that can range-scan can do spatial queries. Every one of those
//! properties is a claim about bit layout, and a bit-layout claim either holds
//! exactly or is worthless. None of them admit a tolerance.
//!
//! So the oracles here are **structural identities**, checked exhaustively where
//! the space is small enough to allow it:
//!
//! - `contains(x, y)` must agree with `range_min(x) <= y <= range_max(x)`, because
//!   the second is what callers will actually use and the first is what the docs
//!   promise. Two ways of asking the same question, one of them a comparison on
//!   `u64`.
//! - `parent(child(c)) == c` at every level, and `level(parent(c, n)) == n`.
//! - The four children of a cell partition it: their ranges tile the parent's with
//!   no gap and no overlap.
//! - `next`/`prev` walk Hilbert order, so `prev(next(c)) == c` away from the ends.
//!
//! The **six cube faces** are where this is most likely to break. S2 projects the
//! sphere onto a cube, and a face boundary is a seam in the parameterization even
//! though it is nothing at all on the sphere. A point exactly on a face edge or at
//! a face corner belongs to one face by convention, and the neighbour queries have
//! to cross the seam correctly — a corner is adjacent to three faces at once.
//! Nothing about the sphere marks these places, which is exactly why an
//! implementation can be wrong there and right everywhere else.

use crate::space::s2::S2CellId;
use hausdorff_base::earth::GeoPoint;

/// The deepest S2 level: a leaf cell, about 1 cm across.
const MAX_LEVEL: u32 = 30;

fn g(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

/// Points sitting on the seams of the cube parameterization and of lat/lon.
fn seam_points() -> Vec<(&'static str, GeoPoint)> {
    vec![
        ("origin", g(0.0, 0.0)),
        ("north pole", g(90.0, 0.0)),
        ("south pole", g(-90.0, 0.0)),
        ("antimeridian east", g(0.0, 180.0)),
        ("antimeridian west", g(0.0, -180.0)),
        ("face centre +x", g(0.0, 0.0)),
        ("face centre +y", g(0.0, 90.0)),
        ("face centre -x", g(0.0, 180.0)),
        ("face centre -y", g(0.0, -90.0)),
        // Cube corners sit where three faces meet: lat = atan(1/sqrt(2)).
        ("cube corner NE", g(35.264_389_682_754_654, 45.0)),
        ("cube corner NW", g(35.264_389_682_754_654, -45.0)),
        ("cube corner SE", g(-35.264_389_682_754_654, 45.0)),
        ("cube corner SW", g(-35.264_389_682_754_654, -45.0)),
        ("face edge x/y", g(0.0, 45.0)),
        ("near the pole", g(89.999_999, 0.0)),
        ("negative zero", g(-0.0, -0.0)),
    ]
}

#[test]
fn the_six_faces_are_valid_distinct_and_self_identifying() {
    let mut seen = std::collections::HashSet::new();
    for face in 0..6u8 {
        let id = S2CellId::from_face(face);
        assert!(id.is_valid(), "face {face} must be a valid cell id");
        assert!(
            id.is_face(),
            "face {face} must report itself as a face cell"
        );
        assert_eq!(id.level(), 0, "a face cell is at level 0");
        assert_eq!(
            id.face(),
            face,
            "face {face} must round-trip through face()"
        );
        assert!(!id.is_leaf(), "a face cell is not a leaf");
        assert!(
            seen.insert(id.as_u64()),
            "face {face} duplicates another face"
        );
    }
    // The six faces partition the sphere, so every point lands on exactly one and
    // each face contains the cell of any point assigned to it.
    for (name, p) in seam_points() {
        let leaf = S2CellId::from_geo_point(p);
        assert!(leaf.is_valid(), "{name}: point produced an invalid cell");
        let owners: Vec<u8> = (0..6u8)
            .filter(|&f| S2CellId::from_face(f).contains(leaf))
            .collect();
        assert_eq!(
            owners.len(),
            1,
            "{name}: should belong to exactly one face, got {owners:?}"
        );
        assert_eq!(
            owners[0],
            leaf.face(),
            "{name}: containment and face() must agree"
        );
    }
}

#[test]
fn contains_agrees_with_the_range_comparison() {
    // The two ways of asking containment must be the same question. Callers use
    // the range form (it is an integer comparison, so it pushes into a query
    // engine); the docs promise the method form.
    for (name, p) in seam_points() {
        let leaf = S2CellId::from_geo_point(p);
        for level in 0..=MAX_LEVEL {
            let ancestor = leaf.parent(level);
            assert!(
                ancestor.contains(leaf),
                "{name}: the level-{level} ancestor must contain its descendant"
            );
            let (lo, hi) = (ancestor.range_min(), ancestor.range_max());
            assert!(
                lo.as_u64() <= leaf.as_u64() && leaf.as_u64() <= hi.as_u64(),
                "{name}: level {level} range [{}, {}] excludes {}",
                lo.as_u64(),
                hi.as_u64(),
                leaf.as_u64()
            );
            // And the range endpoints are themselves inside.
            assert!(ancestor.contains(lo) && ancestor.contains(hi));
            assert!(
                ancestor.intersects(leaf),
                "{name}: containment implies intersection"
            );
        }
    }
}

#[test]
fn parent_and_child_are_mutual_inverses_at_every_level() {
    for (name, p) in seam_points() {
        let leaf = S2CellId::from_geo_point(p);
        assert_eq!(
            leaf.level(),
            MAX_LEVEL,
            "{name}: from_geo_point gives a leaf"
        );
        assert!(leaf.is_leaf());

        for level in 0..MAX_LEVEL {
            let here = leaf.parent(level);
            assert_eq!(
                here.level(),
                level,
                "{name}: parent(level) must be at level"
            );
            // The child containing our leaf, found among the four.
            let child = here
                .children()
                .into_iter()
                .find(|c| c.contains(leaf))
                .unwrap_or_else(|| panic!("{name}: no child of level {level} contains the leaf"));
            assert_eq!(child.level(), level + 1);
            assert_eq!(
                child.immediate_parent().as_u64(),
                here.as_u64(),
                "{name}: immediate_parent must undo the descent at level {level}"
            );
            assert_eq!(
                child.parent(level).as_u64(),
                here.as_u64(),
                "{name}: parent(level) must undo the descent at level {level}"
            );
        }
    }
}

#[test]
fn the_four_children_tile_their_parent_exactly() {
    // A partition: no gap, no overlap, and the union is the parent's whole range.
    // A gap would silently lose geometry in a covering; an overlap would
    // double-count it.
    for face in 0..6u8 {
        let mut cell = S2CellId::from_face(face);
        for level in 0..8u32 {
            let kids = cell.children();
            assert_eq!(kids.len(), 4, "an S2 cell always has four children");
            let mut ranges: Vec<(u64, u64)> = kids
                .iter()
                .map(|c| (c.range_min().as_u64(), c.range_max().as_u64()))
                .collect();
            ranges.sort_unstable();
            // Contiguous and non-overlapping. Note the ranges are expressed in
            // *leaf* ids, and a leaf id always carries its trailing one bit at
            // position 0 — so leaf ids are odd and consecutive leaves differ by 2,
            // not 1. Stepping with `next()` rather than `+ 1` says that without
            // hardcoding the stride.
            for w in ranges.windows(2) {
                let end_of_previous = S2CellId::from_u64(w[0].1);
                assert!(
                    end_of_previous.is_leaf(),
                    "range_max should be a leaf cell id"
                );
                assert_eq!(
                    end_of_previous.next().as_u64(),
                    w[1].0,
                    "face {face} level {level}: the leaf after one child's range_max \
                     must be the next child's range_min"
                );
            }
            // Together they are the parent's range.
            assert_eq!(
                ranges[0].0,
                cell.range_min().as_u64(),
                "face {face} level {level}: children must start where the parent does"
            );
            assert_eq!(
                ranges[3].1,
                cell.range_max().as_u64(),
                "face {face} level {level}: children must end where the parent does"
            );
            for (i, k) in kids.iter().enumerate() {
                assert_eq!(k.level(), level + 1);
                assert!(cell.contains(*k), "child {i} must be inside its parent");
                assert_eq!(
                    k.child_position(level + 1) as usize,
                    // position within the parent, 0..4
                    kids.iter().position(|x| x.as_u64() == k.as_u64()).unwrap(),
                    "child_position must match the enumeration order"
                );
            }
            cell = kids[1];
        }
    }
}

#[test]
fn hilbert_order_walks_forward_and_back() {
    // `next`/`prev` step along the space-filling curve. Being mutual inverses is
    // what makes a range scan well-defined.
    for face in 0..6u8 {
        let start = S2CellId::from_face(face).parent(0).children()[0].parent(1);
        let mut cell = start;
        for step in 0..64u32 {
            let forward = cell.next();
            if !forward.is_valid() {
                break;
            }
            assert_eq!(
                forward.prev().as_u64(),
                cell.as_u64(),
                "face {face} step {step}: prev must undo next"
            );
            assert_eq!(
                forward.level(),
                cell.level(),
                "next must stay at the same level"
            );
            assert!(
                forward.as_u64() > cell.as_u64(),
                "next must increase the id, since numeric order is Hilbert order"
            );
            cell = forward;
        }
    }
}

#[test]
fn hilbert_order_keeps_neighbours_numerically_close() {
    // The property that makes S2 a partitioning tool rather than just a geometry
    // library. Without it a range scan is not a spatial query.
    let level = 14;
    let here = S2CellId::from_geo_point_at_level(g(51.5074, -0.1278), level);
    let near = S2CellId::from_geo_point_at_level(g(51.5075, -0.1279), level);
    let far = S2CellId::from_geo_point_at_level(g(-33.8688, 151.2093), level);
    let gap_near = here.as_u64().abs_diff(near.as_u64());
    let gap_far = here.as_u64().abs_diff(far.as_u64());
    assert!(
        gap_near < gap_far,
        "a neighbouring cell must be numerically closer than an antipodal one \
         ({gap_near} vs {gap_far})"
    );
    // Adjacent cells should share a long bit prefix, i.e. a shallow common ancestor.
    let common = here
        .common_ancestor_level(near)
        .expect("two cells on the same sphere share an ancestor");
    assert!(
        common >= 8,
        "cells ~10 m apart at level {level} should share an ancestor at level 8 or \
         deeper, got {common}"
    );
}

#[test]
fn neighbours_are_symmetric_and_adjacent_across_face_seams() {
    // A cube corner touches three faces, so this is where a neighbour query has to
    // cross a seam that exists only in the parameterization.
    for (name, p) in seam_points() {
        for level in [1u32, 5, 10, 20] {
            let cell = S2CellId::from_geo_point_at_level(p, level);
            let edge = cell.edge_neighbors();
            assert_eq!(
                edge.len(),
                4,
                "{name}: a quad cell has four edge neighbours"
            );
            for (i, n) in edge.iter().enumerate() {
                assert!(n.is_valid(), "{name}: edge neighbour {i} is invalid");
                assert_eq!(n.level(), level, "{name}: a neighbour is at the same level");
                assert_ne!(
                    n.as_u64(),
                    cell.as_u64(),
                    "{name}: a cell is not its own neighbour"
                );
                // Symmetry: adjacency is a symmetric relation, even across a seam.
                assert!(
                    n.edge_neighbors()
                        .iter()
                        .any(|m| m.as_u64() == cell.as_u64()),
                    "{name}: {cell:?} lists {n:?} as a neighbour but not conversely"
                );
            }
            // The four must be distinct.
            let unique: std::collections::HashSet<u64> = edge.iter().map(|c| c.as_u64()).collect();
            assert_eq!(
                unique.len(),
                4,
                "{name}: the four edge neighbours must be distinct at level {level}"
            );
            // All-neighbours is a superset of the edge neighbours.
            let all: std::collections::HashSet<u64> = cell
                .all_neighbors(level)
                .iter()
                .map(|c| c.as_u64())
                .collect();
            for n in &edge {
                assert!(
                    all.contains(&n.as_u64()),
                    "{name}: all_neighbors must include every edge neighbour"
                );
            }
            assert!(
                !all.contains(&cell.as_u64()),
                "{name}: all_neighbors must not include the cell itself"
            );
        }
    }
}

#[test]
fn a_cell_contains_its_own_centre_and_vertices_bound_it() {
    for (name, p) in seam_points() {
        for level in [0u32, 1, 5, 15, 25, MAX_LEVEL] {
            let cell = S2CellId::from_geo_point_at_level(p, level);
            // The centre must land back in the same cell — the basic sanity of the
            // cube-to-sphere projection and its inverse.
            let centre = cell.center();
            assert_eq!(
                S2CellId::from_geo_point_at_level(centre, level).as_u64(),
                cell.as_u64(),
                "{name}: the centre of a level-{level} cell must re-encode to it"
            );
            // Four vertices, all unit vectors.
            let vs = cell.vertices();
            assert_eq!(vs.len(), 4);
            for (i, v) in vs.iter().enumerate() {
                let n = v.as_vector().norm();
                assert!(
                    (n - 1.0).abs() < 1e-12,
                    "{name}: vertex {i} of a level-{level} cell is not a unit vector (norm {n})"
                );
            }
        }
    }
}

#[test]
fn tokens_round_trip_and_reject_garbage() {
    for (name, p) in seam_points() {
        for level in [0u32, 7, 18, MAX_LEVEL] {
            let cell = S2CellId::from_geo_point_at_level(p, level);
            let token = cell.to_token();
            assert!(!token.is_empty(), "{name}: a valid cell has a token");
            assert_eq!(
                S2CellId::from_token(&token).map(|c| c.as_u64()),
                Some(cell.as_u64()),
                "{name}: token {token} must round-trip at level {level}"
            );
        }
    }
    for bad in [
        "",
        " ",
        "zzzzzzzzzzzzzzzzz",
        "not-hex",
        "0x1234",
        "🙂",
        "-1",
        "G",
    ] {
        assert!(
            S2CellId::from_token(bad).is_none(),
            "{bad:?} is not a valid cell token"
        );
    }
}

#[test]
fn invalid_bit_patterns_are_rejected() {
    // A `u64` is not a cell id. Zero has no trailing one bit to mark a level, and
    // a face index above 5 does not exist.
    assert!(
        !S2CellId::from_u64(0).is_valid(),
        "zero is not a valid cell"
    );
    assert!(
        !S2CellId::from_u64(u64::MAX).is_valid(),
        "all-ones is not a valid cell"
    );
    // Walk the high bits: faces 6 and 7 are out of range for a 3-bit face field.
    for face in 6u64..8 {
        let id = S2CellId::from_u64((face << 61) | 1);
        assert!(!id.is_valid(), "face index {face} does not exist");
    }
    // A valid leaf stays valid through a u64 round trip.
    let leaf = S2CellId::from_geo_point(g(12.0, 34.0));
    assert!(S2CellId::from_u64(leaf.as_u64()).is_valid());
    assert_eq!(S2CellId::from_u64(leaf.as_u64()).as_u64(), leaf.as_u64());
}

#[test]
fn cells_at_the_same_level_never_overlap() {
    // Disjointness within a level is what makes a covering a partition. Sampled
    // across the sphere rather than exhaustively, since a level-10 enumeration is
    // six million cells.
    let level = 10;
    let mut cells = Vec::new();
    let mut lat = -80.0;
    while lat <= 80.0 {
        let mut lon = -180.0;
        while lon < 180.0 {
            cells.push(S2CellId::from_geo_point_at_level(g(lat, lon), level));
            lon += 37.0;
        }
        lat += 23.0;
    }
    for (i, a) in cells.iter().enumerate() {
        for b in &cells[i + 1..] {
            if a.as_u64() == b.as_u64() {
                continue;
            }
            assert!(
                !a.intersects(*b),
                "distinct level-{level} cells {a:?} and {b:?} must not intersect"
            );
            assert!(!a.contains(*b) && !b.contains(*a));
        }
    }
}
