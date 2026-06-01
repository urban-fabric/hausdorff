//! Tests for the S2 cell-indexing port: structural invariants, reference
//! tokens, and lat/lon round-trips.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use hausdorff::geo::GeoPoint;

use hausdorff::s2::S2CellId;

const MAX_LEVEL: u32 = 30;

fn gp(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

#[test]
fn face_cells_are_valid_level_zero() {
    for f in 0..6u8 {
        let c = S2CellId::from_face(f);
        assert!(c.is_valid());
        assert_eq!(c.face(), f);
        assert_eq!(c.level(), 0);
        assert!(c.is_face());
    }
}

#[test]
fn zero_and_max_are_invalid() {
    assert!(!S2CellId::from_u64(0).is_valid());
    assert!(!S2CellId::from_u64(u64::MAX).is_valid());
}

#[test]
fn face_tokens_match_reference() {
    // Level-0 face cells are `(2f+1) << 60`, i.e. tokens "1","3",…,"b".
    let expected = ["1", "3", "5", "7", "9", "b"];
    for (f, tok) in expected.iter().enumerate() {
        assert_eq!(S2CellId::from_face(f as u8).to_token(), *tok);
    }
}

#[test]
fn point_lands_on_expected_face() {
    // The +X,+Y,+Z,-X,-Y,-Z axes map to faces 0..5 respectively.
    assert_eq!(S2CellId::from_geo_point(gp(0.0, 0.0)).face(), 0); // +X
    assert_eq!(S2CellId::from_geo_point(gp(0.0, 90.0)).face(), 1); // +Y
    assert_eq!(S2CellId::from_geo_point(gp(90.0, 0.0)).face(), 2); // +Z (north pole)
    assert_eq!(S2CellId::from_geo_point(gp(0.0, 180.0)).face(), 3); // -X
    assert_eq!(S2CellId::from_geo_point(gp(0.0, -90.0)).face(), 4); // -Y
    assert_eq!(S2CellId::from_geo_point(gp(-90.0, 0.0)).face(), 5); // -Z (south pole)
}

#[test]
fn lat_lng_maps_to_leaf_cell() {
    let c = S2CellId::from_geo_point(gp(37.7749, -122.4194));
    assert!(c.is_valid());
    assert!(c.is_leaf());
    assert_eq!(c.level(), MAX_LEVEL);
}

#[test]
fn center_round_trips_within_leaf_resolution() {
    // A leaf cell is sub-centimeter, so its center is within a few 1e-7 deg of
    // any point it contains.
    for &(lat, lon) in &[
        (0.0, 0.0),
        (37.7749, -122.4194),
        (-33.8688, 151.2093),
        (51.5074, -0.1278),
        (-89.0, 179.0),
        (89.0, -179.0),
    ] {
        let p = gp(lat, lon);
        let back = S2CellId::from_geo_point(p).center();
        assert!((back.lat - lat).abs() < 1e-5, "lat {lat} -> {}", back.lat);
        assert!((back.lon - lon).abs() < 1e-5, "lon {lon} -> {}", back.lon);
    }
}

#[test]
fn token_round_trips() {
    for &(lat, lon) in &[(0.0, 0.0), (37.77, -122.41), (-12.3, 45.6)] {
        let c = S2CellId::from_geo_point(gp(lat, lon));
        assert_eq!(S2CellId::from_token(&c.to_token()), Some(c));
    }
    assert_eq!(S2CellId::from_token(""), None);
    assert_eq!(S2CellId::from_token("zzz"), None);
    assert_eq!(S2CellId::from_token("00000000000000000"), None); // > 16 chars
}

#[test]
fn parent_chain_contains_and_tracks_level() {
    let leaf = S2CellId::from_geo_point(gp(37.7749, -122.4194));
    let mut prev = leaf;
    for level in (0..MAX_LEVEL).rev() {
        let parent = leaf.parent(level);
        assert_eq!(parent.level(), level);
        assert!(parent.contains(leaf));
        assert_eq!(prev.immediate_parent(), parent);
        prev = parent;
    }
    assert!(prev.is_face());
}

#[test]
fn children_partition_parent() {
    let parent = S2CellId::from_geo_point_at_level(gp(10.0, 20.0), 10);
    let kids = parent.children();
    for (pos, k) in kids.iter().enumerate() {
        assert_eq!(k.level(), 11);
        assert!(parent.contains(*k));
        assert_eq!(k.parent(10), parent);
        assert_eq!(k.child_position(11) as usize, pos);
    }
    for a in 0..4 {
        for b in (a + 1)..4 {
            assert_ne!(kids[a], kids[b]);
        }
    }
}

#[test]
fn edge_neighbors_are_same_level_and_distinct() {
    let c = S2CellId::from_geo_point_at_level(gp(37.7749, -122.4194), 12);
    let n = c.edge_neighbors();
    for nb in n {
        assert!(nb.is_valid());
        assert_eq!(nb.level(), 12);
        assert_ne!(nb, c);
        assert!(!c.contains(nb));
    }
    for a in 0..4 {
        for b in (a + 1)..4 {
            assert_ne!(n[a], n[b]);
        }
    }
}

#[test]
fn all_neighbors_include_edge_neighbors() {
    let c = S2CellId::from_geo_point_at_level(gp(37.7749, -122.4194), 12);
    let all = c.all_neighbors(12);
    for nb in &all {
        assert!(nb.is_valid());
        assert_eq!(nb.level(), 12);
    }
    for nb in &c.edge_neighbors() {
        assert!(all.contains(nb), "edge neighbor {nb} missing from all_neighbors");
    }
    assert!(all_neighbors_empty_when_level_too_coarse(c));
}

fn all_neighbors_empty_when_level_too_coarse(c: S2CellId) -> bool {
    // No neighbors exist at a level coarser than the cell's own level.
    c.all_neighbors(c.level().saturating_sub(1)).is_empty()
}

#[test]
fn display_renders_face_and_child_path() {
    // S2's canonical "face/child-path" form keeps the slash even for a face
    // cell, so a level-0 cell renders as "1/" (matching the Go reference).
    assert_eq!(S2CellId::from_face(1).to_string(), "1/");
    assert_eq!(S2CellId::from_face(1).children()[2].to_string(), "1/2");
}

#[test]
fn common_ancestor_level_of_self_is_own_level() {
    let a = S2CellId::from_geo_point_at_level(gp(10.0, 20.0), 20);
    assert_eq!(a.common_ancestor_level(a), Some(20));
    let sibling = a.parent(10).children()[0];
    let ancestor = a.common_ancestor_level(sibling);
    assert!(ancestor.is_some_and(|l| l <= 10));
}
