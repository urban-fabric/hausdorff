use super::*;

#[test]
fn parent_matches_golden() {
    let mut f = Fnv::new();
    for raw in sample_cells() {
        let ours = H3Cell::try_new(raw).expect("ours");
        let res = ours.resolution();
        for pres in 0..=res {
            let ou_p = ours.parent(pres).expect("our parent");
            // Metamorphic: the parent is coarser and the cell rebuilds from it.
            assert_eq!(ou_p.resolution(), pres, "parent res {raw:#x} -> r{pres}");
            assert!(
                ours.parent(res).expect("self parent") == ours || pres < res,
                "parent at own res must be self for {raw:#x}"
            );
            f.u64(ou_p.as_u64());
        }
        assert_eq!(ours.parent(res).map(H3Cell::as_u64), Some(raw));
        assert!(ours.parent(res + 1).is_none());
    }
    assert_eq!(f.finish(), golden::PARENT, "parent golden");
}

#[test]
fn children_match_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");
        let pres = ours.resolution();
        let tres = (pres + 2).min(15);
        let children = ours.children_at_resolution(tres);
        // Metamorphic: count matches, every child's parent is this cell, and
        // the center child is the first child.
        assert_eq!(
            ours.children_count(tres),
            children.len() as u64,
            "children_count {raw:#x}"
        );
        for child in &children {
            assert_eq!(child.resolution(), tres, "child res {raw:#x}");
            assert_eq!(child.parent(pres), Some(ours), "child parent {raw:#x}");
        }
        assert_eq!(
            ours.center_child(tres).map(H3Cell::as_u64),
            children.first().map(|c| c.as_u64()),
            "center child {raw:#x} -> r{tres}"
        );
        for child in &children {
            f.u64(child.as_u64());
        }
    }
    assert_eq!(f.finish(), golden::CHILDREN, "children golden");
}

#[test]
fn child_position_and_child_at_match_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");
        let parent_res = ours.resolution();
        let target_res = (parent_res + 2).min(15);

        let children = ours.children_at_resolution(target_res);
        for (position, child) in children.iter().copied().enumerate() {
            // Self-consistency: child_at inverts the position enumeration, and
            // child_position recovers the index.
            assert_eq!(
                ours.child_at(position as u64, target_res),
                Some(child),
                "our child_at {raw:#x} pos={position} r{target_res}"
            );
            assert_eq!(
                child.child_position(parent_res),
                Some(position as u64),
                "child_position {raw:#x} child={:#x}",
                child.as_u64()
            );
            f.u64(child.as_u64());
            f.opt_u64(
                ours.child_at(position as u64, target_res)
                    .map(H3Cell::as_u64),
            );
            f.opt_u64(child.child_position(parent_res));
        }

        assert!(ours.child_at(children.len() as u64, target_res).is_none());
        if parent_res > 0 {
            assert!(ours.child_at(0, parent_res - 1).is_none());
        }
        assert_eq!(ours.child_position(parent_res), Some(0));
        assert_eq!(ours.child_position(target_res.saturating_add(1)), None);
    }
    assert_eq!(f.finish(), golden::CHILD_POS, "child_position golden");
}

#[test]
fn compact_and_uncompact_round_trip() {
    let parent = H3Cell::for_point(40.758, -73.985, 8).expect("parent");
    let mut cells = parent.children_at_resolution(10);

    H3Cell::compact_cells(&mut cells).expect("compact ours");

    // The complete set of r10 descendants of one r8 cell compacts back to
    // exactly that single r8 parent.
    assert_eq!(cells, vec![parent]);

    let uncompacted = H3Cell::uncompact_cells(cells.iter().copied(), 10).expect("uncompact");
    let expected = parent.children_at_resolution(10);
    assert_eq!(uncompacted, expected);
    assert_eq!(
        H3Cell::uncompact_size(cells.iter().copied(), 10),
        Some(expected.len() as u64)
    );

    let mut duplicate = vec![parent, parent];
    assert_eq!(
        H3Cell::compact_cells(&mut duplicate),
        Err(H3CompactionError::DuplicateInput)
    );

    let mut heterogeneous = vec![parent, parent.center_child(9).expect("child")];
    assert_eq!(
        H3Cell::compact_cells(&mut heterogeneous),
        Err(H3CompactionError::HeterogeneousResolution)
    );

    assert!(H3Cell::uncompact_cells([parent], 7).is_none());
    assert_eq!(H3Cell::uncompact_size([parent], 16), None);
}

#[test]
fn pentagons_returns_twelve_per_resolution() {
    for res in [0, 1, 2, 5, 10, 15] {
        let pentagons = H3Cell::pentagons(res);
        assert_eq!(pentagons.len(), 12, "res {res}");
        for cell in pentagons {
            assert_eq!(cell.resolution(), res);
            assert!(cell.is_pentagon());
        }
    }
    assert!(H3Cell::pentagons(16).is_empty());
}
