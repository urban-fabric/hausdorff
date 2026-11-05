use super::*;

#[test]
fn grid_disk_matches_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");
        for k in 0..=4u32 {
            let mut m: Vec<u64> = ours.grid_disk(k).into_iter().map(|c| c.as_u64()).collect();
            // Metamorphic: the origin is always in its own disk, the disk grows
            // monotonically with k, and entries are unique.
            assert!(m.contains(&raw), "grid_disk {raw:#x} k={k} missing origin");
            m.sort_unstable();
            let dedup_len = {
                let mut d = m.clone();
                d.dedup();
                d.len()
            };
            assert_eq!(
                dedup_len,
                m.len(),
                "grid_disk {raw:#x} k={k} has duplicates"
            );
            for x in m {
                f.u64(x);
            }
        }
    }
    assert_eq!(f.finish(), golden::GRID_DISK, "grid_disk golden");
}

#[test]
fn grid_disk_distances_and_ring_are_consistent() {
    let cell = H3Cell::for_point(40.758, -73.985, 9).expect("cell");
    for k in 0..=4 {
        let disk = cell.grid_disk_distances(k);
        assert_eq!(disk.len() as u64, H3Cell::max_grid_disk_size(k));
        assert!(
            disk.iter()
                .any(|(candidate, distance)| { *candidate == cell && *distance == 0 })
        );
        for (candidate, distance) in &disk {
            assert!(*distance <= k);
            assert_eq!(
                cell.grid_distance(*candidate),
                Some(*distance as i32),
                "{:#x} -> {:#x}",
                cell.as_u64(),
                candidate.as_u64()
            );
        }

        let ring = cell.grid_ring(k);
        for candidate in &ring {
            assert_eq!(cell.grid_distance(*candidate), Some(k as i32));
        }
    }
}

#[test]
fn neighbor_check_matches_grid_disk_one() {
    let cell = H3Cell::for_point(40.758, -73.985, 9).expect("cell");
    for neighbor in cell.grid_ring(1) {
        assert!(cell.is_neighbor_with(neighbor));
        assert!(neighbor.is_neighbor_with(cell));
    }
    assert!(!cell.is_neighbor_with(cell));
}
