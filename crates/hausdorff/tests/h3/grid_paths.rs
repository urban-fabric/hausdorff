use super::*;

#[test]
fn grid_distance_matches_golden() {
    let mut f = Fnv::new();
    let anchors = small_sample();
    for &raw in anchors.iter().take(400) {
        let ours = H3Cell::try_new(raw).expect("ours");
        for tgt in ours.grid_disk(3) {
            let got = ours.grid_distance(tgt);
            // Metamorphic: anything inside grid_disk(3) is at distance <= 3, and
            // the distance to self is 0.
            if let Some(d) = got {
                assert!(
                    (0..=3).contains(&d),
                    "grid_distance {raw:#x} out of disk: {d}"
                );
            }
            assert_eq!(ours.grid_distance(ours), Some(0));
            f.opt_u64(got.map(|x| x as u32 as u64));
        }
    }
    assert_eq!(f.finish(), golden::GRID_DISTANCE, "grid_distance golden");
}

#[test]
fn local_ij_matches_golden() {
    let mut f = Fnv::new();
    for raw in small_sample().into_iter().take(300) {
        let origin = H3Cell::try_new(raw).expect("origin");
        for cell in origin.grid_disk(2) {
            match cell.to_local_ij(origin) {
                Some(ij) => {
                    assert_eq!(ij.anchor, origin, "local IJ anchor {raw:#x}");
                    f.byte(1);
                    f.u64(ij.anchor.as_u64());
                    f.i32(ij.i);
                    f.i32(ij.j);
                }
                None => f.byte(0),
            }
        }
    }
    assert_eq!(f.finish(), golden::LOCAL_IJ, "local IJ golden");
}

#[test]
fn grid_path_cells_are_shortest_neighbor_paths() {
    let mut f = Fnv::new();
    for raw in small_sample().into_iter().take(200) {
        let origin = H3Cell::try_new(raw).expect("origin");
        for destination in origin.grid_ring(3) {
            let path = origin.grid_path_cells(destination);

            // Pin the path size to the golden, and assert it agrees with the
            // standalone size accessor (self-consistency).
            f.opt_u64(path.as_ref().map(|cells| cells.len() as u64));
            assert_eq!(
                origin.grid_path_cells_size(destination).map(|s| s as usize),
                path.as_ref().map(Vec::len),
                "path size accessor {raw:#x} -> {:#x}",
                destination.as_u64()
            );

            let Some(path) = path else {
                continue;
            };
            assert_eq!(path.first().copied(), Some(origin));
            assert_eq!(path.last().copied(), Some(destination));
            assert_eq!(
                origin.grid_path_cells_size(destination),
                Some(path.len() as u32)
            );
            for pair in path.windows(2) {
                assert!(
                    pair[0].is_neighbor_with(pair[1]),
                    "path neighbor step {:#x} -> {:#x}",
                    pair[0].as_u64(),
                    pair[1].as_u64()
                );
            }
        }
    }
    assert_eq!(f.finish(), golden::GRID_PATH, "grid_path size golden");
}

#[test]
fn boundary_bbox_matches_golden() {
    let mut f = Fnv::new();
    for raw in sample_cells() {
        let ours = H3Cell::try_new(raw).expect("ours");

        let (mut m_min_lat, mut m_min_lon, mut m_max_lat, mut m_max_lon) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for (lat, lon) in ours.boundary_lat_lng() {
            m_min_lat = m_min_lat.min(lat);
            m_min_lon = m_min_lon.min(lon);
            m_max_lat = m_max_lat.max(lat);
            m_max_lon = m_max_lon.max(lon);
        }

        // The bounding box is well-formed (min <= max) before pinning it.
        assert!(
            m_min_lat <= m_max_lat && m_min_lon <= m_max_lon,
            "bbox {raw:#x}"
        );
        f.f64q(m_min_lat);
        f.f64q(m_min_lon);
        f.f64q(m_max_lat);
        f.f64q(m_max_lon);
    }
    assert_eq!(f.finish(), golden::BBOX, "boundary bbox golden");
}

#[test]
fn exact_cell_area_matches_golden() {
    let mut f = Fnv::new();
    for raw in sample_cells() {
        let ours = H3Cell::try_new(raw).expect("ours");
        // Self-consistency among the area accessors.
        assert!(
            ours.area_rads2() > 0.0 && ours.area_m2() > 0.0,
            "area sign {raw:#x}"
        );
        assert_eq!(ours.area().value(), ours.area_m2());
        f.f64q(ours.area_rads2());
        f.f64q(ours.area_km2());
        f.f64q(ours.area_m2());
    }
    assert_eq!(f.finish(), golden::AREA, "cell area golden");
}
