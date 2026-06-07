//! Differential tests for the in-house H3 codec against the h3o oracle.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::str::FromStr;

use hausdorff::index::h3::{
    H3Cell, H3CompactionError, H3DirectedEdge, H3Line, H3LineError, H3Polygon, H3PolygonError,
    H3Vertex,
};

/// A spread of real cells: every base cell (res 0), plus cells from a lat/lon
/// grid at several resolutions (including r10, the perf-test source resolution).
fn sample_cells() -> Vec<u64> {
    let mut cells: Vec<u64> = h3o::CellIndex::base_cells().map(u64::from).collect();
    let mut lat = -85i32;
    while lat <= 85 {
        let mut lon = -175i32;
        while lon <= 175 {
            for res in [0u8, 1, 2, 5, 7, 9, 10, 12, 15] {
                let r = h3o::Resolution::try_from(res).expect("resolution");
                let ll = h3o::LatLng::new(f64::from(lat), f64::from(lon)).expect("latlng");
                cells.push(u64::from(ll.to_cell(r)));
            }
            lon += 11;
        }
        lat += 7;
    }
    cells.sort_unstable();
    cells.dedup();
    cells
}

#[test]
fn codec_accessors_match_h3o() {
    for raw in sample_cells() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle valid");
        let ours = H3Cell::try_new(raw).expect("ours valid");

        assert_eq!(ours.as_u64(), raw);
        assert_eq!(
            ours.resolution(),
            u8::from(oracle.resolution()),
            "resolution for {raw:#x}"
        );
        assert_eq!(
            ours.base_cell(),
            u8::from(oracle.base_cell()),
            "base_cell for {raw:#x}"
        );
        assert_eq!(
            ours.is_pentagon(),
            oracle.is_pentagon(),
            "is_pentagon for {raw:#x}"
        );
        assert_eq!(ours.to_string(), oracle.to_string(), "token for {raw:#x}");
        assert_eq!(
            H3Cell::from_str(&oracle.to_string())
                .expect("parse")
                .as_u64(),
            raw
        );
    }
}

#[test]
fn raw_parts_reconstruct_sample_cells() {
    for raw in sample_cells() {
        let cell = H3Cell::try_new(raw).expect("valid sample");
        let digits = (1..=cell.resolution())
            .map(|res| cell.digit(res).expect("digit"))
            .collect::<Vec<_>>();
        let rebuilt =
            H3Cell::from_raw_parts(cell.resolution(), cell.base_cell(), &digits).expect("rebuilt");
        assert_eq!(rebuilt, cell, "rebuild {raw:#x}");
        assert!(H3Cell::is_valid(raw));
    }

    assert!(H3Cell::from_raw_parts(16, 0, &[]).is_none());
    assert!(H3Cell::from_raw_parts(1, 122, &[0]).is_none());
    assert!(H3Cell::from_raw_parts(1, 0, &[7]).is_none());
}

#[test]
fn resolution_metadata_matches_h3_reference_constants() {
    assert_eq!(H3Cell::num_cells(0), Some(122));
    assert_eq!(H3Cell::num_cells(1), Some(842));
    assert_eq!(H3Cell::num_cells(15), Some(2 + 120 * 7u64.pow(15)));
    assert_eq!(H3Cell::num_cells(16), None);

    assert_eq!(H3Cell::max_grid_disk_size(0), 1);
    assert_eq!(H3Cell::max_grid_disk_size(1), 7);
    assert_eq!(H3Cell::max_grid_disk_size(4), 61);

    assert_eq!(H3Cell::average_hexagon_edge_length_m(9), Some(200.7861476));
    assert_eq!(
        H3Cell::average_hexagon_area_m2(15),
        Some(8.953115907605802e-01)
    );
    assert_eq!(H3Cell::average_hexagon_area_km2(16), None);
    assert_eq!(H3Cell::average_hexagon_edge_length_km(16), None);

    assert_eq!(
        H3Cell::average_hexagon_edge_length(9).map(|l| l.as_meters()),
        H3Cell::average_hexagon_edge_length_m(9),
    );
    assert_eq!(
        H3Cell::average_hexagon_area(15).map(|a| a.as_square_meters()),
        H3Cell::average_hexagon_area_m2(15),
    );
    assert_eq!(H3Cell::average_hexagon_area(16), None);
}

#[test]
fn validity_matches_h3o_on_random_bits() {
    let mut state = 0x1234_5678_9abc_def0u64;
    let mut valid_seen = 0u32;
    for _ in 0..500_000 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let raw = state;
        let ours = H3Cell::try_new(raw).is_some();
        let oracle = h3o::CellIndex::try_from(raw).is_ok();
        assert_eq!(ours, oracle, "validity mismatch for {raw:#x}");
        if oracle {
            valid_seen += 1;
        }
    }
    let _ = valid_seen;
}

#[test]
fn center_matches_h3o() {
    for raw in sample_cells() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ll: h3o::LatLng = oracle.into();
        let (lat, lon) = H3Cell::try_new(raw).expect("ours").center_lat_lng();
        let dlat = (lat - ll.lat()).abs();
        let dlon = (lon - ll.lng()).abs();
        assert!(
            dlat < 1e-9 && dlon < 1e-9,
            "center {raw:#x}: ours=({lat}, {lon}) h3o=({}, {}) d=({dlat}, {dlon})",
            ll.lat(),
            ll.lng()
        );
    }
}

#[test]
fn for_point_rejects_out_of_range_and_nonfinite() {
    assert!(H3Cell::for_point(52.52, 13.405, 9).is_some());
    assert!(H3Cell::for_point(90.001, 0.0, 9).is_none());
    assert!(H3Cell::for_point(-90.001, 0.0, 9).is_none());
    assert!(H3Cell::for_point(0.0, 180.001, 9).is_none());
    assert!(H3Cell::for_point(0.0, -180.001, 9).is_none());
    assert!(H3Cell::for_point(f64::NAN, 0.0, 9).is_none());
    assert!(H3Cell::for_point(0.0, f64::INFINITY, 9).is_none());
    assert!(H3Cell::for_point(0.0, 0.0, 16).is_none());
    assert!(H3Cell::for_point(90.0, 180.0, 0).is_some());
    assert!(H3Cell::for_point(-90.0, -180.0, 0).is_some());
}

#[test]
fn for_point_matches_h3o() {
    let mut lat = -85i32;
    while lat <= 85 {
        let mut lon = -179i32;
        while lon <= 179 {
            for res in [0u8, 1, 2, 5, 7, 9, 10, 12, 15] {
                let r = h3o::Resolution::try_from(res).expect("res");
                let oracle = h3o::LatLng::new(f64::from(lat), f64::from(lon))
                    .expect("latlng")
                    .to_cell(r);
                let ours = H3Cell::for_point(f64::from(lat), f64::from(lon), res).expect("ours");
                assert_eq!(
                    ours.as_u64(),
                    u64::from(oracle),
                    "for_point lat={lat} lon={lon} res={res}: ours={:#x} h3o={:#x}",
                    ours.as_u64(),
                    u64::from(oracle)
                );
            }
            lon += 13;
        }
        lat += 6;
    }
}

/// A smaller spread for the O(k^2)/O(7^d) ops: every base cell + a coarse
/// lat/lon grid at a few resolutions (includes pentagons).
fn small_sample() -> Vec<u64> {
    let mut v: Vec<u64> = h3o::CellIndex::base_cells().map(u64::from).collect();
    let mut lat = -80i32;
    while lat <= 80 {
        let mut lon = -170i32;
        while lon <= 170 {
            for res in [5u8, 7, 9, 10] {
                let r = h3o::Resolution::try_from(res).expect("res");
                let ll = h3o::LatLng::new(f64::from(lat), f64::from(lon)).expect("latlng");
                v.push(u64::from(ll.to_cell(r)));
            }
            lon += 37;
        }
        lat += 23;
    }
    v.sort_unstable();
    v.dedup();
    v
}

#[test]
fn parent_matches_h3o() {
    for raw in sample_cells() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");
        let res = u8::from(oracle.resolution());
        for pres in 0..=res {
            let or_p = oracle
                .parent(h3o::Resolution::try_from(pres).unwrap())
                .expect("oracle parent");
            let ou_p = ours.parent(pres).expect("our parent");
            assert_eq!(ou_p.as_u64(), u64::from(or_p), "parent {raw:#x} -> r{pres}");
        }
        assert!(ours.parent(res + 1).is_none());
    }
}

#[test]
fn children_match_h3o() {
    for raw in small_sample() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");
        let pres = u8::from(oracle.resolution());
        let tres = (pres + 2).min(15);
        let o: Vec<u64> = oracle
            .children(h3o::Resolution::try_from(tres).unwrap())
            .map(u64::from)
            .collect();
        let m: Vec<u64> = ours
            .children_at_resolution(tres)
            .into_iter()
            .map(|c| c.as_u64())
            .collect();
        assert_eq!(m, o, "children {raw:#x} -> r{tres}");
        assert_eq!(ours.children_count(tres), o.len() as u64);
        assert_eq!(
            ours.center_child(tres).map(H3Cell::as_u64),
            o.first().copied(),
            "center child {raw:#x} -> r{tres}"
        );
    }
}

#[test]
fn child_position_and_child_at_match_h3o() {
    for raw in small_sample() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");
        let parent_res = ours.resolution();
        let target_res = (parent_res + 2).min(15);
        let target_resolution = h3o::Resolution::try_from(target_res).expect("target res");
        let parent_resolution = h3o::Resolution::try_from(parent_res).expect("parent res");

        let children = ours.children_at_resolution(target_res);
        for (position, child) in children.iter().copied().enumerate() {
            let oracle_child = oracle
                .child_at(position as u64, target_resolution)
                .expect("oracle child");
            assert_eq!(
                child.as_u64(),
                u64::from(oracle_child),
                "child_at {raw:#x} pos={position} r{target_res}"
            );
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
            assert_eq!(
                child.child_position(parent_res),
                oracle_child.child_position(parent_resolution),
                "oracle child_position {raw:#x} child={:#x}",
                child.as_u64()
            );
        }

        assert!(ours.child_at(children.len() as u64, target_res).is_none());
        if parent_res > 0 {
            assert!(ours.child_at(0, parent_res - 1).is_none());
        }
        assert_eq!(ours.child_position(parent_res), Some(0));
        assert_eq!(ours.child_position(target_res.saturating_add(1)), None);
    }
}

#[test]
fn compact_and_uncompact_match_h3o() {
    let parent = H3Cell::for_point(40.758, -73.985, 8).expect("parent");
    let mut cells = parent.children_at_resolution(10);
    let mut oracle = cells
        .iter()
        .map(|cell| h3o::CellIndex::try_from(cell.as_u64()).expect("oracle cell"))
        .collect::<Vec<_>>();

    H3Cell::compact_cells(&mut cells).expect("compact ours");
    h3o::CellIndex::compact(&mut oracle).expect("compact oracle");

    let oracle_raw = oracle.into_iter().map(u64::from).collect::<Vec<_>>();
    let ours_raw = cells.iter().map(|cell| cell.as_u64()).collect::<Vec<_>>();
    assert_eq!(ours_raw, oracle_raw);
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

#[test]
fn grid_disk_matches_h3o() {
    for raw in small_sample() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");
        for k in 0..=4u32 {
            let mut o: Vec<u64> = oracle
                .grid_disk::<Vec<_>>(k)
                .into_iter()
                .map(u64::from)
                .collect();
            let mut m: Vec<u64> = ours.grid_disk(k).into_iter().map(|c| c.as_u64()).collect();
            o.sort_unstable();
            m.sort_unstable();
            assert_eq!(m, o, "grid_disk {raw:#x} k={k}");
        }
    }
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

#[test]
fn directed_edges_match_h3o() {
    for raw in small_sample() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");

        let oracle_edges = oracle.edges().map(u64::from).collect::<Vec<_>>();
        let ours_edges = ours
            .edges()
            .into_iter()
            .map(H3DirectedEdge::as_u64)
            .collect::<Vec<_>>();
        assert_eq!(ours_edges, oracle_edges, "origin edges {raw:#x}");
        assert_eq!(ours_edges.len(), if ours.is_pentagon() { 5 } else { 6 });

        for edge in ours.edges() {
            let oracle_edge =
                h3o::DirectedEdgeIndex::try_from(edge.as_u64()).expect("oracle directed edge");
            let destination = edge.destination();
            let (origin, destination_from_cells) = edge.cells();

            assert!(H3DirectedEdge::is_valid(edge.as_u64()));
            assert_eq!(edge.to_string(), oracle_edge.to_string());
            assert_eq!(H3DirectedEdge::from_str(&edge.to_string()).ok(), Some(edge));
            assert_eq!(origin, ours);
            assert_eq!(destination_from_cells, destination);
            assert_eq!(edge.origin().as_u64(), u64::from(oracle_edge.origin()));
            assert_eq!(
                destination.as_u64(),
                u64::from(oracle_edge.destination()),
                "edge destination {:#x}",
                edge.as_u64()
            );
            assert_eq!(
                edge.reverse().as_u64(),
                u64::from(oracle_edge.reverse()),
                "reverse edge {:#x}",
                edge.as_u64()
            );
            assert_eq!(ours.edge(destination), Some(edge));
            assert_eq!(H3DirectedEdge::new(ours, destination), Some(edge));
        }

        for neighbor in ours.grid_ring(1) {
            let oracle_neighbor =
                h3o::CellIndex::try_from(neighbor.as_u64()).expect("oracle neighbor");
            assert_eq!(
                ours.edge(neighbor).map(H3DirectedEdge::as_u64),
                oracle.edge(oracle_neighbor).map(u64::from),
                "edge {raw:#x} -> {:#x}",
                neighbor.as_u64()
            );
        }
    }
}

#[test]
fn directed_edge_geometry_matches_h3o() {
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");
        for edge in ours.edges() {
            let oracle = h3o::DirectedEdgeIndex::try_from(edge.as_u64()).expect("oracle edge");
            let oracle_boundary = oracle.boundary();
            let ours_boundary = edge.boundary_lat_lng();
            assert_eq!(
                ours_boundary.len(),
                oracle_boundary.len(),
                "edge boundary length {:#x}",
                edge.as_u64()
            );
            for (ours_point, oracle_point) in ours_boundary.iter().zip(oracle_boundary.iter()) {
                assert!(
                    (ours_point.0 - oracle_point.lat()).abs() < 1e-9
                        && (ours_point.1 - oracle_point.lng()).abs() < 1e-9,
                    "edge boundary {:#x}: ours={ours_point:?} oracle=({}, {})",
                    edge.as_u64(),
                    oracle_point.lat(),
                    oracle_point.lng()
                );
            }
            assert!(
                (edge.length_rads() - oracle.length_rads()).abs() < 1e-15,
                "edge length rads {:#x}",
                edge.as_u64()
            );
            assert!(
                (edge.length_km() - oracle.length_km()).abs() < 1e-10,
                "edge length km {:#x}",
                edge.as_u64()
            );
            assert!(
                (edge.length_m() - oracle.length_m()).abs() < 1e-7,
                "edge length m {:#x}",
                edge.as_u64()
            );
            assert_eq!(edge.length().as_meters(), edge.length_m());
        }
    }
}

#[test]
fn directed_edges_reject_invalid_inputs() {
    let cell = H3Cell::for_point(40.758, -73.985, 9).expect("cell");
    let far = H3Cell::for_point(34.0522, -118.2437, 9).expect("far");
    let pentagon = H3Cell::pentagons(1).pop().expect("pentagon");

    assert!(H3DirectedEdge::try_new(0).is_none());
    assert!(!H3DirectedEdge::is_valid(cell.as_u64()));
    assert!(H3DirectedEdge::new(cell, cell).is_none());
    assert!(H3DirectedEdge::new(cell, far).is_none());
    assert!(H3DirectedEdge::from_origin_and_direction(cell, 0).is_none());
    assert!(H3DirectedEdge::from_origin_and_direction(cell, 7).is_none());
    assert!(H3DirectedEdge::from_origin_and_direction(pentagon, 1).is_none());
}

#[test]
fn vertices_match_h3o() {
    for raw in small_sample() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");

        let oracle_vertices = oracle.vertexes().map(u64::from).collect::<Vec<_>>();
        let ours_vertices = ours
            .vertexes()
            .into_iter()
            .map(H3Vertex::as_u64)
            .collect::<Vec<_>>();
        assert_eq!(ours_vertices, oracle_vertices, "cell vertexes {raw:#x}");
        assert_eq!(
            ours.vertex_count(),
            if ours.is_pentagon() { 5 } else { 6 },
            "vertex count {raw:#x}"
        );

        for vertex in 0..ours.vertex_count() {
            let h3o_vertex = h3o::Vertex::try_from(vertex).expect("h3o vertex");
            let oracle_vertex = oracle.vertex(h3o_vertex).expect("oracle vertex");
            let ours_vertex = ours.vertex(vertex).expect("our vertex");
            let owner = ours_vertex.owner();
            let ll: h3o::LatLng = oracle_vertex.into();
            let ours_ll = ours_vertex.lat_lng();

            assert_eq!(
                ours_vertex.as_u64(),
                u64::from(oracle_vertex),
                "vertex {raw:#x} v{vertex}"
            );
            assert!(H3Vertex::is_valid(ours_vertex.as_u64()));
            assert_eq!(H3Vertex::try_new(ours_vertex.as_u64()), Some(ours_vertex));
            assert_eq!(
                H3Vertex::from_str(&ours_vertex.to_string()).ok(),
                Some(ours_vertex)
            );
            assert_eq!(owner.as_u64(), u64::from(oracle_vertex.owner()));
            assert_eq!(ours_vertex.vertex(), u8::from(oracle_vertex.vertex()));
            assert!(
                (ours_ll.0 - ll.lat()).abs() < 1e-9 && (ours_ll.1 - ll.lng()).abs() < 1e-9,
                "vertex latlng {raw:#x} v{vertex}: ours={ours_ll:?} oracle=({}, {})",
                ll.lat(),
                ll.lng()
            );
        }

        assert!(ours.vertex(ours.vertex_count()).is_none());
    }
}

#[test]
fn vertices_reject_invalid_inputs() {
    let cell = H3Cell::for_point(40.758, -73.985, 9).expect("cell");
    let pentagon = H3Cell::pentagons(1).pop().expect("pentagon");

    assert!(H3Vertex::try_new(0).is_none());
    assert!(!H3Vertex::is_valid(cell.as_u64()));
    assert!(H3Vertex::new(cell, 6).is_none());
    assert!(H3Vertex::new(pentagon, 5).is_none());
}

#[test]
fn grid_distance_matches_h3o() {
    let anchors = small_sample();
    for &raw in anchors.iter().take(400) {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");
        for tgt in ours.grid_disk(3) {
            let o_tgt = h3o::CellIndex::try_from(tgt.as_u64()).unwrap();
            let want = oracle.grid_distance(o_tgt).ok();
            let got = ours.grid_distance(tgt);
            assert_eq!(got, want, "grid_distance {raw:#x} -> {:#x}", tgt.as_u64());
        }
    }
}

#[test]
fn local_ij_matches_h3o() {
    for raw in small_sample().into_iter().take(300) {
        let oracle_origin = h3o::CellIndex::try_from(raw).expect("oracle origin");
        let origin = H3Cell::try_new(raw).expect("origin");
        for cell in origin.grid_disk(2) {
            let oracle_cell = h3o::CellIndex::try_from(cell.as_u64()).expect("oracle cell");
            let oracle_ij = oracle_cell.to_local_ij(oracle_origin).ok();
            let ours_ij = cell.to_local_ij(origin);
            assert_eq!(
                ours_ij.map(|ij| (ij.anchor.as_u64(), ij.i, ij.j)),
                oracle_ij.map(|ij| (u64::from(ij.anchor), ij.coord.i, ij.coord.j)),
                "local IJ {raw:#x} -> {:#x}",
                cell.as_u64()
            );
        }
    }
}

#[test]
fn grid_path_cells_are_shortest_neighbor_paths() {
    for raw in small_sample().into_iter().take(200) {
        let oracle_origin = h3o::CellIndex::try_from(raw).expect("oracle origin");
        let origin = H3Cell::try_new(raw).expect("origin");
        for destination in origin.grid_ring(3) {
            let oracle_destination =
                h3o::CellIndex::try_from(destination.as_u64()).expect("oracle destination");
            let oracle_size = oracle_origin.grid_path_cells_size(oracle_destination).ok();
            let path = origin.grid_path_cells(destination);

            assert_eq!(
                path.as_ref().map(|cells| cells.len() as i32),
                oracle_size,
                "path size {raw:#x} -> {:#x}",
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
}

#[test]
fn boundary_bbox_matches_h3o() {
    for raw in sample_cells() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");

        let (mut o_min_lat, mut o_min_lon, mut o_max_lat, mut o_max_lon) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for v in oracle.boundary().iter() {
            o_min_lat = o_min_lat.min(v.lat());
            o_min_lon = o_min_lon.min(v.lng());
            o_max_lat = o_max_lat.max(v.lat());
            o_max_lon = o_max_lon.max(v.lng());
        }

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

        for (a, b, what) in [
            (m_min_lat, o_min_lat, "min_lat"),
            (m_min_lon, o_min_lon, "min_lon"),
            (m_max_lat, o_max_lat, "max_lat"),
            (m_max_lon, o_max_lon, "max_lon"),
        ] {
            assert!(
                (a - b).abs() < 1e-7,
                "bbox {what} for {raw:#x}: ours={a} h3o={b}"
            );
        }
    }
}

#[test]
fn exact_cell_area_matches_h3o() {
    for raw in sample_cells() {
        let oracle = h3o::CellIndex::try_from(raw).expect("oracle");
        let ours = H3Cell::try_new(raw).expect("ours");
        assert!(
            (ours.area_rads2() - oracle.area_rads2()).abs() < 1e-15,
            "area rads2 {raw:#x}: ours={} oracle={}",
            ours.area_rads2(),
            oracle.area_rads2()
        );
        assert!(
            (ours.area_km2() - oracle.area_km2()).abs() < 1e-8,
            "area km2 {raw:#x}: ours={} oracle={}",
            ours.area_km2(),
            oracle.area_km2()
        );
        assert!(
            (ours.area_m2() - oracle.area_m2()).abs() < 1e-2,
            "area m2 {raw:#x}: ours={} oracle={}",
            ours.area_m2(),
            oracle.area_m2()
        );
        assert_eq!(ours.area().as_square_meters(), ours.area_m2());
    }
}

#[test]
fn polygon_contains_points_and_holes() {
    let polygon = H3Polygon::new(
        vec![
            (40.70, -74.05),
            (40.70, -73.90),
            (40.85, -73.90),
            (40.85, -74.05),
            (40.70, -74.05),
        ],
        vec![vec![
            (40.74, -74.01),
            (40.74, -73.97),
            (40.78, -73.97),
            (40.78, -74.01),
        ]],
    )
    .expect("polygon");

    assert!(polygon.contains_lat_lng(40.82, -73.98));
    assert!(!polygon.contains_lat_lng(40.76, -73.99));
    assert!(!polygon.contains_lat_lng(41.0, -73.98));
}

#[test]
fn polygon_to_cells_uses_centroid_containment() {
    let polygon = H3Polygon::new(
        vec![
            (40.70, -74.05),
            (40.70, -73.90),
            (40.85, -73.90),
            (40.85, -74.05),
        ],
        Vec::<Vec<(f64, f64)>>::new(),
    )
    .expect("polygon");

    let cells = polygon.to_cells(7).expect("cells");
    assert!(!cells.is_empty());
    assert!(cells.windows(2).all(|pair| pair[0] < pair[1]));
    for cell in &cells {
        let (lat, lng) = cell.center_lat_lng();
        assert!(
            polygon.contains_lat_lng(lat, lng),
            "cell center outside polygon: {cell} ({lat}, {lng})"
        );
    }

    assert_eq!(
        polygon.to_cells(16),
        Err(H3PolygonError::InvalidResolution(16))
    );
}

#[test]
fn cells_to_boundary_rings_dissolves_internal_edges() {
    let cell = H3Cell::for_point(40.758, -73.985, 8).expect("cell");
    let rings = H3Cell::cells_to_boundary_rings([cell]).expect("single ring");
    assert_eq!(rings.len(), 1);
    assert_eq!(rings[0].points().len(), cell.boundary_lat_lng().len());
    assert!((rings[0].area_rads2() - cell.area_rads2()).abs() < 1e-15);

    let neighbor = cell.grid_ring(1).into_iter().next().expect("neighbor");
    let rings = H3Cell::cells_to_boundary_rings([cell, neighbor]).expect("merged ring");
    assert_eq!(rings.len(), 1);
    assert!(rings[0].points().len() < cell.boundary_lat_lng().len() * 2);
    let area = cell.area_rads2() + neighbor.area_rads2();
    assert!((rings[0].area_rads2() - area).abs() < 1e-12);

    assert_eq!(
        H3Cell::cells_to_boundary_rings([cell, cell]),
        Err(H3PolygonError::DuplicateInput)
    );
    assert_eq!(
        H3Cell::cells_to_boundary_rings([cell, cell.center_child(9).expect("child")]),
        Err(H3PolygonError::HeterogeneousResolution)
    );
}

#[test]
fn polygon_rejects_invalid_rings() {
    assert_eq!(
        H3Polygon::new(vec![(0.0, 0.0), (1.0, 1.0)], Vec::<Vec<(f64, f64)>>::new()),
        Err(H3PolygonError::RingTooShort)
    );
    assert_eq!(
        H3Polygon::new(
            vec![(0.0, -179.0), (0.0, 179.0), (1.0, 179.0)],
            Vec::<Vec<(f64, f64)>>::new()
        ),
        Err(H3PolygonError::TransmeridianUnsupported)
    );
    assert_eq!(
        H3Polygon::new(
            vec![(0.0, 200.0), (0.0, 201.0), (1.0, 200.5)],
            Vec::<Vec<(f64, f64)>>::new()
        ),
        Err(H3PolygonError::InvalidCoordinate)
    );
}

#[test]
fn rejects_non_cell_indices() {
    assert!(H3Cell::try_new(0).is_none());
    let valid = h3o::LatLng::new(40.7, -74.0)
        .unwrap()
        .to_cell(h3o::Resolution::Nine);
    let raw = u64::from(valid);
    let edge_mode = (raw & !(0xf << 59)) | (2 << 59);
    assert!(H3Cell::try_new(edge_mode).is_none());
    assert!(H3Cell::try_new(raw).is_some());
}

#[test]
fn line_cover_endpoints_are_vertex_cells_and_path_is_connected() {
    let resolution = 9;
    let vertices = vec![(40.70, -74.00), (40.70, -73.95), (40.72, -73.90)];
    let line = H3Line::new(vertices.clone()).expect("line");
    let cells = line.to_cells(resolution).expect("cells");

    let first = H3Cell::for_point(vertices[0].0, vertices[0].1, resolution).expect("first cell");
    let &(last_lat, last_lng) = vertices.last().expect("last vertex");
    let last = H3Cell::for_point(last_lat, last_lng, resolution).expect("last cell");
    assert_eq!(cells.first().copied(), Some(first));
    assert_eq!(cells.last().copied(), Some(last));

    for &(lat, lng) in &vertices {
        let cell = H3Cell::for_point(lat, lng, resolution).expect("vertex cell");
        assert!(
            cells.contains(&cell),
            "vertex cell {cell} missing from cover"
        );
    }

    let mut unique = cells.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(unique.len(), cells.len(), "cover has duplicate cells");
    for pair in cells.windows(2) {
        assert!(
            pair[0].is_neighbor_with(pair[1]),
            "cover step {} -> {} is not a grid neighbour",
            pair[0],
            pair[1]
        );
    }

    assert_eq!(H3Cell::line_to_cells(&line, resolution), Ok(cells));
}

#[test]
fn line_cover_within_one_cell_is_that_cell() {
    let line = H3Line::new(vec![(40.700, -74.000), (40.701, -74.001)]).expect("line");
    let cell = H3Cell::for_point(40.700, -74.000, 1).expect("cell");
    assert_eq!(line.to_cells(1), Ok(vec![cell]));

    let dot = H3Line::new(vec![(51.5, -0.12)]).expect("dot");
    let dot_cell = H3Cell::for_point(51.5, -0.12, 9).expect("dot cell");
    assert_eq!(dot.to_cells(9), Ok(vec![dot_cell]));
}

#[test]
fn line_cover_matches_concatenated_segment_grid_paths() {
    let resolution = 7;
    let vertices = vec![(40.00, -74.00), (40.00, -73.80), (40.15, -73.65)];
    let line = H3Line::new(vertices.clone()).expect("line");
    let cells = line.to_cells(resolution).expect("cells");

    let vertex_cells: Vec<H3Cell> = vertices
        .iter()
        .map(|&(lat, lng)| H3Cell::for_point(lat, lng, resolution).expect("vertex cell"))
        .collect();
    let mut seen = std::collections::BTreeSet::new();
    let mut expected = Vec::new();
    if seen.insert(vertex_cells[0]) {
        expected.push(vertex_cells[0]);
    }
    for pair in vertex_cells.windows(2) {
        if pair[0] == pair[1] {
            continue;
        }
        for cell in pair[0].grid_path_cells(pair[1]).expect("grid path") {
            if seen.insert(cell) {
                expected.push(cell);
            }
        }
    }
    assert_eq!(cells, expected);
}

#[test]
fn line_cover_validates_input() {
    assert_eq!(
        H3Line::new(Vec::<(f64, f64)>::new()),
        Err(H3LineError::Empty)
    );
    assert_eq!(
        H3Line::new(vec![(40.0, -74.0), (f64::NAN, 0.0)]),
        Err(H3LineError::InvalidCoordinate)
    );
    assert_eq!(
        H3Line::new(vec![(91.0, 0.0)]),
        Err(H3LineError::InvalidCoordinate)
    );
    let line = H3Line::new(vec![(40.0, -74.0), (41.0, -73.0)]).expect("line");
    assert_eq!(line.to_cells(16), Err(H3LineError::InvalidResolution(16)));
}

#[test]
fn for_points_batch_matches_scalar() {
    let lats = [0.0, 12.3, -45.6, 80.0, 91.0, f64::NAN, 51.5];
    let lngs = [0.0, 34.5, -120.0, 179.9, 10.0, 0.0, -0.12];
    for res in [0u8, 5, 9, 15] {
        let batch = H3Cell::for_points(&lats, &lngs, res);
        let scalar: Vec<_> = lats
            .iter()
            .zip(&lngs)
            .map(|(&la, &lo)| H3Cell::for_point(la, lo, res))
            .collect();
        assert_eq!(batch, scalar, "batch != scalar at res {res}");
    }
    assert_eq!(H3Cell::for_points(&lats, &lngs, 16), vec![None; lats.len()]);
    assert_eq!(H3Cell::for_points(&lats, &lngs[..2], 9).len(), 2);
}

#[test]
fn for_points_approx_is_exact_or_edge_neighbour() {
    let mut lats = Vec::new();
    let mut lngs = Vec::new();
    let mut lat = -85.0;
    while lat <= 85.0 {
        let mut lng = -175.0;
        while lng <= 175.0 {
            lats.push(lat);
            lngs.push(lng);
            lng += 7.3;
        }
        lat += 5.1;
    }

    for res in [0u8, 5, 9, 12, 15] {
        let approx = H3Cell::for_points_approx(&lats, &lngs, res);
        let exact = H3Cell::for_points(&lats, &lngs, res);
        for ((a, e), (&la, &lo)) in approx.iter().zip(&exact).zip(lats.iter().zip(&lngs)) {
            match (a, e) {
                (Some(a), Some(e)) => assert!(
                    a == e || a.is_neighbor_with(*e),
                    "approx {a} vs exact {e} at ({lo},{la}) res {res}: not equal or adjacent"
                ),
                (None, None) => {}
                _ => panic!("approx/exact disagree on validity at ({lo},{la}) res {res}"),
            }
        }
    }

    assert_eq!(
        H3Cell::for_points_approx(&lats, &lngs, 16),
        vec![None; lats.len()]
    );
}
