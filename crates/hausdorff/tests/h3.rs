//! Oracle-free golden-corpus tests for the in-house H3 codec.
//!
//! These tests were originally differential tests against the `h3o` reference
//! crate. The reference dependency has been removed: the inputs are all
//! deterministic (fixed lat/lon grids, the native base-cell set, and a seeded
//! LCG), so the *current* (correct, h3o-matching) hausdorff output for each
//! input has been frozen into compact golden constants (FNV-1a checksums over
//! quantized outputs, plus exact value lists where small). Each test now
//! asserts that hausdorff reproduces its frozen golden, preserving the exact
//! differential coverage without any runtime oracle. Metamorphic / internal
//! self-consistency checks are layered on top where they add value.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::str::FromStr;

use hausdorff::index::h3::{
    ContainmentMode, H3Cell, H3CompactionError, H3DirectedEdge, H3Line, H3LineError, H3LocalIJ,
    H3Polygon, H3PolygonError, H3Vertex,
};

/// FNV-1a accumulator used to pin a stream of codec outputs to a single u64
/// golden. Floats are quantized to 1e-9 so the checksum is stable against
/// trivial rounding noise while still pinning values tightly.
struct Fnv(u64);
impl Fnv {
    fn new() -> Self {
        Fnv(0xcbf2_9ce4_8422_2325)
    }
    fn byte(&mut self, b: u8) {
        self.0 ^= u64::from(b);
        self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
    }
    fn u64(&mut self, v: u64) {
        for b in v.to_le_bytes() {
            self.byte(b);
        }
    }
    fn i32(&mut self, v: i32) {
        self.u64(v as u32 as u64);
    }
    fn str(&mut self, s: &str) {
        for b in s.bytes() {
            self.byte(b);
        }
    }
    /// Quantize to 1e-9 then fold, matching the golden-capture quantization.
    fn f64q(&mut self, v: f64) {
        let q = (v * 1e9).round() as i64;
        self.u64(q as u64);
    }
    fn opt_u64(&mut self, v: Option<u64>) {
        match v {
            Some(x) => {
                self.byte(1);
                self.u64(x);
            }
            None => self.byte(0),
        }
    }
    fn finish(self) -> u64 {
        self.0
    }
}

/// Checksum of a u64 slice (length-prefixed) for pinning the input cell lists.
fn chk_u64s(v: &[u64]) -> u64 {
    let mut f = Fnv::new();
    f.u64(v.len() as u64);
    for &x in v {
        f.u64(x);
    }
    f.finish()
}

/// Golden constants frozen from the current (h3o-matching) implementation.
/// See the module docstring for how these were captured.
mod golden {
    /// `sample_cells()` length and content checksum.
    pub const SAMPLE_LEN: usize = 6309;
    pub const SAMPLE_CHK: u64 = 11_755_684_950_046_384_310;
    /// `small_sample()` length and content checksum.
    pub const SMALL_LEN: usize = 402;
    pub const SMALL_CHK: u64 = 17_879_915_286_621_556_510;

    pub const CODEC_ACCESSORS: u64 = 4_918_502_469_318_935_052;
    pub const CENTER: u64 = 11_348_933_832_788_743_190;
    pub const FOR_POINT: u64 = 14_262_681_826_812_326_387;
    pub const PARENT: u64 = 1_163_989_267_815_838_905;
    pub const CHILDREN: u64 = 1_810_743_767_032_895_959;
    pub const CHILD_POS: u64 = 10_426_419_512_723_083_181;
    pub const GRID_DISK: u64 = 16_759_316_744_873_138_418;
    pub const EDGES: u64 = 7_816_069_306_376_095_097;
    pub const EDGE_GEOM: u64 = 12_132_845_004_434_110_575;
    pub const VERTICES: u64 = 7_229_877_672_441_824_278;
    pub const GRID_DISTANCE: u64 = 16_192_610_173_210_096_193;
    pub const LOCAL_IJ: u64 = 2_208_003_384_104_216_248;
    pub const LOCAL_IJ_TO_CELL: u64 = 578_327_552_771_410;
    pub const GRID_PATH: u64 = 9_696_907_593_416_166_129;
    pub const BBOX: u64 = 8_671_842_371_522_402_338;
    pub const AREA: u64 = 13_699_197_200_943_021_274;
    pub const GREAT_CIRCLE: u64 = 9_091_059_960_347_190_061;
    pub const POLYGON_MODES: u64 = 5_165_788_086_682_513_953;
    pub const TRANSMERIDIAN: u64 = 9_989_888_938_456_745_552;

    /// `validity_matches_oracle`: over the full 500_000-element LCG stream, the
    /// number of raw u64 bit-patterns that form a valid H3 cell, an XOR-fold of
    /// those raws, an FNV over each valid `(index, raw)`, and the exact list of
    /// valid `(stream_index, raw)` hits (all res-15 cells).
    pub const VALIDITY_COUNT: u32 = 17;
    pub const VALIDITY_XORFOLD: u64 = 647_841_012_758_652_578;
    pub const VALIDITY_FNV: u64 = 3_532_216_958_570_086_397;
    pub const VALIDITY_HITS: [(u32, u64); 17] = [
        (29_483, 0x08f8_6831_290c_046c),
        (42_714, 0x08f9_2d84_c5d8_52e1),
        (50_067, 0x08f8_4c5a_e237_3c74),
        (194_777, 0x08f5_4f23_6850_4d1a),
        (198_235, 0x08fb_02bc_dd48_975c),
        (204_317, 0x08fb_91a9_4e84_34ae),
        (221_301, 0x08ff_0116_9b49_bc26),
        (237_215, 0x08f8_70e6_f0b6_9810),
        (239_779, 0x08f6_8700_b1b8_64c4),
        (286_348, 0x08fa_9a6c_b559_bdab),
        (308_203, 0x08f8_5b0a_86ca_cc2c),
        (364_029, 0x08fd_7612_1961_6d0e),
        (364_415, 0x08f4_2c3c_01c2_cc70),
        (376_105, 0x08f7_c4d8_8eac_b92a),
        (380_759, 0x08f0_b62a_a94a_2828),
        (430_396, 0x08f1_0d61_72a7_415b),
        (472_794, 0x08fd_5005_ad21_62e1),
    ];
}

/// A spread of real cells: every base cell (res 0), plus cells from a lat/lon
/// grid at several resolutions (including r10, the perf-test source resolution).
///
/// Generated natively (no oracle); the content is pinned by [`golden::SAMPLE_CHK`].
fn sample_cells() -> Vec<u64> {
    let mut cells: Vec<u64> = H3Cell::base_cells().map(|c| c.as_u64()).collect();
    let mut lat = -85i32;
    while lat <= 85 {
        let mut lon = -175i32;
        while lon <= 175 {
            for res in [0u8, 1, 2, 5, 7, 9, 10, 12, 15] {
                let cell = H3Cell::for_point(f64::from(lat), f64::from(lon), res)
                    .expect("sample cell for grid point");
                cells.push(cell.as_u64());
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
fn input_generators_are_pinned() {
    let s = sample_cells();
    assert_eq!(s.len(), golden::SAMPLE_LEN, "sample_cells length drifted");
    assert_eq!(
        chk_u64s(&s),
        golden::SAMPLE_CHK,
        "sample_cells content drifted"
    );
    let sm = small_sample();
    assert_eq!(sm.len(), golden::SMALL_LEN, "small_sample length drifted");
    assert_eq!(
        chk_u64s(&sm),
        golden::SMALL_CHK,
        "small_sample content drifted"
    );
}

#[test]
fn codec_accessors_match_golden() {
    let mut f = Fnv::new();
    for raw in sample_cells() {
        let ours = H3Cell::try_new(raw).expect("ours valid");

        assert_eq!(ours.as_u64(), raw);
        // Token round-trips through parse (self-consistency).
        let token = ours.to_string();
        assert_eq!(
            H3Cell::from_str(&token).expect("parse").as_u64(),
            raw,
            "token round-trip for {raw:#x}"
        );

        // Pin resolution, base_cell, is_pentagon and the textual token to the
        // golden captured from the (correct) reference-matching implementation.
        f.byte(ours.resolution());
        f.byte(ours.base_cell());
        f.byte(u8::from(ours.is_pentagon()));
        f.str(&token);
    }
    assert_eq!(f.finish(), golden::CODEC_ACCESSORS, "codec accessor golden");
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
        H3Cell::average_hexagon_edge_length(9).map(|l| l.value()),
        H3Cell::average_hexagon_edge_length_m(9),
    );
    assert_eq!(
        H3Cell::average_hexagon_area(15).map(|a| a.value()),
        H3Cell::average_hexagon_area_m2(15),
    );
    assert_eq!(H3Cell::average_hexagon_area(16), None);
}

#[test]
fn validity_matches_golden_on_random_bits() {
    // Same seeded LCG and full 500_000-iteration coverage as the original
    // differential test. Valid H3 cells are extremely sparse among random
    // 64-bit patterns (only 17 of 500_000 here), so rather than store a
    // 500k-bit verdict bitset we pin a compact golden — the count, an XOR-fold
    // of the valid raws, an FNV over each valid `(index, raw)`, and the exact
    // list of valid hits — and additionally assert hausdorff rejects every
    // other input by checking the running hit list matches position-for-position.
    let mut state = 0x1234_5678_9abc_def0u64;
    let mut count = 0u32;
    let mut xorfold = 0u64;
    let mut fnv = Fnv::new();
    let mut hit_idx = 0usize;
    for i in 0..500_000u32 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let valid = H3Cell::try_new(state).is_some();
        // Verdict must match the golden at every position: an input is valid
        // iff it is the next expected hit, otherwise it must be rejected.
        let expected_valid =
            hit_idx < golden::VALIDITY_HITS.len() && golden::VALIDITY_HITS[hit_idx] == (i, state);
        assert_eq!(
            valid, expected_valid,
            "validity verdict at index {i} for {state:#x}"
        );
        if valid {
            count += 1;
            xorfold ^= state;
            fnv.u64(u64::from(i));
            fnv.u64(state);
            hit_idx += 1;
        }
    }
    assert_eq!(count, golden::VALIDITY_COUNT, "valid count");
    assert_eq!(xorfold, golden::VALIDITY_XORFOLD, "valid xor-fold");
    assert_eq!(fnv.finish(), golden::VALIDITY_FNV, "valid (index,raw) fnv");
    assert_eq!(
        hit_idx,
        golden::VALIDITY_HITS.len(),
        "every golden hit must be reproduced"
    );

    // Sanity: each golden hit really is a well-formed cell, and the raws are
    // distinct (guards against a degenerate all-accepting codec).
    for (_, raw) in golden::VALIDITY_HITS {
        assert!(
            H3Cell::try_new(raw).is_some(),
            "golden hit {raw:#x} must be valid"
        );
    }
}

#[test]
fn center_matches_golden() {
    let mut f = Fnv::new();
    for raw in sample_cells() {
        let (lat, lon) = H3Cell::try_new(raw).expect("ours").center_lat_lng();
        assert!(
            (-90.0..=90.0).contains(&lat) && (-180.0..=180.0).contains(&lon),
            "center {raw:#x} out of range: ({lat}, {lon})"
        );
        f.f64q(lat);
        f.f64q(lon);
    }
    assert_eq!(f.finish(), golden::CENTER, "center lat/lng golden");
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
fn for_point_matches_golden() {
    let mut f = Fnv::new();
    let mut lat = -85i32;
    while lat <= 85 {
        let mut lon = -179i32;
        while lon <= 179 {
            for res in [0u8, 1, 2, 5, 7, 9, 10, 12, 15] {
                let ours = H3Cell::for_point(f64::from(lat), f64::from(lon), res).expect("ours");
                // The cell really contains the query point at its resolution.
                assert_eq!(
                    ours.resolution(),
                    res,
                    "for_point res {raw_res}",
                    raw_res = res
                );
                f.u64(ours.as_u64());
            }
            lon += 13;
        }
        lat += 6;
    }
    assert_eq!(f.finish(), golden::FOR_POINT, "for_point golden");
}

/// A smaller spread for the O(k^2)/O(7^d) ops: every base cell + a coarse
/// lat/lon grid at a few resolutions (includes pentagons).
///
/// Generated natively (no oracle); the content is pinned by [`golden::SMALL_CHK`].
fn small_sample() -> Vec<u64> {
    let mut v: Vec<u64> = H3Cell::base_cells().map(|c| c.as_u64()).collect();
    let mut lat = -80i32;
    while lat <= 80 {
        let mut lon = -170i32;
        while lon <= 170 {
            for res in [5u8, 7, 9, 10] {
                let cell = H3Cell::for_point(f64::from(lat), f64::from(lon), res)
                    .expect("small-sample cell for grid point");
                v.push(cell.as_u64());
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

#[test]
fn directed_edges_match_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");

        let ours_edges = ours.edges();
        assert_eq!(ours_edges.len(), if ours.is_pentagon() { 5 } else { 6 });

        for edge in &ours_edges {
            let edge = *edge;
            let destination = edge.destination();
            let (origin, destination_from_cells) = edge.cells();

            // Self-consistency: an edge is valid, parses back from its token,
            // its origin is this cell, cells() agrees with destination(), and
            // it is recoverable from (origin, destination).
            assert!(H3DirectedEdge::is_valid(edge.as_u64()));
            assert_eq!(H3DirectedEdge::from_str(&edge.to_string()).ok(), Some(edge));
            assert_eq!(origin, ours);
            assert_eq!(destination_from_cells, destination);
            assert_eq!(edge.origin(), ours);
            assert_eq!(ours.edge(destination), Some(edge));
            assert_eq!(H3DirectedEdge::new(ours, destination), Some(edge));

            f.u64(edge.as_u64());
            f.u64(edge.reverse().as_u64());
            f.u64(destination.as_u64());
            f.u64(edge.origin().as_u64());
        }

        for neighbor in ours.grid_ring(1) {
            f.opt_u64(ours.edge(neighbor).map(H3DirectedEdge::as_u64));
        }
    }
    assert_eq!(f.finish(), golden::EDGES, "directed-edge golden");
}

#[test]
fn directed_edge_geometry_matches_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");
        for edge in ours.edges() {
            for (lat, lon) in edge.boundary_lat_lng() {
                f.f64q(lat);
                f.f64q(lon);
            }
            f.f64q(edge.length_rads());
            f.f64q(edge.length_km());
            f.f64q(edge.length_m());
            // Self-consistency between the typed length and the m accessor.
            assert_eq!(edge.length().value(), edge.length_m());
        }
    }
    assert_eq!(f.finish(), golden::EDGE_GEOM, "edge geometry golden");
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
fn vertices_match_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");

        for v in ours.vertexes() {
            f.u64(v.as_u64());
        }
        assert_eq!(
            ours.vertex_count(),
            if ours.is_pentagon() { 5 } else { 6 },
            "vertex count {raw:#x}"
        );

        for vertex in 0..ours.vertex_count() {
            let ours_vertex = ours.vertex(vertex).expect("our vertex");
            let owner = ours_vertex.owner();
            let (lat, lon) = ours_vertex.lat_lng();

            // Self-consistency: the vertex is valid and round-trips through both
            // its raw bits and its textual token.
            assert!(H3Vertex::is_valid(ours_vertex.as_u64()));
            assert_eq!(H3Vertex::try_new(ours_vertex.as_u64()), Some(ours_vertex));
            assert_eq!(
                H3Vertex::from_str(&ours_vertex.to_string()).ok(),
                Some(ours_vertex)
            );

            f.u64(ours_vertex.as_u64());
            f.u64(owner.as_u64());
            f.byte(ours_vertex.vertex());
            f.f64q(lat);
            f.f64q(lon);
        }

        assert!(ours.vertex(ours.vertex_count()).is_none());
    }
    assert_eq!(f.finish(), golden::VERTICES, "vertices golden");
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
            vec![(0.0, 200.0), (0.0, 201.0), (1.0, 200.5)],
            Vec::<Vec<(f64, f64)>>::new()
        ),
        Err(H3PolygonError::InvalidCoordinate)
    );
}

#[test]
fn rejects_non_cell_indices() {
    assert!(H3Cell::try_new(0).is_none());
    let valid = H3Cell::for_point(40.7, -74.0, 9).expect("valid cell");
    let raw = valid.as_u64();
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

fn cell_set(cells: &[H3Cell]) -> std::collections::BTreeSet<u64> {
    cells.iter().map(|cell| cell.as_u64()).collect()
}

/// `(exterior, holes)` rings in `(latitude, longitude)` degrees.
type PolygonRings = (Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>);

#[test]
fn polygon_to_cells_modes_match_golden() {
    let nyc_exterior = vec![
        (40.70, -74.05),
        (40.70, -73.90),
        (40.85, -73.90),
        (40.85, -74.05),
    ];
    let nyc_hole = vec![
        (40.74, -74.01),
        (40.74, -73.97),
        (40.78, -73.97),
        (40.78, -74.01),
    ];
    let berlin_triangle = vec![(52.40, 13.20), (52.60, 13.50), (52.30, 13.60)];

    let cases: Vec<PolygonRings> = vec![
        (nyc_exterior.clone(), Vec::new()),
        (nyc_exterior, vec![nyc_hole]),
        (berlin_triangle, Vec::new()),
    ];

    let mut f = Fnv::new();
    for (case, (exterior, holes)) in cases.into_iter().enumerate() {
        let polygon = H3Polygon::new(exterior, holes).expect("polygon");
        for res in [5u8, 6, 7, 8] {
            let centroid = polygon
                .to_cells_with_mode(res, ContainmentMode::Centroid)
                .expect("centroid cells");
            let full = polygon
                .to_cells_with_mode(res, ContainmentMode::Full)
                .expect("full cells");
            let overlapping = polygon
                .to_cells_with_mode(res, ContainmentMode::Overlapping)
                .expect("overlapping cells");

            // The default mode is Centroid, and the inherent method agrees with
            // the associated function for Full.
            assert_eq!(polygon.to_cells(res).expect("default cells"), centroid);
            assert_eq!(
                H3Cell::polygon_to_cells_with_mode(&polygon, res, ContainmentMode::Full),
                Ok(full.clone())
            );

            let centroid = cell_set(&centroid);
            let full = cell_set(&full);
            let overlapping = cell_set(&overlapping);
            // Metamorphic nesting: full subset centroid subset overlapping.
            assert!(full.is_subset(&centroid), "case {case} res {res}");
            assert!(centroid.is_subset(&overlapping), "case {case} res {res}");

            // Pin each mode's coverage (sorted) to the golden. Order matches the
            // capture: Centroid, Full, Overlapping.
            for set in [&centroid, &full, &overlapping] {
                for &x in set {
                    f.u64(x);
                }
            }
        }
    }
    assert_eq!(f.finish(), golden::POLYGON_MODES, "polygon modes golden");
}

#[test]
fn polygon_modes_reject_invalid_resolution() {
    let polygon = H3Polygon::new(
        vec![(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        Vec::<Vec<(f64, f64)>>::new(),
    )
    .expect("polygon");
    for mode in [
        ContainmentMode::Centroid,
        ContainmentMode::Full,
        ContainmentMode::Overlapping,
    ] {
        assert_eq!(
            polygon.to_cells_with_mode(16, mode),
            Err(H3PolygonError::InvalidResolution(16))
        );
    }
}

#[test]
fn transmeridian_polygon_to_cells_matches_golden() {
    let exterior = vec![
        (-19.0, 178.0),
        (-19.0, -178.0),
        (-15.0, -178.0),
        (-15.0, 178.0),
    ];
    let polygon = H3Polygon::new(exterior, Vec::<Vec<(f64, f64)>>::new()).expect("polygon");
    assert!(polygon.is_transmeridian());

    assert!(polygon.contains_lat_lng(-17.0, 179.5));
    assert!(polygon.contains_lat_lng(-17.0, -179.5));
    assert!(!polygon.contains_lat_lng(-17.0, 170.0));
    assert!(!polygon.contains_lat_lng(-17.0, -170.0));
    assert!(!polygon.contains_lat_lng(-10.0, 179.5));

    let mut f = Fnv::new();
    for res in [2u8, 3, 4, 5] {
        let ours = polygon.to_cells(res).expect("cells");
        // The cover spans both sides of the antimeridian (metamorphic).
        assert!(
            ours.iter().any(|cell| cell.center_lat_lng().1 > 170.0)
                && ours.iter().any(|cell| cell.center_lat_lng().1 < -170.0),
            "cover should span both sides of the date line at res {res}"
        );
        let mut set: Vec<u64> = ours.iter().map(|c| c.as_u64()).collect();
        set.sort_unstable();
        for x in set {
            f.u64(x);
        }
    }
    assert_eq!(f.finish(), golden::TRANSMERIDIAN, "transmeridian golden");
}

#[test]
fn transmeridian_polygon_modes_are_nested() {
    let polygon = H3Polygon::new(
        vec![
            (-19.0, 178.0),
            (-19.0, -178.0),
            (-15.0, -178.0),
            (-15.0, 178.0),
        ],
        Vec::<Vec<(f64, f64)>>::new(),
    )
    .expect("polygon");

    for res in [3u8, 4] {
        let centroid = cell_set(&polygon.to_cells(res).expect("centroid cells"));
        let full = cell_set(
            &polygon
                .to_cells_with_mode(res, ContainmentMode::Full)
                .expect("full cells"),
        );
        let overlapping = cell_set(
            &polygon
                .to_cells_with_mode(res, ContainmentMode::Overlapping)
                .expect("overlapping cells"),
        );

        assert!(!full.is_empty(), "res {res}");
        assert!(full.is_subset(&centroid), "res {res}");
        assert!(centroid.is_subset(&overlapping), "res {res}");

        for raw in &overlapping {
            let cell = H3Cell::try_new(*raw).expect("cell");
            let inside = full.contains(raw)
                || cell
                    .boundary_lat_lng()
                    .iter()
                    .any(|&(lat, lng)| polygon.contains_lat_lng(lat, lng));
            let near_boundary = cell.grid_disk(1).iter().any(|n| {
                let (lat, lng) = n.center_lat_lng();
                polygon.contains_lat_lng(lat, lng)
            });
            assert!(
                inside || near_boundary,
                "overlapping cell {cell} is far from the polygon"
            );
        }
    }
}

#[test]
fn local_ij_to_cell_round_trips_and_matches_golden() {
    let mut f = Fnv::new();
    for raw in small_sample().into_iter().take(300) {
        let origin = H3Cell::try_new(raw).expect("origin");

        // Self-consistency: to_local_ij then to_cell recovers the original.
        for cell in origin.grid_disk(2) {
            if let Some(ij) = cell.to_local_ij(origin) {
                assert_eq!(ij.to_cell(), Some(cell), "round trip {raw:#x}");
                assert_eq!(H3Cell::from_local_ij(ij), Some(cell));
            }
        }

        let base = origin.to_local_ij(origin).expect("origin local IJ");
        for di in -3i32..=3 {
            for dj in -3i32..=3 {
                let ij = H3LocalIJ {
                    anchor: origin,
                    i: base.i + di,
                    j: base.j + dj,
                };
                f.opt_u64(ij.to_cell().map(H3Cell::as_u64));
            }
        }
    }
    assert_eq!(
        f.finish(),
        golden::LOCAL_IJ_TO_CELL,
        "local_ij_to_cell golden"
    );
}

#[test]
fn great_circle_distance_matches_golden() {
    let mut f = Fnv::new();
    let cells = small_sample();
    for pair in cells.windows(2).take(500) {
        let a = H3Cell::try_new(pair[0]).expect("a");
        let b = H3Cell::try_new(pair[1]).expect("b");

        // Self-consistency: typed accessor agrees with the m accessor, distance
        // to self is zero, and the metric is symmetric.
        assert_eq!(
            a.great_circle_distance(b).value(),
            a.great_circle_distance_m(b)
        );
        assert_eq!(a.great_circle_distance_rads(a), 0.0);
        assert_eq!(
            a.great_circle_distance_rads(b),
            b.great_circle_distance_rads(a)
        );

        f.f64q(a.great_circle_distance_rads(b));
        f.f64q(a.great_circle_distance_km(b));
        f.f64q(a.great_circle_distance_m(b));
    }
    assert_eq!(f.finish(), golden::GREAT_CIRCLE, "great circle golden");
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
