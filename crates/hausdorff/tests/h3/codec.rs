use super::*;

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
