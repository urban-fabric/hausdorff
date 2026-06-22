use super::*;

#[test]
fn beta2007_grid_parses_with_expected_extent() {
    let g = beta2007();
    assert_eq!(g.subgrids.len(), 1); // single-subgrid file
    let s = &g.subgrids[0];
    // Germany: 47–55.3°N, 5.5–15.67°E at 0.1°×(1/6)° spacing, 84×62 nodes.
    assert_eq!((s.n_rows, s.n_cols), (84, 62));
    assert!((s.lat0 - 47.0).abs() < 1e-9 && (s.lon0 - 5.5).abs() < 1e-9);
    assert!((s.dlat - 0.1).abs() < 1e-9);
    assert!(g.contains(52.0, 10.0));
    assert!(!g.contains(40.0, 10.0)); // outside (south of Germany)
}

/// Build a synthetic **two-subgrid** NTv2 buffer: a coarse parent (50–52°N,
/// 0–2°E, 1° nodes, +10″ lat shift) with a finer child nested inside it
/// (50.5–51°N, 0.5–1°E, 0.25° nodes, +20″ lat shift).
fn two_subgrid_gsb() -> Vec<u8> {
    let mut b = Vec::new();
    let rec = |buf: &mut Vec<u8>, tag: &[u8; 8], val: &[u8]| {
        buf.extend_from_slice(tag);
        let mut v = val.to_vec();
        v.resize(8, 0);
        buf.extend_from_slice(&v);
    };
    // Overview header (NUM_OREC = 11).
    rec(&mut b, b"NUM_OREC", &11i32.to_le_bytes());
    rec(&mut b, b"NUM_SREC", &11i32.to_le_bytes());
    rec(&mut b, b"NUM_FILE", &2i32.to_le_bytes());
    rec(&mut b, b"GS_TYPE ", b"SECONDS ");
    rec(&mut b, b"VERSION ", b"NTv2.0  ");
    rec(&mut b, b"SYSTEM_F", b"FROM    ");
    rec(&mut b, b"SYSTEM_T", b"TO      ");
    rec(&mut b, b"MAJOR_F ", &6_378_137.0f64.to_le_bytes());
    rec(&mut b, b"MINOR_F ", &6_356_752.0f64.to_le_bytes());
    rec(&mut b, b"MAJOR_T ", &6_378_137.0f64.to_le_bytes());
    rec(&mut b, b"MINOR_T ", &6_356_752.0f64.to_le_bytes());
    // One subgrid: header (11 records) + `gs_count` constant nodes.
    let subgrid = |b: &mut Vec<u8>,
                   name: &[u8; 8],
                   parent: &[u8; 8],
                   s_lat: f64,
                   n_lat: f64,
                   e_long: f64,
                   w_long: f64,
                   inc: f64,
                   lat_shift: f32| {
        let gs = (((n_lat - s_lat) / inc).round() as i32 + 1)
            * (((w_long - e_long) / inc).round() as i32 + 1);
        rec(b, b"SUB_NAME", name);
        rec(b, b"PARENT  ", parent);
        rec(b, b"CREATED ", b"00000000");
        rec(b, b"UPDATED ", b"00000000");
        rec(b, b"S_LAT   ", &s_lat.to_le_bytes());
        rec(b, b"N_LAT   ", &n_lat.to_le_bytes());
        rec(b, b"E_LONG  ", &e_long.to_le_bytes());
        rec(b, b"W_LONG  ", &w_long.to_le_bytes());
        rec(b, b"LAT_INC ", &inc.to_le_bytes());
        rec(b, b"LONG_INC", &inc.to_le_bytes());
        rec(b, b"GS_COUNT", &gs.to_le_bytes());
        for _ in 0..gs {
            b.extend_from_slice(&lat_shift.to_le_bytes()); // lat shift
            b.extend_from_slice(&0f32.to_le_bytes()); // lon shift
            b.extend_from_slice(&0f32.to_le_bytes()); // lat acc
            b.extend_from_slice(&0f32.to_le_bytes()); // lon acc
        }
    };
    // Parent (0–2°E in positive-west seconds: E edge 2°E = −7200, W edge 0).
    subgrid(
        &mut b,
        b"PARENT00",
        b"NONE    ",
        180_000.0,
        187_200.0,
        -7200.0,
        0.0,
        3600.0,
        10.0,
    );
    // Child nested inside (0.5–1°E: E edge −3600, W edge −1800).
    subgrid(
        &mut b,
        b"CHILD000",
        b"PARENT00",
        181_800.0,
        183_600.0,
        -3600.0,
        -1800.0,
        900.0,
        20.0,
    );
    b
}

#[test]
fn multi_subgrid_selects_the_finest_covering_grid() {
    let g = GridShift::parse_ntv2(&two_subgrid_gsb());
    assert_eq!(g.subgrids.len(), 2);
    // Inside the child (and the parent): the finer child's +20″ shift wins.
    let child = g.forward(GeoPoint::from_lat_lon(50.7, 0.7));
    assert!(
        (child.lat - (50.7 + 20.0 / 3600.0)).abs() < 1e-9,
        "child shift: {}",
        child.lat
    );
    // Inside the parent only: the parent's +10″ shift.
    let parent = g.forward(GeoPoint::from_lat_lon(50.1, 0.1));
    assert!(
        (parent.lat - (50.1 + 10.0 / 3600.0)).abs() < 1e-9,
        "parent shift: {}",
        parent.lat
    );
    // Outside both subgrids: unchanged.
    assert_eq!(
        g.forward(GeoPoint::from_lat_lon(40.0, 0.5)),
        GeoPoint::from_lat_lon(40.0, 0.5)
    );
}

#[test]
fn matches_proj_beta2007() {
    // PROJ 9.4.0 DHDN (EPSG:4314) → WGS84 (4326), which uses this same grid.
    // (source_lat, source_lon, proj_target_lat, proj_target_lon)
    let refs = [
        (52.0, 10.0, 51.998_628_928, 9.998_786_383),
        (48.5, 9.2, 48.499_021_746, 9.198_951_791),
        (54.0, 13.0, 53.998_410_697, 12.998_281_314),
        (49.0, 7.0, 48.998_951_376, 6.999_247_781),
        (50.5, 11.5, 50.498_811_482, 11.498_582_009),
    ];
    check_against_proj(beta2007(), &refs);
}

#[test]
fn matches_proj_ntf_r93() {
    // PROJ NTF (EPSG:4275) → WGS84 (4326), which uses ntf_r93.gsb.
    check_against_proj(
        ntf_r93(),
        &[
            (48.85, 2.35, 48.849_933_563, 2.349_295_594),
            (43.6, 1.44, 43.599_980_695, 1.439_310_962),
            (45.75, 4.85, 45.749_971_904, 4.849_417_525),
            (47.2, -1.55, 47.199_929_067, -1.550_870_300),
        ],
    );
}

#[test]
fn matches_proj_nzgd2000() {
    // PROJ NZGD49 (EPSG:4272) → WGS84 (4326), which uses nzgd2kgrid0005.gsb.
    check_against_proj(
        nzgd2000(),
        &[
            (-36.85, 174.76, -36.848_196_691, 174.760_191_647),
            (-41.29, 174.78, -41.288_275_516, 174.780_190_614),
            (-43.53, 172.64, -43.528_327_299, 172.640_130_644),
            (-45.03, 168.66, -45.028_369_020, 168.660_071_307),
        ],
    );
}

/// Assert the grid's forward shift matches PROJ to sub-mm at every reference.
fn check_against_proj(g: &GridShift, refs: &[(f64, f64, f64, f64)]) {
    let mut max_deg = 0.0_f64;
    for &(lat, lon, rlat, rlon) in refs {
        let w = g.forward(GeoPoint::from_lat_lon(lat, lon));
        max_deg = max_deg.max((w.lat - rlat).abs()).max((w.lon - rlon).abs());
    }
    assert!(
        max_deg < 1e-7,
        "grid shift vs PROJ off by {max_deg}° (~{:.3} m)",
        max_deg * 111_320.0
    );
}

#[test]
fn forward_inverse_round_trips() {
    let g = beta2007();
    for &(lat, lon) in &[(52.0, 10.0), (48.5, 9.2), (54.0, 13.0), (49.0, 7.0)] {
        let p = GeoPoint::from_lat_lon(lat, lon);
        let back = g.inverse(g.forward(p));
        assert!((back.lat - lat).abs() < 1e-9 && (back.lon - lon).abs() < 1e-9);
    }
}
