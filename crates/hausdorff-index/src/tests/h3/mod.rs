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

use crate::index::h3::{
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
    // Moved once, deliberately: the ring-area reduction changed from a hand-rolled
    // plain Kahan adder to euclid's `compensated::sum`, which is
    // Kahan–Babuška–Neumaier. The two differ by at most 7.5e-14 relative over
    // these cells, and Neumaier is the correct one of the pair — on the
    // discriminating sequence `[1, 1e100, 1, -1e100]` plain Kahan returns 0 where
    // the exact answer, which Neumaier gives, is 2.
    //
    // Worth knowing if this moves again: `f64q` quantises at an *absolute* 1e-9,
    // and `area_m2` runs to ~1e12, so `v * 1e9` saturates `i64` and this hash is
    // far more sensitive than a 1e-9 tolerance suggests. A last-ULP change in any
    // area moves it.
    pub const AREA: u64 = 17_555_332_431_189_521_501;
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

mod codec;

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

mod hierarchy;

mod grid;

mod directed_edges;

mod vertices;

mod grid_paths;

mod polygon;

mod line_cover;

mod polygon_modes;

mod batch;

mod pentagons;
