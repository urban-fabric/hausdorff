//! `gridshift` — high-accuracy datum transforms via **NTv2 grid-shift files**
//! (the `.gsb` industry standard), embedded so there is **no runtime file
//! access**.
//!
//! A 7-parameter [`Helmert`](crate::earth::crs::Helmert) shift is a single rigid
//! transform for a whole datum — good to ~1 m. National mapping agencies publish
//! **NTv2 grids** that model the residual distortion node-by-node, reaching
//! ~cm accuracy (e.g. OSTN15 for OSGB, BETA2007 for the German DHDN). An NTv2
//! grid is just a regular lattice of `(lat_shift, lon_shift)` values plus a small
//! header (bounds, spacing, dimensions); a point is shifted by **bilinear
//! interpolation** of its four surrounding nodes.
//!
//! The standard `.gsb` file is committed and pulled in with
//! [`include_bytes!`] — so it is baked into the compiled library **at compile
//! time** and parsed in memory; at run time it is a `&'static [u8]`, never a file
//! read (no PROJ data dir, no paths, no mtimes — deterministic). The committed
//! artifact stays the *verifiable, official* grid (diff it against the agency
//! distribution; drop in a newer one), not a lossy custom blob.

use std::sync::OnceLock;

use crate::earth::GeoPoint;

/// One **subgrid** of an NTv2 file: a south-west-origin lattice of
/// `(lat_shift, lon_shift)` node values (arc-seconds) over a rectangular extent.
#[derive(Clone, Debug, PartialEq)]
struct SubGrid {
    /// South edge latitude (degrees).
    lat0: f64,
    /// West edge longitude (east-positive degrees).
    lon0: f64,
    /// Node spacing (degrees).
    dlat: f64,
    dlon: f64,
    /// Latitude nodes (south→north) / longitude nodes (west→east).
    n_rows: usize,
    n_cols: usize,
    /// Canonical lattice `[lat_shift, lon_shift]` (arc-seconds), row-major from
    /// the SW corner, latitude increasing then longitude increasing.
    shifts: Vec<[f32; 2]>,
}

/// An NTv2 grid-shift dataset: the **subgrids** of a `.gsb` file (a single grid,
/// or a coarse parent with finer nested children). A point is shifted by the
/// **finest subgrid that contains it** — bilinear interpolation of its four
/// surrounding nodes.
#[derive(Clone, Debug, PartialEq)]
pub struct GridShift {
    subgrids: Vec<SubGrid>,
}

#[inline]
fn i32_le(b: &[u8], off: usize) -> i32 {
    i32::from_le_bytes(b[off..off + 4].try_into().unwrap())
}
#[inline]
fn f64_le(b: &[u8], off: usize) -> f64 {
    f64::from_le_bytes(b[off..off + 8].try_into().unwrap())
}
#[inline]
fn f32_le(b: &[u8], off: usize) -> f32 {
    f32::from_le_bytes(b[off..off + 4].try_into().unwrap())
}

impl SubGrid {
    /// Parse one subgrid — its header (`num_srec` records) and node block — at byte
    /// `off`; returns the subgrid and the offset of the **next** subgrid. Each
    /// header record is an 8-byte ASCII tag + 8-byte value; the nodes are
    /// `GS_COUNT` records of four `f32` (`lat_shift, lon_shift, lat_acc, lon_acc`,
    /// seconds). NTv2 longitude is **positive-west** and nodes run from the SE
    /// corner west then north — both normalised here into the SW-origin lattice.
    fn parse(b: &[u8], off: usize, num_srec: usize) -> (SubGrid, usize) {
        let (mut s_lat, mut n_lat) = (0.0, 0.0);
        let (mut e_long, mut w_long) = (0.0, 0.0);
        let (mut lat_inc, mut long_inc) = (0.0, 0.0);
        let mut gs_count = 0usize;
        for r in 0..num_srec {
            let o = off + r * 16;
            let v = o + 8;
            match &b[o..o + 8] {
                b"S_LAT   " => s_lat = f64_le(b, v),
                b"N_LAT   " => n_lat = f64_le(b, v),
                b"E_LONG  " => e_long = f64_le(b, v),
                b"W_LONG  " => w_long = f64_le(b, v),
                b"LAT_INC " => lat_inc = f64_le(b, v),
                b"LONG_INC" => long_inc = f64_le(b, v),
                b"GS_COUNT" => gs_count = i32_le(b, v) as usize,
                _ => {}
            }
        }

        let node_base = off + num_srec * 16;
        let n_rows = ((n_lat - s_lat) / lat_inc).round() as usize + 1;
        let n_cols = ((w_long - e_long) / long_inc).round() as usize + 1;
        let mut shifts = Vec::with_capacity(n_rows * n_cols);
        for r in 0..n_rows {
            for c in 0..n_cols {
                // canonical col `c` is west→east; NTv2 native col is east→west.
                let j = n_cols - 1 - c;
                let o = node_base + (r * n_cols + j) * 16;
                shifts.push([f32_le(b, o), f32_le(b, o + 4)]);
            }
        }

        let sub = SubGrid {
            lat0: s_lat / 3600.0,
            lon0: -w_long / 3600.0, // east-positive west edge
            dlat: lat_inc / 3600.0,
            dlon: long_inc / 3600.0,
            n_rows,
            n_cols,
            shifts,
        };
        (sub, node_base + gs_count * 16)
    }

    /// True iff `(lat, lon)` (east-positive degrees) is inside this subgrid.
    fn contains(&self, lat: f64, lon: f64) -> bool {
        lat >= self.lat0
            && lat <= self.lat0 + (self.n_rows - 1) as f64 * self.dlat
            && lon >= self.lon0
            && lon <= self.lon0 + (self.n_cols - 1) as f64 * self.dlon
    }

    /// Bilinearly-interpolated `(lat_shift, lon_shift)` in **arc-seconds**.
    fn shift_seconds(&self, lat: f64, lon: f64) -> (f64, f64) {
        let fr = (lat - self.lat0) / self.dlat;
        let fc = (lon - self.lon0) / self.dlon;
        let r0 = (fr.floor() as usize).min(self.n_rows - 2);
        let c0 = (fc.floor() as usize).min(self.n_cols - 2);
        let (tr, tc) = (fr - r0 as f64, fc - c0 as f64);
        let at = |r: usize, c: usize| self.shifts[r * self.n_cols + c];
        let (s00, s01, s10, s11) = (
            at(r0, c0),
            at(r0, c0 + 1),
            at(r0 + 1, c0),
            at(r0 + 1, c0 + 1),
        );
        let bilerp = |a: f32, b: f32, c: f32, d: f32| {
            let south = f64::from(a) * (1.0 - tc) + f64::from(b) * tc;
            let north = f64::from(c) * (1.0 - tc) + f64::from(d) * tc;
            south * (1.0 - tr) + north * tr
        };
        (
            bilerp(s00[0], s01[0], s10[0], s11[0]),
            bilerp(s00[1], s01[1], s10[1], s11[1]),
        )
    }
}

impl GridShift {
    /// Parse a standard NTv2 (`.gsb`) buffer: an overview header of `NUM_OREC`
    /// records (with `NUM_SREC`, `NUM_FILE`), then `NUM_FILE` subgrids each a
    /// header + node block. Multi-subgrid files (parent + nested children) are all
    /// read; selection picks the finest covering subgrid.
    fn parse_ntv2(b: &[u8]) -> GridShift {
        let num_orec = i32_le(b, 8) as usize;
        let num_srec = i32_le(b, 24) as usize;
        let num_file = i32_le(b, 40).max(1) as usize;
        let mut subgrids = Vec::with_capacity(num_file);
        let mut off = num_orec * 16; // first subgrid header
        for _ in 0..num_file {
            let (sub, next) = SubGrid::parse(b, off, num_srec);
            subgrids.push(sub);
            off = next;
        }
        GridShift { subgrids }
    }

    /// The **finest** subgrid containing `(lat, lon)` — a denser nested child wins
    /// over its coarser parent — or `None` if no subgrid covers the point.
    fn select(&self, lat: f64, lon: f64) -> Option<&SubGrid> {
        self.subgrids
            .iter()
            .filter(|s| s.contains(lat, lon))
            .min_by(|a, b| a.dlat.total_cmp(&b.dlat))
    }

    /// True iff `(lat, lon)` (east-positive degrees) lies within any subgrid.
    #[must_use]
    pub fn contains(&self, lat: f64, lon: f64) -> bool {
        self.select(lat, lon).is_some()
    }

    /// Bilinearly-interpolated `(lat_shift, lon_shift)` in **arc-seconds** from the
    /// finest covering subgrid, or `None` if no subgrid covers `(lat, lon)`.
    fn shift_seconds(&self, lat: f64, lon: f64) -> Option<(f64, f64)> {
        Some(self.select(lat, lon)?.shift_seconds(lat, lon))
    }

    /// Shift a point from the grid's **source** datum to its **target** datum
    /// (e.g. DHDN → ETRS89/WGS84). NTv2's longitude shift is positive-west, so the
    /// east-positive longitude decreases by it. Points outside the grid pass
    /// through unchanged (the caller falls back to a Helmert shift).
    #[must_use]
    pub fn forward(&self, p: GeoPoint) -> GeoPoint {
        match self.shift_seconds(p.lat, p.lon) {
            Some((dlat, dlon)) => {
                GeoPoint::from_lat_lon(p.lat + dlat / 3600.0, p.lon - dlon / 3600.0)
            }
            None => p,
        }
    }

    /// Shift a point from the **target** datum back to the **source** datum — the
    /// inverse of [`forward`](Self::forward). The grid is indexed by source
    /// coordinates, so the inverse is found by a few fixed-point iterations.
    #[must_use]
    pub fn inverse(&self, p: GeoPoint) -> GeoPoint {
        let mut guess = p;
        for _ in 0..6 {
            let fwd = self.forward(guess);
            let (dlat, dlon) = (p.lat - fwd.lat, p.lon - fwd.lon);
            guess = GeoPoint::from_lat_lon(guess.lat + dlat, guess.lon + dlon);
            if dlat.abs() < 1e-12 && dlon.abs() < 1e-12 {
                break;
            }
        }
        guess
    }
}

/// The embedded **BETA2007** NTv2 grid (BKG — German Federal Agency for
/// Cartography and Geodesy; DHDN ↔ ETRS89), parsed once on first use from the
/// `.gsb` bytes baked into the library.
#[must_use]
pub fn beta2007() -> &'static GridShift {
    static G: OnceLock<GridShift> = OnceLock::new();
    G.get_or_init(|| GridShift::parse_ntv2(include_bytes!("../../data/BETA2007.gsb")))
}

/// The embedded **ntf_r93** NTv2 grid (IGN France; NTF ↔ RGF93), the official
/// French legacy-datum transform.
#[must_use]
pub fn ntf_r93() -> &'static GridShift {
    static G: OnceLock<GridShift> = OnceLock::new();
    G.get_or_init(|| GridShift::parse_ntv2(include_bytes!("../../data/ntf_r93.gsb")))
}

/// The embedded **nzgd2kgrid0005** NTv2 grid (LINZ New Zealand;
/// NZGD49 ↔ NZGD2000).
#[must_use]
pub fn nzgd2000() -> &'static GridShift {
    static G: OnceLock<GridShift> = OnceLock::new();
    G.get_or_init(|| GridShift::parse_ntv2(include_bytes!("../../data/nzgd2kgrid0005.gsb")))
}

#[cfg(test)]
mod tests {
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
}
