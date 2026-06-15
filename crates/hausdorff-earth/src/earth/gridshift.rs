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

/// An NTv2 grid-shift dataset: the first subgrid of a `.gsb` file, parsed into a
/// south-west-origin lattice of `(lat_shift, lon_shift)` node values (arc-seconds).
#[derive(Clone, Debug, PartialEq)]
pub struct GridShift {
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

impl GridShift {
    /// Parse the first subgrid of a standard NTv2 (`.gsb`) buffer. Each header
    /// record is an 8-byte ASCII tag + 8-byte value; the overview header has
    /// `NUM_OREC` records, the subgrid header `NUM_SREC`, then `GS_COUNT` node
    /// records of four `f32` (`lat_shift, lon_shift, lat_acc, lon_acc`, seconds).
    /// NTv2 longitude is **positive-west** and nodes run from the SE corner west
    /// then north — both normalised here into the SW-origin canonical lattice.
    fn parse_ntv2(b: &[u8]) -> GridShift {
        let num_orec = i32_le(b, 8) as usize;
        let num_srec = i32_le(b, 24) as usize;
        let sub = num_orec * 16; // subgrid header offset

        let (mut s_lat, mut n_lat) = (0.0, 0.0);
        let (mut e_long, mut w_long) = (0.0, 0.0);
        let (mut lat_inc, mut long_inc) = (0.0, 0.0);
        for r in 0..num_srec {
            let o = sub + r * 16;
            let v = o + 8;
            match &b[o..o + 8] {
                b"S_LAT   " => s_lat = f64_le(b, v),
                b"N_LAT   " => n_lat = f64_le(b, v),
                b"E_LONG  " => e_long = f64_le(b, v),
                b"W_LONG  " => w_long = f64_le(b, v),
                b"LAT_INC " => lat_inc = f64_le(b, v),
                b"LONG_INC" => long_inc = f64_le(b, v),
                _ => {}
            }
        }

        let node_base = (num_orec + num_srec) * 16;
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

        GridShift {
            lat0: s_lat / 3600.0,
            lon0: -w_long / 3600.0, // east-positive west edge
            dlat: lat_inc / 3600.0,
            dlon: long_inc / 3600.0,
            n_rows,
            n_cols,
            shifts,
        }
    }

    /// True iff `(lat, lon)` (east-positive degrees) is inside the grid extent.
    #[must_use]
    pub fn contains(&self, lat: f64, lon: f64) -> bool {
        lat >= self.lat0
            && lat <= self.lat0 + (self.n_rows - 1) as f64 * self.dlat
            && lon >= self.lon0
            && lon <= self.lon0 + (self.n_cols - 1) as f64 * self.dlon
    }

    /// Bilinearly-interpolated `(lat_shift, lon_shift)` in **arc-seconds** at
    /// `(lat, lon)`, or `None` if outside the grid.
    fn shift_seconds(&self, lat: f64, lon: f64) -> Option<(f64, f64)> {
        if !self.contains(lat, lon) {
            return None;
        }
        let fr = (lat - self.lat0) / self.dlat;
        let fc = (lon - self.lon0) / self.dlon;
        let r0 = (fr.floor() as usize).min(self.n_rows - 2);
        let c0 = (fc.floor() as usize).min(self.n_cols - 2);
        let (tr, tc) = (fr - r0 as f64, fc - c0 as f64);
        let at = |r: usize, c: usize| self.shifts[r * self.n_cols + c];
        let (s00, s01, s10, s11) = (at(r0, c0), at(r0, c0 + 1), at(r0 + 1, c0), at(r0 + 1, c0 + 1));
        let bilerp = |a: f32, b: f32, c: f32, d: f32| {
            let south = f64::from(a) * (1.0 - tc) + f64::from(b) * tc;
            let north = f64::from(c) * (1.0 - tc) + f64::from(d) * tc;
            south * (1.0 - tr) + north * tr
        };
        Some((
            bilerp(s00[0], s01[0], s10[0], s11[0]),
            bilerp(s00[1], s01[1], s10[1], s11[1]),
        ))
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
        // Germany: 47–55.3°N, 5.5–15.67°E at 0.1°×(1/6)° spacing, 84×62 nodes.
        assert_eq!((g.n_rows, g.n_cols), (84, 62));
        assert!((g.lat0 - 47.0).abs() < 1e-9 && (g.lon0 - 5.5).abs() < 1e-9);
        assert!((g.dlat - 0.1).abs() < 1e-9);
        assert!(g.contains(52.0, 10.0));
        assert!(!g.contains(40.0, 10.0)); // outside (south of Germany)
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
