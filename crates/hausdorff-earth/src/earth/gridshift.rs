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
mod tests;
