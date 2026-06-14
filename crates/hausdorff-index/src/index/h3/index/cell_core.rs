//! [`H3Cell`] construction, accessors, and point projection: `try_new` /
//! `from_raw_parts` and the base-cell / pentagon enumerators, the bit-field
//! accessors (`resolution`, `base_cell`, `digit`, `is_*`), and the
//! latitude/longitude ↔ cell projections (`center_lat_lng`, `for_point(s)`).

use super::*;

impl H3Cell {
    /// True iff `bits` is a valid H3 cell index.
    pub fn is_valid(bits: u64) -> bool {
        is_valid_cell(bits)
    }

    /// Every resolution-0 base cell in canonical base-cell order.
    pub fn base_cells() -> impl Iterator<Item = Self> {
        (0..base_cells::COUNT).map(|base_cell| Self(set_base_cell(DEFAULT_CELL_INDEX, base_cell)))
    }

    /// The 12 pentagon cells at `resolution`.
    pub fn pentagons(resolution: u8) -> Vec<Self> {
        if resolution > MAX_RESOLUTION {
            return Vec::new();
        }
        Self::base_cells()
            .filter(|cell| cell.is_pentagon())
            .filter_map(|cell| cell.center_child(resolution))
            .collect()
    }

    /// Constructs a valid cell from raw parts.
    ///
    /// `digits` must contain exactly `resolution` digits, where each digit is in
    /// `0..=6`. Pentagon deleted-K-axis combinations are rejected.
    pub fn from_raw_parts(resolution: u8, base_cell: u8, digits: &[u8]) -> Option<Self> {
        if resolution > MAX_RESOLUTION
            || base_cell >= base_cells::COUNT
            || digits.len() != usize::from(resolution)
        {
            return None;
        }

        let mut bits = set_base_cell(set_resolution(DEFAULT_CELL_INDEX, resolution), base_cell);
        for (idx, digit) in digits.iter().copied().enumerate() {
            if digit > 6 {
                return None;
            }
            bits = set_direction(bits, digit, (idx + 1) as u8);
        }
        Self::try_new(bits)
    }

    /// Wraps a raw u64 H3 index, validating it is a well-formed cell.
    pub fn try_new(bits: u64) -> Option<Self> {
        is_valid_cell(bits).then_some(Self(bits))
    }

    /// The raw u64 H3 index.
    pub const fn as_u64(self) -> u64 {
        self.0
    }

    /// The cell's resolution (0..=15).
    pub const fn resolution(self) -> u8 {
        get_resolution(self.0)
    }

    /// The cell's base cell (0..=121).
    pub const fn base_cell(self) -> u8 {
        get_base_cell(self.0)
    }

    /// The 1-indexed H3 direction digit at resolution `resolution`.
    ///
    /// Returns `None` when `resolution` is outside `1..=15`. For digit slots
    /// beyond this cell's own resolution the returned value is the H3 unused
    /// digit sentinel `7`.
    pub const fn digit(self, resolution: u8) -> Option<u8> {
        if resolution == 0 || resolution > MAX_RESOLUTION {
            return None;
        }
        Some(get_digit(self.0, resolution))
    }

    /// True iff this cell's resolution is Class III.
    pub const fn is_res_class_iii(self) -> bool {
        self.resolution() & 1 == 1
    }

    /// True iff this cell is a pentagon: its base cell is a pentagon and every
    /// resolution digit is `0` (center) — pentagon-ness only descends through
    /// center children.
    pub fn is_pentagon(self) -> bool {
        base_cells::is_pentagon(self.base_cell()) && !has_nonzero_digit(self.0)
    }

    /// The cell center as geographic `(latitude, longitude)` in **degrees**.
    pub fn center_lat_lng(self) -> (f64, f64) {
        let (lat, lng) = self.center_lat_lng_rad();
        (lat.to_degrees(), lng.to_degrees())
    }

    pub(super) fn center_lat_lng_rad(self) -> (f64, f64) {
        let res = self.resolution() as usize;
        faceijk::FaceIJK::from_cell(self.0, self.base_cell(), res).to_center_lat_lng_rad(res)
    }

    /// The cell at `resolution` containing the given geographic point (degrees).
    /// Returns `None` for non-finite or out-of-range coordinates
    /// (`lat ∈ [-90, 90]`, `lng ∈ [-180, 180]`) or an out-of-range resolution.
    pub fn for_point(lat_deg: f64, lng_deg: f64, resolution: u8) -> Option<Self> {
        if !lat_deg.is_finite()
            || !lng_deg.is_finite()
            || !(-90.0..=90.0).contains(&lat_deg)
            || !(-180.0..=180.0).contains(&lng_deg)
            || resolution > MAX_RESOLUTION
        {
            return None;
        }
        let v = coord::Vec3d::from_lat_lng_rad(lat_deg.to_radians(), lng_deg.to_radians());
        let bits =
            faceijk::FaceIJK::from_vec3d(v, resolution as usize).to_cell_bits(resolution as usize);
        Self::try_new(bits)
    }

    /// The cell at `resolution` containing the given [`GeoPoint`](hausdorff_base::earth::GeoPoint).
    /// Convenience over [`for_point`](Self::for_point).
    pub fn for_geo_point(point: hausdorff_base::earth::GeoPoint, resolution: u8) -> Option<Self> {
        Self::for_point(point.lat, point.lon, resolution)
    }

    /// The cells for a column of `(latitude, longitude)` points in degrees — the
    /// batch form of [`for_point`](Self::for_point) over the struct-of-arrays
    /// layout a column-oriented caller already has. `lats` and `lngs` are zipped
    /// to the shorter length; the per-point rule is identical to `for_point`
    /// (`None` for a non-finite or out-of-range coordinate), and the result is
    /// all-`None` when `resolution` is out of range. The resolution is validated
    /// once and the output is written contiguously; the per-point projection is
    /// the same code as `for_point`, so results are bit-identical.
    #[must_use]
    pub fn for_points(lats: &[f64], lngs: &[f64], resolution: u8) -> Vec<Option<Self>> {
        let n = lats.len().min(lngs.len());
        if resolution > MAX_RESOLUTION {
            return vec![None; n];
        }
        let res = resolution as usize;
        lats.iter()
            .zip(lngs)
            .map(|(&lat, &lng)| {
                if !lat.is_finite()
                    || !lng.is_finite()
                    || !(-90.0..=90.0).contains(&lat)
                    || !(-180.0..=180.0).contains(&lng)
                {
                    return None;
                }
                let v = coord::Vec3d::from_lat_lng_rad(lat.to_radians(), lng.to_radians());
                let bits = faceijk::FaceIJK::from_vec3d(v, res).to_cell_bits(res);
                Self::try_new(bits)
            })
            .collect()
    }

    /// As [`for_points`](Self::for_points) but using the faster, approximate
    /// gnomonic projection: it trades the exact `acos` + `tan` for the identity
    /// `tan(acos(x)) = sqrt(1 - x²)/x` (one `sqrt` + a divide). Results are
    /// identical to [`for_points`](Self::for_points) for points in a cell's
    /// interior and can differ only for a point essentially on a cell boundary,
    /// where it lands in an adjacent — still valid — cell. Prefer it for
    /// bucketing / partitioning where that boundary-epsilon difference is
    /// acceptable; use [`for_points`](Self::for_points) when you need the exact
    /// standard cell.
    #[must_use]
    pub fn for_points_approx(lats: &[f64], lngs: &[f64], resolution: u8) -> Vec<Option<Self>> {
        let n = lats.len().min(lngs.len());
        if resolution > MAX_RESOLUTION {
            return vec![None; n];
        }
        let res = resolution as usize;
        lats.iter()
            .zip(lngs)
            .map(|(&lat, &lng)| {
                if !lat.is_finite()
                    || !lng.is_finite()
                    || !(-90.0..=90.0).contains(&lat)
                    || !(-180.0..=180.0).contains(&lng)
                {
                    return None;
                }
                let v = coord::Vec3d::from_lat_lng_rad(lat.to_radians(), lng.to_radians());
                let bits = faceijk::FaceIJK::from_vec3d_fast(v, res).to_cell_bits(res);
                Self::try_new(bits)
            })
            .collect()
    }
}
