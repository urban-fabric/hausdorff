//! H3 index bit codec + validity.
//!
//! An H3 cell index packs, from the high bits down: a reserved high bit (0), a
//! 4-bit mode (1 = cell), 3 reserved bits (0 for cells), a 4-bit resolution, a
//! 7-bit base cell, then 15 three-bit digits (one per resolution level). Digits
//! beyond the cell's resolution are set to the sentinel `7`.
//!
//! Layout (bit offsets):
//! ```text
//!  63        high bit (must be 0)
//!  59..=62   mode (cell = 1)
//!  56..=58   reserved (0 for cells)
//!  52..=55   resolution (0..=15)
//!  45..=51   base cell (0..=121)
//!  0..=44    15 digits x 3 bits  (digit r at bits (15-r)*3 .. +2)
//! ```

use std::collections::BTreeSet;
use std::fmt;
use std::str::FromStr;

use super::{base_cells, children, coord, faceijk, grid};

const MODE_OFFSET: u64 = 59;
const RESERVED_OFFSET: u64 = 56;
const RESOLUTION_OFFSET: u64 = 52;
const BASE_CELL_OFFSET: u64 = 45;

const CELL_MODE: u64 = 1;
const DIRECTED_EDGE_MODE: u64 = 2;
const VERTEX_MODE: u64 = 4;
pub const MAX_RESOLUTION: u8 = 15;
/// Sentinel digit value stored in unused (beyond-resolution) digit slots.
const INVALID_DIGIT: u8 = 7;
/// The K-axis digit — the "deleted subsequence" leading digit for pentagons.
const K_AXIS_DIGIT: u8 = 1;

const EARTH_RADIUS_KM: f64 = 6371.007180918475;
const NUM_HEX_VERTS: u8 = 6;
const NUM_PENT_VERTS: u8 = 5;
const DIRECTIONS_CCW: [u8; 6] = [2, 3, 1, 5, 4, 6];
const REV_NEIGHBOR_DIRECTIONS_HEX: [u8; 6] = [5, 3, 4, 1, 0, 2];
const VERTEX_TO_DIRECTION_HEX: [u8; 6] = [6, 2, 3, 1, 5, 4];
const VERTEX_TO_DIRECTION_PENT: [u8; 5] = [6, 2, 3, 5, 4];
const DIRECTION_TO_VERTEX_HEX: [u8; 6] = [3, 1, 2, 5, 4, 0];
const DIRECTION_TO_VERTEX_PENT: [u8; 5] = [1, 2, 4, 3, 0];
const HEXAGON_EDGE_ORDER: [u8; 6] = [1, 5, 4, 6, 2, 3];
const PENTAGON_EDGE_ORDER: [u8; 5] = [2, 3, 5, 4, 6];

const AVG_HEXAGON_AREA_KM2: [f64; 16] = [
    4.357449416078383e+06,
    6.097884417941332e+05,
    8.680178039899720e+04,
    1.239343465508816e+04,
    1.770347654491307e+03,
    2.529038581819449e+02,
    3.612906216441245e+01,
    5.161293359717191e+00,
    7.373275975944177e-01,
    1.053325134272067e-01,
    1.504750190766435e-02,
    2.149643129451879e-03,
    3.070918756316060e-04,
    4.387026794728296e-05,
    6.267181135324313e-06,
    8.953115907605790e-07,
];

const AVG_HEXAGON_AREA_M2: [f64; 16] = [
    4.357449416078390e+12,
    6.097884417941339e+11,
    8.680178039899731e+10,
    1.239343465508818e+10,
    1.770347654491309e+09,
    2.529038581819452e+08,
    3.612906216441250e+07,
    5.161293359717198e+06,
    7.373275975944188e+05,
    1.053325134272069e+05,
    1.504750190766437e+04,
    2.149643129451882e+03,
    3.070918756316063e+02,
    4.387026794728301e+01,
    6.267181135324322e+00,
    8.953115907605802e-01,
];

const AVG_HEXAGON_EDGE_LENGTH_KM: [f64; 16] = [
    1281.256011,
    483.0568391,
    182.5129565,
    68.97922179,
    26.07175968,
    9.854090990,
    3.724532667,
    1.406475763,
    0.531414010,
    0.200786148,
    0.075863783,
    0.028663897,
    0.010830188,
    0.004092010,
    0.001546100,
    0.000584169,
];

const AVG_HEXAGON_EDGE_LENGTH_M: [f64; 16] = [
    1281256.011,
    483056.8391,
    182512.9565,
    68979.22179,
    26071.75968,
    9854.090990,
    3724.532667,
    1406.475763,
    531.4140101,
    200.7861476,
    75.86378287,
    28.66389748,
    10.83018784,
    4.092010473,
    1.546099657,
    0.584168630,
];

/// An H3 cell index.
///
/// Always a valid cell (mode 1) — construct via [`H3Cell::try_new`] /
/// [`H3Cell::from_str`], both of which reject non-cell or malformed indices the
/// same way the H3 reference does.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct H3Cell(u64);

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
        let res = self.resolution() as usize;
        let fijk = faceijk::FaceIJK::from_cell(self.0, self.base_cell(), res);
        let (lat, lng) = fijk.to_center_lat_lng_rad(res);
        (lat.to_degrees(), lng.to_degrees())
    }

    /// The cell at `resolution` containing the given geographic point (degrees).
    /// Returns `None` for non-finite coordinates or an out-of-range resolution.
    pub fn for_point(lat_deg: f64, lng_deg: f64, resolution: u8) -> Option<Self> {
        if !lat_deg.is_finite() || !lng_deg.is_finite() || resolution > MAX_RESOLUTION {
            return None;
        }
        let v = coord::Vec3d::from_lat_lng_rad(lat_deg.to_radians(), lng_deg.to_radians());
        let bits =
            faceijk::FaceIJK::from_vec3d(v, resolution as usize).to_cell_bits(resolution as usize);
        Self::try_new(bits)
    }

    /// The cell at `resolution` containing the given [`GeoPoint`](crate::earth::GeoPoint).
    /// Convenience over [`for_point`](Self::for_point).
    pub fn for_geo_point(point: crate::earth::GeoPoint, resolution: u8) -> Option<Self> {
        Self::for_point(point.lat, point.lon, resolution)
    }

    /// The parent cell at `resolution` (coarser-or-equal), or `None`.
    pub fn parent(self, resolution: u8) -> Option<Self> {
        (resolution <= self.resolution()).then(|| {
            // Coarsening preserves validity; safe to construct directly.
            Self(set_unused(set_resolution(self.0, resolution), resolution))
        })
    }

    /// All descendant cells at `resolution` (finer-or-equal); empty if coarser
    /// or out of range.
    pub fn children_at_resolution(self, resolution: u8) -> Vec<Self> {
        if resolution < self.resolution() || resolution > MAX_RESOLUTION {
            return Vec::new();
        }
        children::Children::new(self.0, self.resolution(), self.is_pentagon(), resolution)
            .map(Self)
            .collect()
    }

    /// Number of descendant cells at `resolution`.
    pub fn children_count(self, resolution: u8) -> u64 {
        if resolution > MAX_RESOLUTION {
            return 0;
        }
        children::children_count(self.resolution(), resolution, self.is_pentagon())
    }

    /// The center descendant at `resolution`.
    pub fn center_child(self, resolution: u8) -> Option<Self> {
        if resolution < self.resolution() || resolution > MAX_RESOLUTION {
            return None;
        }
        children::Children::new(self.0, self.resolution(), self.is_pentagon(), resolution)
            .next()
            .map(Self)
    }

    /// Position of this cell within its `parent_resolution` parent's ordered
    /// child list. Returns `None` when `parent_resolution` is finer than this
    /// cell.
    pub fn child_position(self, parent_resolution: u8) -> Option<u64> {
        if parent_resolution > self.resolution() || parent_resolution > MAX_RESOLUTION {
            return None;
        }

        let parent_is_pentagon = self.parent(parent_resolution)?.is_pentagon();
        let mut position = 0u64;

        if parent_is_pentagon {
            for res in (parent_resolution + 1)..=self.resolution() {
                let parent = self.parent(res - 1).expect("valid parent");
                let parent_is_pentagon = parent.is_pentagon();
                let mut digit = get_direction(self.0, res);

                if parent_is_pentagon && digit > 0 {
                    digit -= 1;
                }
                if digit == 0 {
                    continue;
                }

                let diff = self.resolution() - res;
                let hex_count = children::hexagon_children_count(diff);
                let count0 = if parent_is_pentagon {
                    children::pentagon_children_count(diff)
                } else {
                    hex_count
                };
                position += u64::from(digit - 1) * hex_count + count0;
            }
        } else {
            for res in (parent_resolution + 1)..=self.resolution() {
                let diff = self.resolution() - res;
                let hex_count = children::hexagon_children_count(diff);
                let digit = get_direction(self.0, res);
                position += u64::from(digit) * hex_count;
            }
        }

        Some(position)
    }

    /// Child cell at `position` within this cell's ordered children at
    /// `resolution`.
    pub fn child_at(self, mut position: u64, resolution: u8) -> Option<Self> {
        if resolution < self.resolution()
            || resolution > MAX_RESOLUTION
            || position >= self.children_count(resolution)
        {
            return None;
        }

        let mut child = set_resolution(self.0, resolution);
        let mut current_res = self.resolution();

        if self.is_pentagon() {
            for res in (self.resolution() + 1)..=resolution {
                current_res = res;
                let diff = resolution - res;
                let pentagon_count = children::pentagon_children_count(diff);
                if position < pentagon_count {
                    child = set_direction(child, 0, res);
                } else {
                    let count = children::hexagon_children_count(diff);
                    position -= pentagon_count;
                    child = set_direction(child, (position / count) as u8 + 2, res);
                    position %= count;
                    break;
                }
            }
        }

        for res in (current_res + 1)..=resolution {
            let diff = resolution - res;
            let count = children::hexagon_children_count(diff);
            child = set_direction(child, (position / count) as u8, res);
            position %= count;
        }

        Self::try_new(child)
    }

    /// All cells within grid distance `k` (the filled disk, including `self`).
    pub fn grid_disk(self, k: u32) -> Vec<Self> {
        grid::grid_disk(self.0, k).into_iter().map(Self).collect()
    }

    /// All cells within grid distance `k`, paired with distance from `self`.
    pub fn grid_disk_distances(self, k: u32) -> Vec<(Self, u32)> {
        grid::grid_disk_distances(self.0, k)
            .into_iter()
            .map(|(cell, distance)| (Self(cell), distance))
            .collect()
    }

    /// Cells at exactly grid distance `k` from `self`.
    pub fn grid_ring(self, k: u32) -> Vec<Self> {
        self.grid_disk_distances(k)
            .into_iter()
            .filter_map(|(cell, distance)| (distance == k).then_some(cell))
            .collect()
    }

    /// Grid distance to another cell, or `None` if undefined (different
    /// resolution, too far, or across a pentagon).
    pub fn grid_distance(self, other: Self) -> Option<i32> {
        grid::grid_distance(self.0, other.0)
    }

    /// This cell's IJ coordinate in `origin`'s local coordinate system.
    pub fn to_local_ij(self, origin: Self) -> Option<H3LocalIJ> {
        let ij = coord::CoordIJ::from(grid::local_ijk(self.0, origin.0)?);
        Some(H3LocalIJ {
            anchor: origin,
            i: ij.i,
            j: ij.j,
        })
    }

    /// True iff `other` is one grid step away from this cell.
    pub fn is_neighbor_with(self, other: Self) -> bool {
        self.grid_distance(other) == Some(1)
    }

    /// Number of cells in a shortest grid path from this cell to `other`,
    /// including both endpoints.
    pub fn grid_path_cells_size(self, other: Self) -> Option<u32> {
        self.grid_distance(other)
            .and_then(|distance| u32::try_from(distance + 1).ok())
    }

    /// A shortest grid path from this cell to `other`, including both
    /// endpoints. Returns `None` when H3 local-IJ distance is undefined.
    pub fn grid_path_cells(self, other: Self) -> Option<Vec<Self>> {
        let mut remaining = self.grid_distance(other)?;
        let mut current = self;
        let mut path = Vec::with_capacity(usize::try_from(remaining + 1).ok()?);
        path.push(current);

        while remaining > 0 {
            let mut neighbors = current.grid_ring(1);
            neighbors.sort_unstable();
            current = neighbors
                .into_iter()
                .find(|neighbor| neighbor.grid_distance(other) == Some(remaining - 1))?;
            path.push(current);
            remaining -= 1;
        }

        Some(path)
    }

    /// Directed edge from this cell to `destination`, if the cells are
    /// neighbors.
    pub fn edge(self, destination: Self) -> Option<H3DirectedEdge> {
        let direction = grid::direction_for_neighbor(self.0, destination.0)?;
        H3DirectedEdge::from_origin_and_direction(self, direction)
    }

    /// All directed edges originating at this cell.
    pub fn edges(self) -> Vec<H3DirectedEdge> {
        let directions: &[u8] = if self.is_pentagon() {
            &PENTAGON_EDGE_ORDER
        } else {
            &HEXAGON_EDGE_ORDER
        };
        directions
            .iter()
            .copied()
            .filter_map(|direction| H3DirectedEdge::from_origin_and_direction(self, direction))
            .collect()
    }

    /// Number of topological vertices on this cell.
    pub fn vertex_count(self) -> u8 {
        if self.is_pentagon() {
            NUM_PENT_VERTS
        } else {
            NUM_HEX_VERTS
        }
    }

    /// Canonical vertex index for `vertex` (`0..5` for hexagons, `0..4` for
    /// pentagons).
    pub fn vertex(self, vertex: u8) -> Option<H3Vertex> {
        let vertex_count = self.vertex_count();
        if vertex >= vertex_count {
            return None;
        }

        let mut owner = self;
        let mut owner_vertex = vertex;

        if self.resolution() == 0 || get_direction(self.0, self.resolution()) != 0 {
            let left = self.direction_for_vertex(vertex)?;
            let (left_bits, left_rotation) = grid::neighbor_rotations(self.0, left, 0)?;
            let left_cell = Self(left_bits);
            if left_cell < owner {
                owner = left_cell;
            }

            if self.resolution() == 0 || get_direction(left_bits, self.resolution()) != 0 {
                let right_vertex = (vertex + vertex_count - 1) % vertex_count;
                let right = self.direction_for_vertex(right_vertex)?;
                let (right_bits, right_rotation) = grid::neighbor_rotations(self.0, right, 0)?;
                let right_cell = Self(right_bits);

                if right_cell < owner {
                    owner = right_cell;
                    let direction = if owner.is_pentagon() {
                        grid::direction_for_neighbor(owner.0, self.0)?
                    } else {
                        let offset = (REV_NEIGHBOR_DIRECTIONS_HEX[usize::from(right) - 1]
                            + right_rotation)
                            % NUM_HEX_VERTS;
                        DIRECTIONS_CCW[usize::from(offset)]
                    };
                    owner_vertex = owner.vertex_for_direction(direction)?;
                }
            }

            if owner == left_cell {
                let direction = if owner.is_pentagon() {
                    grid::direction_for_neighbor(owner.0, self.0)?
                } else {
                    let offset = (REV_NEIGHBOR_DIRECTIONS_HEX[usize::from(left) - 1]
                        + left_rotation)
                        % NUM_HEX_VERTS;
                    DIRECTIONS_CCW[usize::from(offset)]
                };
                owner_vertex = (owner.vertex_for_direction(direction)? + 1) % owner.vertex_count();
            }
        }

        Some(H3Vertex(set_reserved_bits(
            set_mode(owner.0, VERTEX_MODE),
            owner_vertex,
        )))
    }

    /// All canonical vertex indexes for this cell.
    pub fn vertexes(self) -> Vec<H3Vertex> {
        (0..self.vertex_count())
            .filter_map(|vertex| self.vertex(vertex))
            .collect()
    }

    fn vertex_rotations(self) -> u8 {
        let fijk =
            faceijk::FaceIJK::from_cell(self.0, self.base_cell(), self.resolution() as usize);
        let mut rotations =
            faceijk::FaceIJK::base_cell_rotation_count(self.base_cell(), fijk.face).unwrap_or(0);

        if base_cells::is_pentagon(self.base_cell()) {
            let direction_faces = base_cells::pentagon_direction_faces(self.base_cell());
            let home_face = base_cells::home_face(self.base_cell());

            if fijk.face != home_face
                && (base_cells::is_polar_pentagon(self.base_cell())
                    || fijk.face == direction_faces[3])
            {
                rotations = (rotations + 1) % NUM_HEX_VERTS;
            }

            if first_axe(self.0) == Some(3) && fijk.face == direction_faces[3] {
                rotations = (rotations + 5) % NUM_HEX_VERTS;
            } else if first_axe(self.0) == Some(5) && fijk.face == direction_faces[1] {
                rotations = (rotations + 1) % NUM_HEX_VERTS;
            }
        }

        rotations
    }

    fn direction_for_vertex(self, vertex: u8) -> Option<u8> {
        if vertex >= self.vertex_count() {
            return None;
        }
        let rotations = self.vertex_rotations();
        if self.is_pentagon() {
            let index = (vertex + rotations) % NUM_PENT_VERTS;
            Some(VERTEX_TO_DIRECTION_PENT[usize::from(index)])
        } else {
            let index = (vertex + rotations) % NUM_HEX_VERTS;
            Some(VERTEX_TO_DIRECTION_HEX[usize::from(index)])
        }
    }

    fn vertex_for_direction(self, direction: u8) -> Option<u8> {
        if direction == 0 || direction > 6 || (self.is_pentagon() && direction == K_AXIS_DIGIT) {
            return None;
        }
        let rotations = self.vertex_rotations();
        if self.is_pentagon() {
            let vertex = DIRECTION_TO_VERTEX_PENT[usize::from(direction - 2)];
            Some((vertex + NUM_PENT_VERTS - rotations) % NUM_PENT_VERTS)
        } else {
            let vertex = DIRECTION_TO_VERTEX_HEX[usize::from(direction - 1)];
            Some((vertex + NUM_HEX_VERTS - rotations) % NUM_HEX_VERTS)
        }
    }

    /// Cell boundary as geographic `(latitude, longitude)` points in degrees.
    pub fn boundary_lat_lng(self) -> Vec<(f64, f64)> {
        faceijk::FaceIJK::cell_boundary_lat_lng_rad(
            self.0,
            self.base_cell(),
            self.resolution() as usize,
            self.is_pentagon(),
        )
        .into_iter()
        .map(|(lat, lng)| (lat.to_degrees(), lng.to_degrees()))
        .collect()
    }

    /// Exact spherical cell area in radians squared.
    pub fn area_rads2(self) -> f64 {
        linear_ring_area_rads2(&faceijk::FaceIJK::cell_boundary_lat_lng_rad(
            self.0,
            self.base_cell(),
            self.resolution() as usize,
            self.is_pentagon(),
        ))
    }

    /// Exact spherical cell area in square kilometers.
    pub fn area_km2(self) -> f64 {
        self.area_rads2() * EARTH_RADIUS_KM * EARTH_RADIUS_KM
    }

    /// Exact spherical cell area in square meters.
    pub fn area_m2(self) -> f64 {
        self.area_km2() * 1_000_000.0
    }

    /// Maximum number of cells in a grid disk of radius `k`.
    pub fn max_grid_disk_size(k: u32) -> u64 {
        let k = u64::from(k);
        3 * k * (k + 1) + 1
    }

    /// Average H3 hexagon area in square kilometers for `resolution`.
    pub fn average_hexagon_area_km2(resolution: u8) -> Option<f64> {
        AVG_HEXAGON_AREA_KM2.get(usize::from(resolution)).copied()
    }

    /// Average H3 hexagon area in square meters for `resolution`.
    pub fn average_hexagon_area_m2(resolution: u8) -> Option<f64> {
        AVG_HEXAGON_AREA_M2.get(usize::from(resolution)).copied()
    }

    /// Average H3 hexagon edge length in kilometers for `resolution`.
    pub fn average_hexagon_edge_length_km(resolution: u8) -> Option<f64> {
        AVG_HEXAGON_EDGE_LENGTH_KM
            .get(usize::from(resolution))
            .copied()
    }

    /// Average H3 hexagon edge length in meters for `resolution`.
    pub fn average_hexagon_edge_length_m(resolution: u8) -> Option<f64> {
        AVG_HEXAGON_EDGE_LENGTH_M
            .get(usize::from(resolution))
            .copied()
    }

    /// Number of H3 cells at `resolution`.
    pub fn num_cells(resolution: u8) -> Option<u64> {
        (resolution <= MAX_RESOLUTION).then(|| 2 + 120 * 7u64.pow(u32::from(resolution)))
    }

    /// Compresses a unique, same-resolution cell set by replacing full child
    /// branches with their parent cells.
    pub fn compact_cells(cells: &mut Vec<Self>) -> Result<(), H3CompactionError> {
        let Some(first) = cells.first() else {
            return Ok(());
        };
        let resolution = first.resolution();
        if cells.iter().any(|cell| cell.resolution() != resolution) {
            return Err(H3CompactionError::HeterogeneousResolution);
        }

        cells.sort_unstable();
        if cells.windows(2).any(|window| window[0] == window[1]) {
            return Err(H3CompactionError::DuplicateInput);
        }

        let mut set = cells.iter().copied().collect::<BTreeSet<_>>();
        for res in (1..=resolution).rev() {
            let parents = set
                .iter()
                .filter(|cell| cell.resolution() == res)
                .filter_map(|cell| cell.parent(res - 1))
                .collect::<BTreeSet<_>>();

            for parent in parents {
                let children = parent.children_at_resolution(res);
                if children.iter().all(|child| set.contains(child)) {
                    for child in children {
                        set.remove(&child);
                    }
                    set.insert(parent);
                }
            }
        }

        *cells = set.into_iter().collect();
        Ok(())
    }

    /// Exact output size for [`H3Cell::uncompact_cells`].
    pub fn uncompact_size(
        compacted: impl IntoIterator<Item = Self>,
        resolution: u8,
    ) -> Option<u64> {
        if resolution > MAX_RESOLUTION {
            return None;
        }
        let mut size = 0u64;
        for cell in compacted {
            if cell.resolution() > resolution {
                return None;
            }
            size = size.checked_add(cell.children_count(resolution))?;
        }
        Some(size)
    }

    /// Expands compacted cells to the requested `resolution`.
    pub fn uncompact_cells(
        compacted: impl IntoIterator<Item = Self>,
        resolution: u8,
    ) -> Option<Vec<Self>> {
        let compacted = compacted.into_iter().collect::<Vec<_>>();
        let size = Self::uncompact_size(compacted.iter().copied(), resolution)?;
        let mut out = Vec::with_capacity(usize::try_from(size).ok()?);
        for cell in compacted {
            if cell.resolution() > resolution {
                return None;
            }
            out.extend(cell.children_at_resolution(resolution));
        }
        Some(out)
    }
}

/// A cell coordinate in another cell's local IJ coordinate system.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct H3LocalIJ {
    /// Anchor/origin cell for this local coordinate system.
    pub anchor: H3Cell,
    /// Local I coordinate.
    pub i: i32,
    /// Local J coordinate.
    pub j: i32,
}

/// An H3 directed edge index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct H3DirectedEdge(u64);

impl H3DirectedEdge {
    /// Wraps a raw u64 H3 directed-edge index, validating it is well-formed.
    pub fn try_new(bits: u64) -> Option<Self> {
        is_valid_directed_edge(bits).then_some(Self(bits))
    }

    /// True iff `bits` is a valid H3 directed-edge index.
    pub fn is_valid(bits: u64) -> bool {
        is_valid_directed_edge(bits)
    }

    /// Constructs a directed edge from neighboring cells.
    pub fn new(origin: H3Cell, destination: H3Cell) -> Option<Self> {
        origin.edge(destination)
    }

    /// Constructs the edge from `origin` in H3 direction digit `1..=6`.
    pub fn from_origin_and_direction(origin: H3Cell, direction: u8) -> Option<Self> {
        if !(1..=6).contains(&direction) || (origin.is_pentagon() && direction == K_AXIS_DIGIT) {
            return None;
        }
        let bits = set_reserved_bits(set_mode(origin.0, DIRECTED_EDGE_MODE), direction);
        Self::try_new(bits)
    }

    /// The raw u64 H3 directed-edge index.
    pub const fn as_u64(self) -> u64 {
        self.0
    }

    /// The H3 direction digit (`1..=6`) encoded by this edge.
    pub const fn direction(self) -> u8 {
        get_reserved_bits(self.0)
    }

    /// Origin cell.
    pub fn origin(self) -> H3Cell {
        H3Cell(set_reserved_bits(set_mode(self.0, CELL_MODE), 0))
    }

    /// Destination cell.
    pub fn destination(self) -> H3Cell {
        let (neighbor, _) =
            grid::neighbor_rotations(self.origin().0, self.direction(), 0).expect("valid edge");
        H3Cell(neighbor)
    }

    /// Origin/destination cells.
    pub fn cells(self) -> (H3Cell, H3Cell) {
        (self.origin(), self.destination())
    }

    /// Canonical topological start/end vertices for this directed edge.
    pub fn vertices(self) -> (H3Vertex, H3Vertex) {
        let origin = self.origin();
        let start = origin
            .vertex_for_direction(self.direction())
            .expect("valid edge direction");
        let end = (start + 1) % origin.vertex_count();
        (
            origin.vertex(start).expect("valid start vertex"),
            origin.vertex(end).expect("valid end vertex"),
        )
    }

    /// Directed-edge boundary as geographic `(latitude, longitude)` points in
    /// degrees. Edges that cross an icosahedron face include the H3 distortion
    /// point, so the result has two or three points.
    pub fn boundary_lat_lng(self) -> Vec<(f64, f64)> {
        self.boundary_lat_lng_rad()
            .into_iter()
            .map(|(lat, lng)| (lat.to_degrees(), lng.to_degrees()))
            .collect()
    }

    /// Directed-edge length in radians.
    pub fn length_rads(self) -> f64 {
        let boundary = self.boundary_lat_lng_rad();
        boundary
            .windows(2)
            .map(|points| distance_rads(points[0], points[1]))
            .sum()
    }

    /// Directed-edge length in kilometers.
    pub fn length_km(self) -> f64 {
        self.length_rads() * EARTH_RADIUS_KM
    }

    /// Directed-edge length in meters.
    pub fn length_m(self) -> f64 {
        self.length_km() * 1000.0
    }

    /// Directed edge with origin and destination reversed.
    pub fn reverse(self) -> Self {
        let (origin, destination) = self.cells();
        destination.edge(origin).expect("valid reverse edge")
    }

    fn boundary_lat_lng_rad(self) -> Vec<(f64, f64)> {
        let origin = self.origin();
        let start = origin
            .vertex_for_direction(self.direction())
            .expect("valid edge direction");
        faceijk::FaceIJK::cell_boundary_segment_lat_lng_rad(
            origin.0,
            origin.base_cell(),
            origin.resolution() as usize,
            origin.is_pentagon(),
            start,
            2,
        )
    }
}

/// An H3 vertex index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct H3Vertex(u64);

impl H3Vertex {
    /// Wraps a raw u64 H3 vertex index, validating it is canonical.
    pub fn try_new(bits: u64) -> Option<Self> {
        is_valid_vertex(bits).then_some(Self(bits))
    }

    /// True iff `bits` is a valid canonical H3 vertex index.
    pub fn is_valid(bits: u64) -> bool {
        is_valid_vertex(bits)
    }

    /// Constructs the canonical vertex index for `owner` and topological
    /// `vertex`.
    pub fn new(owner: H3Cell, vertex: u8) -> Option<Self> {
        owner.vertex(vertex)
    }

    /// The raw u64 H3 vertex index.
    pub const fn as_u64(self) -> u64 {
        self.0
    }

    /// Topological vertex number on the canonical owner cell.
    pub const fn vertex(self) -> u8 {
        get_reserved_bits(self.0)
    }

    /// Canonical owner cell.
    pub fn owner(self) -> H3Cell {
        H3Cell(set_reserved_bits(set_mode(self.0, CELL_MODE), 0))
    }

    /// Vertex coordinate as `(latitude, longitude)` in degrees.
    pub fn lat_lng(self) -> (f64, f64) {
        let (lat, lng) = self.lat_lng_rad();
        (lat.to_degrees(), lng.to_degrees())
    }

    fn lat_lng_rad(self) -> (f64, f64) {
        let owner = self.owner();
        faceijk::FaceIJK::cell_boundary_segment_lat_lng_rad(
            owner.0,
            owner.base_cell(),
            owner.resolution() as usize,
            owner.is_pentagon(),
            self.vertex(),
            1,
        )
        .into_iter()
        .next()
        .expect("vertex boundary point")
    }
}

/// Error returned when compacting an invalid input set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum H3CompactionError {
    /// Every input cell must have the same resolution.
    HeterogeneousResolution,
    /// Input cells must be unique.
    DuplicateInput,
}

impl fmt::Display for H3CompactionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::HeterogeneousResolution => f.write_str("heterogeneous H3 cell resolutions"),
            Self::DuplicateInput => f.write_str("duplicate H3 cell input"),
        }
    }
}

impl std::error::Error for H3CompactionError {}

// ---- raw-bits accessors (shared with later phases) -------------------------

#[inline]
pub(super) const fn get_mode(bits: u64) -> u8 {
    ((bits >> MODE_OFFSET) & 0xf) as u8
}

#[inline]
pub(super) const fn get_resolution(bits: u64) -> u8 {
    ((bits >> RESOLUTION_OFFSET) & 0xf) as u8
}

#[inline]
pub(super) const fn get_base_cell(bits: u64) -> u8 {
    ((bits >> BASE_CELL_OFFSET) & 0x7f) as u8
}

#[inline]
pub(super) const fn get_reserved_bits(bits: u64) -> u8 {
    ((bits >> RESERVED_OFFSET) & 0x7) as u8
}

/// Returns the 3-bit digit at resolution `res` (1..=15).
#[inline]
pub(super) const fn get_digit(bits: u64, res: u8) -> u8 {
    debug_assert!(res >= 1 && res <= MAX_RESOLUTION);
    ((bits >> ((MAX_RESOLUTION - res) as u64 * 3)) & 0x7) as u8
}

/// A cell index with mode=cell, resolution 0, base cell 0, and all 15 digits
/// set to the sentinel `7` — the starting point for building an index.
pub(super) const DEFAULT_CELL_INDEX: u64 =
    ((CELL_MODE) << MODE_OFFSET) | ((1u64 << (3 * MAX_RESOLUTION as u64)) - 1);

#[inline]
const fn digit_shift(res: u8) -> u64 {
    (MAX_RESOLUTION - res) as u64 * 3
}

/// `get_digit`, named for the projection code.
#[inline]
pub(super) const fn get_direction(bits: u64, res: u8) -> u8 {
    get_digit(bits, res)
}

#[inline]
pub(super) const fn set_resolution(bits: u64, res: u8) -> u64 {
    (bits & !(0xf << RESOLUTION_OFFSET)) | ((res as u64) << RESOLUTION_OFFSET)
}

#[inline]
const fn set_mode(bits: u64, mode: u64) -> u64 {
    (bits & !(0xf << MODE_OFFSET)) | (mode << MODE_OFFSET)
}

#[inline]
const fn set_reserved_bits(bits: u64, value: u8) -> u64 {
    (bits & !(0x7 << RESERVED_OFFSET)) | ((value as u64) << RESERVED_OFFSET)
}

#[inline]
pub(super) const fn set_base_cell(bits: u64, base_cell: u8) -> u64 {
    (bits & !(0x7f << BASE_CELL_OFFSET)) | ((base_cell as u64) << BASE_CELL_OFFSET)
}

#[inline]
pub(super) const fn set_direction(bits: u64, digit: u8, res: u8) -> u64 {
    let shift = digit_shift(res);
    (bits & !(0x7 << shift)) | ((digit as u64) << shift)
}

#[inline]
pub(super) const fn clr_direction(bits: u64, res: u8) -> u64 {
    bits & !(0x7 << digit_shift(res))
}

/// Rotates a single direction digit `count` times by 60°.
#[inline]
pub(super) const fn rotate_direction(digit: u8, count: usize, ccw: bool) -> u8 {
    rotate_digit(digit, count, ccw)
}

/// Sets all unused (beyond-`res`) digit slots to the sentinel `7`.
#[inline]
pub(super) const fn set_unused(bits: u64, res: u8) -> u64 {
    bits | ((1u64 << digit_shift(res)) - 1)
}

/// The first non-center digit (1..=6) reading from resolution 1, or `None`.
pub(super) const fn first_axe(bits: u64) -> Option<u8> {
    let res = get_resolution(bits);
    let mut r = 1u8;
    while r <= res {
        let d = get_digit(bits, r);
        if d != 0 {
            return Some(d);
        }
        r += 1;
    }
    None
}

/// A single 60° digit rotation (`rotate60_once`).
#[inline]
const fn rotate_digit_once(digit: u8, ccw: bool) -> u8 {
    match digit {
        1 => {
            if ccw {
                5
            } else {
                3
            }
        } // K -> IK / JK
        2 => {
            if ccw {
                3
            } else {
                6
            }
        } // J -> JK / IJ
        3 => {
            if ccw {
                1
            } else {
                2
            }
        } // JK -> K / J
        4 => {
            if ccw {
                6
            } else {
                5
            }
        } // I -> IJ / IK
        5 => {
            if ccw {
                4
            } else {
                1
            }
        } // IK -> I / K
        6 => {
            if ccw {
                2
            } else {
                4
            }
        } // IJ -> J / I
        _ => 0, // Center
    }
}

/// `count` 60° digit rotations.
const fn rotate_digit(digit: u8, count: usize, ccw: bool) -> u8 {
    const CCW_SEQ: [u8; 6] = [1, 5, 4, 6, 2, 3];
    const CW_SEQ: [u8; 6] = [1, 3, 2, 6, 4, 5];
    if count == 0 || digit == 0 {
        return digit;
    }
    let offset = match digit {
        1 => 0,
        2 => {
            if ccw {
                4
            } else {
                2
            }
        }
        3 => {
            if ccw {
                5
            } else {
                1
            }
        }
        4 => {
            if ccw {
                2
            } else {
                4
            }
        }
        5 => {
            if ccw {
                1
            } else {
                5
            }
        }
        6 => 3,
        _ => return digit,
    };
    let index = (count + offset) % 6;
    if ccw { CCW_SEQ[index] } else { CW_SEQ[index] }
}

/// Rotates a whole index `count` times by 60°.
pub(super) fn rotate60(mut bits: u64, count: usize, ccw: bool) -> u64 {
    let res = get_resolution(bits);
    let mut r = 1u8;
    while r <= res {
        let d = get_direction(bits, r);
        let d = if count == 1 {
            rotate_digit_once(d, ccw)
        } else {
            rotate_digit(d, count, ccw)
        };
        bits = set_direction(bits, d, r);
        r += 1;
    }
    bits
}

/// Rotates a pentagon index by 60°, accounting for the missing K-axis sequence.
pub(super) fn pentagon_rotate60(mut bits: u64, ccw: bool) -> u64 {
    let res = get_resolution(bits);
    if res == 0 {
        return bits;
    }
    // Direction that would rotate to the (deleted) K axis.
    let target = if ccw { 3 } else { 5 }; // JK.axe() / IK.axe()
    let double = first_axe(bits) == Some(target);
    let mut r = 1u8;
    while r <= res {
        let d = get_direction(bits, r);
        let d = if double {
            rotate_digit(d, 2, ccw)
        } else {
            rotate_digit_once(d, ccw)
        };
        bits = set_direction(bits, d, r);
        r += 1;
    }
    bits
}

/// True iff any used digit (1..=res) is non-center.
const fn has_nonzero_digit(bits: u64) -> bool {
    let res = get_resolution(bits);
    let mut r = 1u8;
    while r <= res {
        if get_digit(bits, r) != 0 {
            return true;
        }
        r += 1;
    }
    false
}

/// Mirrors H3 `isValidCell`: high bit clear, cell mode, reserved bits clear,
/// base cell in range, used digits 0..=6 (no `7`), the pentagon deleted-K-axis
/// rule, and unused digits all `7`.
fn is_valid_cell(bits: u64) -> bool {
    if (bits >> 63) & 1 != 0 {
        return false;
    }
    if get_mode(bits) as u64 != CELL_MODE {
        return false;
    }
    if (bits >> RESERVED_OFFSET) & 0x7 != 0 {
        return false;
    }
    let base_cell = get_base_cell(bits);
    if base_cell >= base_cells::COUNT {
        return false;
    }
    let res = get_resolution(bits); // 4 bits -> always 0..=15

    let is_pent_base = base_cells::is_pentagon(base_cell);
    let mut found_first_nonzero = false;
    let mut r = 1u8;
    while r <= res {
        let d = get_digit(bits, r);
        if !found_first_nonzero && d != 0 {
            found_first_nonzero = true;
            // Pentagons have no K-axis child (deleted subsequence).
            if is_pent_base && d == K_AXIS_DIGIT {
                return false;
            }
        }
        if d == INVALID_DIGIT {
            return false;
        }
        r += 1;
    }
    // Unused digits must be the sentinel 7.
    while r <= MAX_RESOLUTION {
        if get_digit(bits, r) != INVALID_DIGIT {
            return false;
        }
        r += 1;
    }
    true
}

fn is_valid_directed_edge(bits: u64) -> bool {
    if get_mode(bits) as u64 != DIRECTED_EDGE_MODE {
        return false;
    }
    let direction = get_reserved_bits(bits);
    if !(1..=6).contains(&direction) {
        return false;
    }

    let origin = set_reserved_bits(set_mode(bits, CELL_MODE), 0);
    let Some(origin) = H3Cell::try_new(origin) else {
        return false;
    };

    !(origin.is_pentagon() && direction == K_AXIS_DIGIT)
}

fn is_valid_vertex(bits: u64) -> bool {
    if get_mode(bits) as u64 != VERTEX_MODE {
        return false;
    }

    let owner_bits = set_reserved_bits(set_mode(bits, CELL_MODE), 0);
    let Some(owner) = H3Cell::try_new(owner_bits) else {
        return false;
    };
    let vertex = get_reserved_bits(bits);
    owner.vertex(vertex).map(H3Vertex::as_u64) == Some(bits)
}

fn distance_rads(a: (f64, f64), b: (f64, f64)) -> f64 {
    let sin_lat = ((b.0 - a.0) * 0.5).sin();
    let sin_lng = ((b.1 - a.1) * 0.5).sin();
    let value = sin_lat.mul_add(sin_lat, a.0.cos() * b.0.cos() * sin_lng * sin_lng);
    2.0 * value.sqrt().atan2((1.0 - value).sqrt())
}

fn linear_ring_area_rads2(ring: &[(f64, f64)]) -> f64 {
    let mut adder = FloatAdder::default();
    for (a, b) in ring
        .iter()
        .zip(ring.iter().cycle().skip(1))
        .take(ring.len())
    {
        // Boundary tuples are `(lat, lng)`; Cagnoli's formula takes `(lng, lat)`.
        adder.add(cagnoli((a.1, a.0), (b.1, b.0)));
    }

    let area = adder.value();
    if area < 0.0 {
        area + 4.0 * std::f64::consts::PI
    } else {
        area
    }
}

fn cagnoli(a: (f64, f64), b: (f64, f64)) -> f64 {
    let a_lat = a.1.mul_add(0.5, std::f64::consts::PI * 0.25);
    let b_lat = b.1.mul_add(0.5, std::f64::consts::PI * 0.25);
    let sin_a = a_lat.sin() * b_lat.sin();
    let cos_a = a_lat.cos() * b_lat.cos();

    let delta = b.0 - a.0;
    -2.0 * (sin_a * delta.sin()).atan2(sin_a * delta.cos() + cos_a)
}

#[derive(Default)]
struct FloatAdder {
    sum: f64,
    correction: f64,
}

impl FloatAdder {
    fn add(&mut self, value: f64) {
        let y = value - self.correction;
        let t = self.sum + y;
        self.correction = (t - self.sum) - y;
        self.sum = t;
    }

    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

impl fmt::Display for H3Cell {
    /// Lowercase hex with no leading zeros — matches H3 `h3ToString`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:x}", self.0)
    }
}

impl fmt::Display for H3DirectedEdge {
    /// Lowercase hex with no leading zeros — matches H3 `h3ToString`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:x}", self.0)
    }
}

impl fmt::Display for H3Vertex {
    /// Lowercase hex with no leading zeros — matches H3 `h3ToString`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:x}", self.0)
    }
}

/// Error parsing or validating an H3 cell token.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvalidH3Cell;

impl fmt::Display for InvalidH3Cell {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid H3 cell index")
    }
}

impl std::error::Error for InvalidH3Cell {}

impl FromStr for H3Cell {
    type Err = InvalidH3Cell;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bits = u64::from_str_radix(s.trim(), 16).map_err(|_| InvalidH3Cell)?;
        Self::try_new(bits).ok_or(InvalidH3Cell)
    }
}

/// Error parsing or validating an H3 directed-edge token.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvalidH3DirectedEdge;

impl fmt::Display for InvalidH3DirectedEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid H3 directed edge index")
    }
}

impl std::error::Error for InvalidH3DirectedEdge {}

impl FromStr for H3DirectedEdge {
    type Err = InvalidH3DirectedEdge;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bits = u64::from_str_radix(s.trim(), 16).map_err(|_| InvalidH3DirectedEdge)?;
        Self::try_new(bits).ok_or(InvalidH3DirectedEdge)
    }
}

/// Error parsing or validating an H3 vertex token.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvalidH3Vertex;

impl fmt::Display for InvalidH3Vertex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid H3 vertex index")
    }
}

impl std::error::Error for InvalidH3Vertex {}

impl FromStr for H3Vertex {
    type Err = InvalidH3Vertex;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bits = u64::from_str_radix(s.trim(), 16).map_err(|_| InvalidH3Vertex)?;
        Self::try_new(bits).ok_or(InvalidH3Vertex)
    }
}
