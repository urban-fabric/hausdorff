//! The public H3 index types — [`H3Cell`], [`H3DirectedEdge`], [`H3Vertex`],
//! and [`H3LocalIJ`] — with their parsing / `Display` and the operations the
//! workspace uses: parent/children, child positions, compaction, `grid_disk`,
//! `grid_distance`, local IJ, grid paths, edges, vertices, and exact area.
//!
//! The supporting machinery lives in submodules: the u64 [`bits`] codec
//! (shared across the whole `h3` module), the index [`validity`] rules behind
//! `try_new`, and the spherical [`area`] / edge-length math.

use std::collections::BTreeSet;
use std::fmt;
use std::str::FromStr;

use euclid::{AreaSquareMeters, LengthMeters};

use super::{base_cells, children, coord, faceijk, grid};

mod area;
mod bits;
mod validity;

pub use bits::MAX_RESOLUTION;
// Codec helpers the sibling h3 modules reach via the `index::` path.
use area::{
    AVG_HEXAGON_AREA_KM2, AVG_HEXAGON_AREA_M2, AVG_HEXAGON_EDGE_LENGTH_KM,
    AVG_HEXAGON_EDGE_LENGTH_M, EARTH_RADIUS_KM, distance_rads, linear_ring_area_rads2,
};
use bits::{
    CELL_MODE, DIRECTED_EDGE_MODE, K_AXIS_DIGIT, VERTEX_MODE, get_digit, get_reserved_bits,
    has_nonzero_digit, set_mode, set_reserved_bits, set_unused,
};
pub(super) use bits::{
    DEFAULT_CELL_INDEX, clr_direction, first_axe, get_base_cell, get_direction, get_resolution,
    pentagon_rotate60, rotate_direction, rotate60, set_base_cell, set_direction, set_resolution,
};
use validity::{is_valid_cell, is_valid_directed_edge, is_valid_vertex};

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

    /// The cell at `resolution` containing the given [`GeoPoint`](crate::earth::GeoPoint).
    /// Convenience over [`for_point`](Self::for_point).
    pub fn for_geo_point(point: crate::earth::GeoPoint, resolution: u8) -> Option<Self> {
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

    /// Exact spherical cell area as a typed [`AreaSquareMeters`] — prefer this
    /// over [`area_m2`](Self::area_m2) so the Earth-area unit travels with the
    /// value and can't be confused for a bare scalar.
    pub fn area(self) -> AreaSquareMeters {
        AreaSquareMeters::new(self.area_m2())
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

    /// Average H3 hexagon area for `resolution` as a typed [`AreaSquareMeters`] —
    /// the unit-carrying form of [`average_hexagon_area_m2`](Self::average_hexagon_area_m2).
    pub fn average_hexagon_area(resolution: u8) -> Option<AreaSquareMeters> {
        Self::average_hexagon_area_m2(resolution).map(AreaSquareMeters::new)
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

    /// Average H3 hexagon edge length for `resolution` as a typed [`LengthMeters`]
    /// — the unit-carrying form of
    /// [`average_hexagon_edge_length_m`](Self::average_hexagon_edge_length_m).
    pub fn average_hexagon_edge_length(resolution: u8) -> Option<LengthMeters> {
        Self::average_hexagon_edge_length_m(resolution).map(LengthMeters::new)
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

    /// Directed-edge length as a typed [`LengthMeters`] — prefer this over
    /// [`length_m`](Self::length_m) so the Earth-distance unit travels with the
    /// value (and isn't mistaken for the radian measure [`length_rads`](Self::length_rads)).
    pub fn length(self) -> LengthMeters {
        LengthMeters::new(self.length_m())
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
