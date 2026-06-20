//! The public H3 index types — [`H3Cell`], [`H3DirectedEdge`], [`H3Vertex`],
//! and [`H3LocalIJ`] — with their parsing / `Display` and the operations the
//! workspace uses: parent/children, child positions, compaction, `grid_disk`,
//! `grid_distance`, local IJ, grid paths, edges, vertices, and exact area.
//!
//! The supporting machinery lives in submodules: the u64 [`bits`] codec
//! (shared across the whole `h3` module), the index [`validity`] rules behind
//! `try_new`, and the spherical [`area`] / edge-length math.
//!
//! [`H3Cell`]'s own operations are grouped into sibling submodules —
//! construction/projection ([`cell_core`]), the hierarchy ([`cell_hierarchy`]),
//! grid traversal ([`cell_traversal`]), edges/vertices ([`cell_topology`]), and
//! geometric/tabular measures ([`cell_measures`]) — each adding an `impl
//! H3Cell` block.

use std::fmt;
use std::str::FromStr;

use oresme::LengthMeters;

use super::{base_cells, children, coord, faceijk, grid};

mod area;
mod bits;
mod cell_core;
mod cell_hierarchy;
mod cell_measures;
mod cell_topology;
mod cell_traversal;
mod validity;

use area::{
    AVG_HEXAGON_AREA_KM2, AVG_HEXAGON_AREA_M2, AVG_HEXAGON_EDGE_LENGTH_KM,
    AVG_HEXAGON_EDGE_LENGTH_M, EARTH_RADIUS_KM, distance_rads, linear_ring_area_rads2,
};
pub use bits::MAX_RESOLUTION;
use bits::{
    CELL_MODE, DIRECTED_EDGE_MODE, K_AXIS_DIGIT, VERTEX_MODE, get_digit, get_reserved_bits,
    has_nonzero_digit, set_mode, set_reserved_bits, set_unused,
};
pub(super) use bits::{
    DEFAULT_CELL_INDEX, clr_direction, first_axe, get_base_cell, get_direction, get_resolution,
    pentagon_rotate60, rotate_direction, rotate60, set_base_cell, set_direction, set_resolution,
};
use validity::{is_valid_cell, is_valid_directed_edge, is_valid_vertex};

/// An H3 cell index.
///
/// Always a valid cell (mode 1) — construct via [`H3Cell::try_new`] /
/// [`H3Cell::from_str`], both of which reject non-cell or malformed indices the
/// same way the H3 reference does.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct H3Cell(u64);

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

impl H3LocalIJ {
    /// The cell this local IJ coordinate denotes, at the anchor's resolution —
    /// the inverse of [`H3Cell::to_local_ij`] (mirrors H3 `localIjToCell`).
    /// Returns `None` when the coordinate is outside the anchor's local frame
    /// (too far away or across a pentagon's deleted subsequence).
    pub fn to_cell(self) -> Option<H3Cell> {
        let ijk = coord::CoordIJK::new(self.i, self.j, 0).checked_normalize()?;
        grid::cell_from_local_ijk(self.anchor.0, ijk).and_then(H3Cell::try_new)
    }
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
