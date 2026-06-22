//! Lightweight polygon operations for H3 cells.

use std::collections::{BTreeSet, VecDeque};
use std::fmt;

use super::index::{H3Cell, MAX_RESOLUTION};

mod boundary;
mod ring;

pub use boundary::H3BoundaryRing;
use ring::{
    BBox, normalize_ring, ring_contains_point, ring_edges, ring_is_transmeridian,
    segments_cross_properly, segments_intersect, shift_ring_east,
};

/// A geographic polygon in `(latitude, longitude)` degrees.
///
/// Rings are stored without a duplicated closing point. A polygon whose
/// exterior ring crosses the antimeridian (some edge spans more than 180° of
/// longitude, so the short way runs through ±180) is detected at construction
/// and handled by shifting all ring longitudes into a continuous `(0°, 360°)`
/// frame; query points are shifted the same way, so containment and covers
/// work across the date line.
#[derive(Clone, Debug, PartialEq)]
pub struct H3Polygon {
    exterior: Vec<(f64, f64)>,
    holes: Vec<Vec<(f64, f64)>>,
    transmeridian: bool,
    bbox: BBox,
}

/// How a cell must relate to a polygon to be part of its cell cover — mirrors
/// the H3 v4 / h3o polygon-to-cells containment modes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ContainmentMode {
    /// The cell's center point is inside the polygon. Cells are uniquely
    /// assigned (adjacent polygons get disjoint covers), but the cover can
    /// both overshoot and undershoot the polygon area.
    Centroid,
    /// The cell is fully contained inside the polygon. Never overshoots, but
    /// leaves boundary-straddling area uncovered.
    Full,
    /// The cell intersects the polygon at all, even partially — a guaranteed
    /// superset of the polygon area (no false negatives), at the cost of
    /// overshooting and of adjacent polygons sharing boundary cells.
    Overlapping,
}

impl H3Polygon {
    /// Constructs a polygon from an exterior ring and zero or more hole rings.
    pub fn new(
        exterior: impl Into<Vec<(f64, f64)>>,
        holes: impl Into<Vec<Vec<(f64, f64)>>>,
    ) -> Result<Self, H3PolygonError> {
        let mut exterior = normalize_ring(exterior.into())?;
        let mut holes = holes
            .into()
            .into_iter()
            .map(normalize_ring)
            .collect::<Result<Vec<_>, _>>()?;

        let transmeridian = ring_is_transmeridian(&exterior);
        if transmeridian {
            shift_ring_east(&mut exterior);
            for hole in &mut holes {
                shift_ring_east(hole);
            }
        }
        let bbox = BBox::from_rings(std::iter::once(&exterior).chain(holes.iter()))?;

        Ok(Self {
            exterior,
            holes,
            transmeridian,
            bbox,
        })
    }

    /// Exterior ring, without a duplicated closing point. For an
    /// antimeridian-crossing polygon the longitudes are in the continuous
    /// `(0°, 360°)` frame (west-hemisphere longitudes shifted by `+360`).
    pub fn exterior(&self) -> &[(f64, f64)] {
        &self.exterior
    }

    /// Hole rings, without duplicated closing points. Longitudes are framed as
    /// in [`exterior`](Self::exterior).
    pub fn holes(&self) -> &[Vec<(f64, f64)>] {
        &self.holes
    }

    /// True iff this polygon crosses the antimeridian.
    pub fn is_transmeridian(&self) -> bool {
        self.transmeridian
    }

    /// True iff the point `(latitude, longitude)` in degrees is inside the
    /// polygon and outside all holes.
    pub fn contains_lat_lng(&self, lat: f64, lng: f64) -> bool {
        self.contains_framed((lat, self.frame_lng(lng)))
    }

    /// H3 cells at `resolution` whose centers are inside this polygon —
    /// [`to_cells_with_mode`](Self::to_cells_with_mode) in
    /// [`ContainmentMode::Centroid`].
    pub fn to_cells(&self, resolution: u8) -> Result<Vec<H3Cell>, H3PolygonError> {
        self.to_cells_with_mode(resolution, ContainmentMode::Centroid)
    }

    /// H3 cells at `resolution` related to this polygon as `mode` requires.
    pub fn to_cells_with_mode(
        &self,
        resolution: u8,
        mode: ContainmentMode,
    ) -> Result<Vec<H3Cell>, H3PolygonError> {
        if resolution > MAX_RESOLUTION {
            return Err(H3PolygonError::InvalidResolution(resolution));
        }

        let mut seen = BTreeSet::new();
        let mut queue = VecDeque::new();
        for seed in self.seed_cells(resolution) {
            queue.push_back(seed);
        }

        let mut out = Vec::new();
        while let Some(cell) = queue.pop_front() {
            if !seen.insert(cell) {
                continue;
            }

            let boundary = self.framed_cell_boundary(cell);
            let Some(cell_bbox) = BBox::from_points(boundary.iter().copied()) else {
                continue;
            };
            if !self.bbox.intersects(&cell_bbox) {
                continue;
            }

            let matched = match mode {
                ContainmentMode::Centroid => {
                    let (lat, lng) = cell.center_lat_lng();
                    self.contains_lat_lng(lat, lng)
                }
                ContainmentMode::Full => self.covers_cell(&boundary),
                ContainmentMode::Overlapping => self.overlaps_cell(&boundary),
            };
            if matched {
                out.push(cell);
            }

            for neighbor in cell.grid_ring(1) {
                if !seen.contains(&neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }

        out.sort_unstable();
        Ok(out)
    }

    /// Shifts a query longitude into this polygon's continuous frame.
    fn frame_lng(&self, lng: f64) -> f64 {
        if self.transmeridian && lng < 0.0 {
            lng + 360.0
        } else {
            lng
        }
    }

    /// Point-in-polygon for a point already in this polygon's frame.
    fn contains_framed(&self, point: (f64, f64)) -> bool {
        self.bbox.contains(point)
            && ring_contains_point(&self.exterior, point)
            && !self
                .holes
                .iter()
                .any(|hole| ring_contains_point(hole, point))
    }

    /// Cell boundary in `(latitude, longitude)` degrees, shifted into this
    /// polygon's frame.
    fn framed_cell_boundary(&self, cell: H3Cell) -> Vec<(f64, f64)> {
        cell.boundary_lat_lng()
            .into_iter()
            .map(|(lat, lng)| (lat, self.frame_lng(lng)))
            .collect()
    }

    /// All rings, exterior first.
    fn rings(&self) -> impl Iterator<Item = &[(f64, f64)]> {
        std::iter::once(self.exterior.as_slice()).chain(self.holes.iter().map(Vec::as_slice))
    }

    /// True iff the framed cell boundary polygon is fully contained in this
    /// polygon: every cell vertex is inside, no polygon edge properly crosses
    /// a cell edge, and no hole is nested inside the cell.
    fn covers_cell(&self, cell: &[(f64, f64)]) -> bool {
        cell.iter().all(|&vertex| self.contains_framed(vertex))
            && !self.rings().any(|ring| {
                ring_edges(ring).any(|(a, b)| {
                    ring_edges(cell).any(|(c, d)| segments_cross_properly(a, b, c, d))
                })
            })
            && !self
                .holes
                .iter()
                .any(|hole| hole.iter().any(|&vertex| ring_contains_point(cell, vertex)))
    }

    /// True iff the framed cell boundary polygon intersects this polygon at
    /// all: a cell vertex inside the polygon, a ring vertex inside the cell,
    /// or any edge intersection (the standard exact polygon-overlap test for
    /// plumb-line edges).
    fn overlaps_cell(&self, cell: &[(f64, f64)]) -> bool {
        cell.iter().any(|&vertex| self.contains_framed(vertex))
            || self
                .rings()
                .any(|ring| ring.iter().any(|&vertex| ring_contains_point(cell, vertex)))
            || self.rings().any(|ring| {
                ring_edges(ring)
                    .any(|(a, b)| ring_edges(cell).any(|(c, d)| segments_intersect(a, b, c, d)))
            })
    }

    fn seed_cells(&self, resolution: u8) -> Vec<H3Cell> {
        let mut seeds = Vec::new();
        let mut push_seed = |point: (f64, f64)| {
            let lng = if point.1 > 180.0 {
                point.1 - 360.0
            } else {
                point.1
            };
            if let Some(cell) = H3Cell::for_point(point.0, lng, resolution) {
                seeds.push(cell);
            }
        };

        push_seed(self.bbox.center());
        for &point in &self.exterior {
            push_seed(point);
        }
        for hole in &self.holes {
            for &point in hole {
                push_seed(point);
            }
        }
        seeds.sort_unstable();
        seeds.dedup();
        seeds
    }
}

impl H3Cell {
    /// Convenience wrapper for [`H3Polygon::to_cells`].
    pub fn polygon_to_cells(
        polygon: &H3Polygon,
        resolution: u8,
    ) -> Result<Vec<Self>, H3PolygonError> {
        polygon.to_cells(resolution)
    }

    /// Convenience wrapper for [`H3Polygon::to_cells_with_mode`].
    pub fn polygon_to_cells_with_mode(
        polygon: &H3Polygon,
        resolution: u8,
        mode: ContainmentMode,
    ) -> Result<Vec<Self>, H3PolygonError> {
        polygon.to_cells_with_mode(resolution, mode)
    }
}

/// Polygon-operation error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum H3PolygonError {
    /// A ring had fewer than three vertices after removing a duplicated close.
    RingTooShort,
    /// Coordinates must be finite, with latitude in `[-90, 90]` and longitude in
    /// `[-180, 180]`.
    InvalidCoordinate,
    /// Resolution must be in `0..=15`.
    InvalidResolution(u8),
    /// Cell dissolution requires a same-resolution input set.
    HeterogeneousResolution,
    /// Cell dissolution requires unique input cells.
    DuplicateInput,
    /// Boundary arcs did not form closed rings.
    OpenBoundary,
}

impl fmt::Display for H3PolygonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RingTooShort => f.write_str("polygon ring has fewer than three vertices"),
            Self::InvalidCoordinate => f.write_str("invalid polygon coordinate"),
            Self::InvalidResolution(resolution) => {
                write!(f, "invalid H3 resolution {resolution}")
            }
            Self::HeterogeneousResolution => f.write_str("heterogeneous H3 cell resolutions"),
            Self::DuplicateInput => f.write_str("duplicate H3 cell input"),
            Self::OpenBoundary => f.write_str("H3 boundary arcs did not close"),
        }
    }
}

impl std::error::Error for H3PolygonError {}
