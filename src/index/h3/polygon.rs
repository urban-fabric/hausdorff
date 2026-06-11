//! Lightweight polygon operations for H3 cells.

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;

use super::index::{H3Cell, H3DirectedEdge, H3Vertex, MAX_RESOLUTION};

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

    /// Dissolves a same-resolution cell set into external boundary rings.
    ///
    /// Ring points are `(latitude, longitude)` degrees and do not repeat the
    /// first point at the end. Distortion vertices are preserved.
    pub fn cells_to_boundary_rings(
        cells: impl IntoIterator<Item = Self>,
    ) -> Result<Vec<H3BoundaryRing>, H3PolygonError> {
        let cells = cells.into_iter().collect::<Vec<_>>();
        let Some(first) = cells.first().copied() else {
            return Ok(Vec::new());
        };
        if cells
            .iter()
            .any(|cell| cell.resolution() != first.resolution())
        {
            return Err(H3PolygonError::HeterogeneousResolution);
        }

        let mut set = BTreeSet::new();
        for cell in cells {
            if !set.insert(cell) {
                return Err(H3PolygonError::DuplicateInput);
            }
        }

        let mut arcs = Vec::new();
        for cell in &set {
            for edge in cell.edges() {
                if !set.contains(&edge.destination()) {
                    arcs.push(Arc::from_edge(edge));
                }
            }
        }
        rings_from_arcs(arcs)
    }
}

/// A dissolved H3 cell-set boundary ring.
#[derive(Clone, Debug, PartialEq)]
pub struct H3BoundaryRing {
    points: Vec<(f64, f64)>,
}

impl H3BoundaryRing {
    /// Ring points as `(latitude, longitude)` degrees, without a duplicated
    /// closing point.
    pub fn points(&self) -> &[(f64, f64)] {
        &self.points
    }

    /// Area enclosed by the ring in radians squared.
    pub fn area_rads2(&self) -> f64 {
        linear_ring_area_rads2(&self.points)
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

#[derive(Clone, Copy, Debug, PartialEq)]
struct BBox {
    min_lat: f64,
    max_lat: f64,
    min_lng: f64,
    max_lng: f64,
}

impl BBox {
    fn from_rings<'a>(
        rings: impl IntoIterator<Item = &'a Vec<(f64, f64)>>,
    ) -> Result<Self, H3PolygonError> {
        let mut points = Vec::new();
        for ring in rings {
            points.extend(ring.iter().copied());
        }
        Self::from_points(points).ok_or(H3PolygonError::RingTooShort)
    }

    fn from_points(points: impl IntoIterator<Item = (f64, f64)>) -> Option<Self> {
        let mut bbox = None::<Self>;
        for (lat, lng) in points {
            bbox = Some(match bbox {
                Some(current) => Self {
                    min_lat: current.min_lat.min(lat),
                    max_lat: current.max_lat.max(lat),
                    min_lng: current.min_lng.min(lng),
                    max_lng: current.max_lng.max(lng),
                },
                None => Self {
                    min_lat: lat,
                    max_lat: lat,
                    min_lng: lng,
                    max_lng: lng,
                },
            });
        }
        bbox
    }

    fn contains(self, point: (f64, f64)) -> bool {
        point.0 >= self.min_lat
            && point.0 <= self.max_lat
            && point.1 >= self.min_lng
            && point.1 <= self.max_lng
    }

    fn intersects(self, other: &Self) -> bool {
        self.min_lat <= other.max_lat
            && self.max_lat >= other.min_lat
            && self.min_lng <= other.max_lng
            && self.max_lng >= other.min_lng
    }

    fn center(self) -> (f64, f64) {
        (
            (self.min_lat + self.max_lat) * 0.5,
            (self.min_lng + self.max_lng) * 0.5,
        )
    }
}

#[derive(Clone, Debug)]
struct Arc {
    edge: H3DirectedEdge,
    start: H3Vertex,
    end: H3Vertex,
    points: Vec<(f64, f64)>,
}

impl Arc {
    fn from_edge(edge: H3DirectedEdge) -> Self {
        let (start, end) = edge.vertices();
        Self {
            edge,
            start,
            end,
            points: edge.boundary_lat_lng(),
        }
    }
}

fn rings_from_arcs(arcs: Vec<Arc>) -> Result<Vec<H3BoundaryRing>, H3PolygonError> {
    let mut by_start: BTreeMap<H3Vertex, Vec<usize>> = BTreeMap::new();
    for (idx, arc) in arcs.iter().enumerate() {
        by_start.entry(arc.start).or_default().push(idx);
    }

    let mut visited = vec![false; arcs.len()];
    let mut rings = Vec::new();
    for idx in 0..arcs.len() {
        if visited[idx] {
            continue;
        }

        let start = arcs[idx].start;
        let mut current = idx;
        let mut points = Vec::new();

        loop {
            if visited[current] {
                return Err(H3PolygonError::OpenBoundary);
            }
            visited[current] = true;

            if points.is_empty() {
                points.extend_from_slice(&arcs[current].points);
            } else {
                points.extend(arcs[current].points.iter().copied().skip(1));
            }

            let next_vertex = arcs[current].end;
            if next_vertex == start {
                break;
            }

            let Some(candidates) = by_start.get(&next_vertex) else {
                return Err(H3PolygonError::OpenBoundary);
            };
            current = candidates
                .iter()
                .copied()
                .filter(|candidate| !visited[*candidate])
                .min_by_key(|candidate| arcs[*candidate].edge)
                .ok_or(H3PolygonError::OpenBoundary)?;
        }

        if points.len() > 1 && points.first() == points.last() {
            points.pop();
        }
        rings.push(H3BoundaryRing { points });
    }

    rings.sort_by(|a, b| {
        b.area_rads2()
            .partial_cmp(&a.area_rads2())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(rings)
}

fn normalize_ring(mut ring: Vec<(f64, f64)>) -> Result<Vec<(f64, f64)>, H3PolygonError> {
    if ring.first() == ring.last() {
        ring.pop();
    }
    if ring.len() < 3 {
        return Err(H3PolygonError::RingTooShort);
    }
    for &(lat, lng) in &ring {
        if !lat.is_finite()
            || !lng.is_finite()
            || !(-90.0..=90.0).contains(&lat)
            || !(-180.0..=180.0).contains(&lng)
        {
            return Err(H3PolygonError::InvalidCoordinate);
        }
    }
    Ok(ring)
}

/// True iff some ring edge spans more than 180° of longitude, meaning the ring
/// takes the short way across the antimeridian (mirrors the H3/h3o heuristic).
fn ring_is_transmeridian(ring: &[(f64, f64)]) -> bool {
    ring_edges(ring).any(|(a, b)| (a.1 - b.1).abs() > 180.0)
}

/// Shifts west-hemisphere longitudes by `+360` so an antimeridian-crossing
/// ring becomes continuous.
fn shift_ring_east(ring: &mut [(f64, f64)]) {
    for point in ring {
        if point.1 < 0.0 {
            point.1 += 360.0;
        }
    }
}

/// Ring edges as consecutive point pairs, closing back to the first point.
fn ring_edges(ring: &[(f64, f64)]) -> impl Iterator<Item = ((f64, f64), (f64, f64))> {
    ring.iter()
        .copied()
        .zip(ring.iter().copied().cycle().skip(1))
        .take(ring.len())
}

/// Twice the signed area of triangle `abc` in the lat/lng plane.
fn orientation(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0).mul_add(c.1 - a.1, -((b.1 - a.1) * (c.0 - a.0)))
}

/// True iff segments `ab` and `cd` cross at an interior point of both —
/// endpoint touches and collinear overlaps do not count, so a shared boundary
/// does not break full containment.
fn segments_cross_properly(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64)) -> bool {
    let d1 = orientation(c, d, a);
    let d2 = orientation(c, d, b);
    let d3 = orientation(a, b, c);
    let d4 = orientation(a, b, d);
    ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
}

/// True iff segments `ab` and `cd` intersect at all, including endpoint
/// touches and collinear overlaps.
fn segments_intersect(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64)) -> bool {
    if segments_cross_properly(a, b, c, d) {
        return true;
    }
    let d1 = orientation(c, d, a);
    let d2 = orientation(c, d, b);
    let d3 = orientation(a, b, c);
    let d4 = orientation(a, b, d);
    (d1 == 0.0 && collinear_point_on_segment(c, d, a))
        || (d2 == 0.0 && collinear_point_on_segment(c, d, b))
        || (d3 == 0.0 && collinear_point_on_segment(a, b, c))
        || (d4 == 0.0 && collinear_point_on_segment(a, b, d))
}

/// Bounds check for a point already known to be collinear with segment `ab`.
fn collinear_point_on_segment(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> bool {
    p.0 >= a.0.min(b.0) && p.0 <= a.0.max(b.0) && p.1 >= a.1.min(b.1) && p.1 <= a.1.max(b.1)
}

fn ring_contains_point(ring: &[(f64, f64)], point: (f64, f64)) -> bool {
    let (lat, lng) = point;
    let mut inside = false;

    for (&a, &b) in ring
        .iter()
        .zip(ring.iter().cycle().skip(1))
        .take(ring.len())
    {
        if point_on_segment(a, b, point) {
            return true;
        }

        let (lat_a, lng_a) = a;
        let (lat_b, lng_b) = b;
        let crosses = (lat_a > lat) != (lat_b > lat);
        if crosses {
            let lng_at_lat = (lng_b - lng_a).mul_add((lat - lat_a) / (lat_b - lat_a), lng_a);
            if lng < lng_at_lat {
                inside = !inside;
            }
        }
    }

    inside
}

fn point_on_segment(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> bool {
    const EPS: f64 = 1e-12;
    let cross = (p.0 - a.0).mul_add(b.1 - a.1, -((p.1 - a.1) * (b.0 - a.0)));
    if cross.abs() > EPS {
        return false;
    }
    let min_lat = a.0.min(b.0) - EPS;
    let max_lat = a.0.max(b.0) + EPS;
    let min_lng = a.1.min(b.1) - EPS;
    let max_lng = a.1.max(b.1) + EPS;
    p.0 >= min_lat && p.0 <= max_lat && p.1 >= min_lng && p.1 <= max_lng
}

fn linear_ring_area_rads2(ring: &[(f64, f64)]) -> f64 {
    let mut adder = FloatAdder::default();
    for (a, b) in ring
        .iter()
        .zip(ring.iter().cycle().skip(1))
        .take(ring.len())
    {
        adder.add(cagnoli(
            (a.1.to_radians(), a.0.to_radians()),
            (b.1.to_radians(), b.0.to_radians()),
        ));
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
