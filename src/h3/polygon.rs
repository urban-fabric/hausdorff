//! Lightweight polygon operations for H3 cells.

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;

use super::index::{H3Cell, H3DirectedEdge, H3Vertex, MAX_RESOLUTION};

/// A geographic polygon in `(latitude, longitude)` degrees.
///
/// Rings are stored without a duplicated closing point. The current cover
/// implementation supports ordinary non-transmeridian polygons and uses H3
/// centroid containment semantics.
#[derive(Clone, Debug, PartialEq)]
pub struct H3Polygon {
    exterior: Vec<(f64, f64)>,
    holes: Vec<Vec<(f64, f64)>>,
    bbox: BBox,
}

impl H3Polygon {
    /// Constructs a polygon from an exterior ring and zero or more hole rings.
    pub fn new(
        exterior: impl Into<Vec<(f64, f64)>>,
        holes: impl Into<Vec<Vec<(f64, f64)>>>,
    ) -> Result<Self, H3PolygonError> {
        let exterior = normalize_ring(exterior.into())?;
        let holes = holes
            .into()
            .into_iter()
            .map(normalize_ring)
            .collect::<Result<Vec<_>, _>>()?;
        let bbox = BBox::from_rings(std::iter::once(&exterior).chain(holes.iter()))?;
        if bbox.transmeridian() {
            return Err(H3PolygonError::TransmeridianUnsupported);
        }

        Ok(Self {
            exterior,
            holes,
            bbox,
        })
    }

    /// Exterior ring, without a duplicated closing point.
    pub fn exterior(&self) -> &[(f64, f64)] {
        &self.exterior
    }

    /// Hole rings, without duplicated closing points.
    pub fn holes(&self) -> &[Vec<(f64, f64)>] {
        &self.holes
    }

    /// True iff the point `(latitude, longitude)` in degrees is inside the
    /// polygon and outside all holes.
    pub fn contains_lat_lng(&self, lat: f64, lng: f64) -> bool {
        let point = (lat, lng);
        self.bbox.contains(point)
            && ring_contains_point(&self.exterior, point)
            && !self
                .holes
                .iter()
                .any(|hole| ring_contains_point(hole, point))
    }

    /// H3 cells at `resolution` whose centers are inside this polygon.
    pub fn to_cells(&self, resolution: u8) -> Result<Vec<H3Cell>, H3PolygonError> {
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

            let Some(cell_bbox) = BBox::from_points(cell.boundary_lat_lng()) else {
                continue;
            };
            if !self.bbox.intersects(&cell_bbox) {
                continue;
            }

            let (lat, lng) = cell.center_lat_lng();
            if self.contains_lat_lng(lat, lng) {
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

    fn seed_cells(&self, resolution: u8) -> Vec<H3Cell> {
        let mut seeds = Vec::new();
        let mut push_seed = |point: (f64, f64)| {
            if let Some(cell) = H3Cell::for_point(point.0, point.1, resolution) {
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
    /// Coordinates must be finite and latitudes must be in `[-90, 90]`.
    InvalidCoordinate,
    /// First pass does not yet split rings crossing the antimeridian.
    TransmeridianUnsupported,
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
            Self::TransmeridianUnsupported => {
                f.write_str("transmeridian polygons are not yet supported")
            }
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

    fn transmeridian(self) -> bool {
        self.max_lng - self.min_lng > 180.0
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
        if !lat.is_finite() || !lng.is_finite() || !(-90.0..=90.0).contains(&lat) {
            return Err(H3PolygonError::InvalidCoordinate);
        }
    }
    Ok(ring)
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
