use std::collections::{BTreeMap, BTreeSet};

use super::super::index::{H3Cell, H3DirectedEdge, H3Vertex};
use super::H3PolygonError;
use super::ring::linear_ring_area_rads2;

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

impl H3Cell {
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
