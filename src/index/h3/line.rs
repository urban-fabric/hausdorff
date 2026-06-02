//! Lightweight polyline operations for H3 cells — the line counterpart to
//! [`super::polygon`].

use std::collections::BTreeSet;
use std::fmt;

use super::index::{H3Cell, MAX_RESOLUTION};

/// A geographic polyline in `(latitude, longitude)` degrees, sharing the
/// coordinate convention of [`super::polygon::H3Polygon`].
///
/// [`to_cells`](H3Line::to_cells) covers the line by snapping each vertex to its
/// cell and walking the H3 grid line between consecutive vertex cells (the
/// `gridPathCells` idiom). The result is a *connected* chain of cells from the
/// first vertex to the last — a guaranteed grid-path traversal rather than an
/// unguaranteed point sample, which is why a polyline cover cannot be faked by
/// snapping a handful of sampled points.
#[derive(Clone, Debug, PartialEq)]
pub struct H3Line {
    points: Vec<(f64, f64)>,
}

impl H3Line {
    /// Constructs a polyline from its `(latitude, longitude)` degree vertices.
    ///
    /// Requires at least one vertex; every coordinate must be finite, with
    /// `lat ∈ [-90, 90]` and `lng ∈ [-180, 180]`.
    pub fn new(points: impl Into<Vec<(f64, f64)>>) -> Result<Self, H3LineError> {
        let points = points.into();
        if points.is_empty() {
            return Err(H3LineError::Empty);
        }
        for &(lat, lng) in &points {
            if !lat.is_finite()
                || !lng.is_finite()
                || !(-90.0..=90.0).contains(&lat)
                || !(-180.0..=180.0).contains(&lng)
            {
                return Err(H3LineError::InvalidCoordinate);
            }
        }
        Ok(Self { points })
    }

    /// The polyline vertices as `(latitude, longitude)` degrees.
    pub fn points(&self) -> &[(f64, f64)] {
        &self.points
    }

    /// H3 cells at `resolution` covering the line, ordered from the first vertex
    /// to the last with duplicates removed.
    ///
    /// Each vertex is snapped to its containing cell, and consecutive vertex
    /// cells are joined by the H3 grid line ([`H3Cell::grid_path_cells`]); for a
    /// vertex pair landing in the same cell that step is empty. A single-vertex
    /// line yields just its containing cell. Returns
    /// [`PathUndefined`](H3LineError::PathUndefined) when the grid path between
    /// two vertex cells is undefined (too far apart, or separated by a pentagon),
    /// mirroring H3's `gridPathCells`.
    pub fn to_cells(&self, resolution: u8) -> Result<Vec<H3Cell>, H3LineError> {
        if resolution > MAX_RESOLUTION {
            return Err(H3LineError::InvalidResolution(resolution));
        }

        let mut seen = BTreeSet::new();
        let mut out = Vec::new();
        let mut prev: Option<H3Cell> = None;
        for &(lat, lng) in &self.points {
            let cell =
                H3Cell::for_point(lat, lng, resolution).ok_or(H3LineError::UnmappablePoint)?;
            match prev {
                // Walk the grid line to this vertex's cell; the path includes
                // `previous` (already emitted) and `cell`, and `seen` drops the
                // shared joints and any cell revisited by a self-crossing line.
                Some(previous) if previous != cell => {
                    let path = previous
                        .grid_path_cells(cell)
                        .ok_or(H3LineError::PathUndefined)?;
                    for step in path {
                        if seen.insert(step) {
                            out.push(step);
                        }
                    }
                }
                _ => {
                    if seen.insert(cell) {
                        out.push(cell);
                    }
                }
            }
            prev = Some(cell);
        }

        Ok(out)
    }
}

impl H3Cell {
    /// Convenience wrapper for [`H3Line::to_cells`].
    pub fn line_to_cells(line: &H3Line, resolution: u8) -> Result<Vec<Self>, H3LineError> {
        line.to_cells(resolution)
    }
}

/// Polyline-operation error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum H3LineError {
    /// A polyline needs at least one vertex.
    Empty,
    /// Coordinates must be finite, with `lat ∈ [-90, 90]` and `lng ∈ [-180, 180]`.
    InvalidCoordinate,
    /// Resolution must be in `0..=15`.
    InvalidResolution(u8),
    /// A vertex could not be mapped to an H3 cell.
    UnmappablePoint,
    /// The grid path between two consecutive vertex cells is undefined (the
    /// cells are too far apart or separated by a pentagon).
    PathUndefined,
}

impl fmt::Display for H3LineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("polyline has no vertices"),
            Self::InvalidCoordinate => f.write_str("invalid polyline coordinate"),
            Self::InvalidResolution(resolution) => {
                write!(f, "invalid H3 resolution {resolution}")
            }
            Self::UnmappablePoint => f.write_str("polyline vertex could not be mapped to a cell"),
            Self::PathUndefined => f.write_str("grid path between polyline vertices is undefined"),
        }
    }
}

impl std::error::Error for H3LineError {}
