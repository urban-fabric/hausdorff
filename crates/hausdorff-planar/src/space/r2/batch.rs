//! `space::r2::batch` — Structure-of-Arrays planar point batches.
//!
//! A columnar (`xs[]`, `ys[]`) counterpart to the Array-of-Structs [`Coord`] for
//! allocation-friendly bulk scans (bounding rectangles, and a base for batched
//! planar kernels). The core geometry model is unchanged; this is an optional
//! view callers opt into for large point sets.

use super::{Coord, Rect};

/// A batch of planar points stored column-wise (`xs`, `ys`, equal length).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PointBatch {
    xs: Vec<f64>,
    ys: Vec<f64>,
}

impl PointBatch {
    /// Build from parallel `xs`/`ys` columns, or `None` if their lengths differ.
    #[must_use]
    pub fn from_columns(xs: Vec<f64>, ys: Vec<f64>) -> Option<Self> {
        (xs.len() == ys.len()).then_some(Self { xs, ys })
    }

    /// Build the columnar batch from an Array-of-Structs slice of coordinates.
    #[must_use]
    pub fn from_coords(coords: &[Coord]) -> Self {
        Self {
            xs: coords.iter().map(|c| c.x).collect(),
            ys: coords.iter().map(|c| c.y).collect(),
        }
    }

    /// Number of points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.xs.len()
    }

    /// Whether the batch has no points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.xs.is_empty()
    }

    /// The x column.
    #[must_use]
    pub fn xs(&self) -> &[f64] {
        &self.xs
    }

    /// The y column.
    #[must_use]
    pub fn ys(&self) -> &[f64] {
        &self.ys
    }

    /// The `i`-th point as a [`Coord`].
    #[must_use]
    pub fn coord(&self, i: usize) -> Coord {
        Coord::new(self.xs[i], self.ys[i])
    }

    /// Planar (Euclidean) distance from every point to a fixed `target` — the
    /// one-to-many distance fan over the columnar layout. Runs a runtime-selected
    /// SIMD kernel (AVX on supporting x86-64 CPUs, scalar otherwise); the result
    /// is identical either way.
    #[must_use]
    pub fn distances_to(&self, target: Coord) -> Vec<f64> {
        let mut out = vec![0.0_f64; self.len()];
        crate::simd::euclidean_distances(&self.xs, &self.ys, target.x, target.y, &mut out);
        out
    }

    /// The axis-aligned bounding rectangle of the batch, or `None` if empty.
    #[must_use]
    pub fn bounding_rect(&self) -> Option<Rect> {
        if self.is_empty() {
            return None;
        }
        let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
        for (&x, &y) in self.xs.iter().zip(&self.ys) {
            min_x = min_x.min(x);
            max_x = max_x.max(x);
            min_y = min_y.min(y);
            max_y = max_y.max(y);
        }
        Some(Rect::new(
            Coord::new(min_x, min_y),
            Coord::new(max_x, max_y),
        ))
    }
}
