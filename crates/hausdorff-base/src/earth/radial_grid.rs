//! Radial sampling: concentric rings of equally-spaced points around a center.
//!
//! A [`RadialGrid`] describes `ring_count` concentric rings spaced `step` apart,
//! each carrying `points_per_ring` points spread evenly around the ring. It is a
//! cheap `Copy` descriptor — points are computed on demand, so constructing a
//! grid is free and you only pay for the points you query.
//!
//! ```
//! use euclid::LengthMeters;
//! use hausdorff_base::earth::{GeoPoint, RadialGrid};
//!
//! // 15 rings, +100 m radius per ring, 32 points on each ring.
//! let grid = RadialGrid::new(
//!     GeoPoint::from_lat_lon(52.52, 13.405),
//!     15,
//!     LengthMeters::new(100.0),
//!     32,
//! );
//! assert_eq!(grid.len(), 15 * 32);
//!
//! // Query the whole thing, a single ring, or one point.
//! let all: Vec<_> = grid.points().collect();      // flat, ring-major, with metadata
//! let outer = grid.ring(14);                       // Vec<GeoPoint> for the outermost ring
//! let p = grid.point(0, 0);                         // ring 0, spoke 0 (due North)
//! # let _ = (all, outer, p);
//! ```

use super::{Bearing, GeoPoint, destination_m};
use euclid::LengthMeters;

/// One sampled point on a [`RadialGrid`], carrying its polar coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RingPoint {
    /// Ring index, 0-based; ring 0 is the innermost, at `first_radius` from the center.
    pub ring: u32,
    /// Spoke index, 0-based, increasing clockwise from the bearing offset.
    pub spoke: u32,
    /// Radius of this ring, measured from the center.
    pub radius: LengthMeters,
    /// Compass bearing of this spoke (0 = North, clockwise).
    pub bearing: Bearing,
    /// The sampled geographic point.
    pub point: GeoPoint,
}

/// Concentric rings of equally-spaced sample points around a center point.
///
/// The radius schedule is linear with a configurable starting offset: ring `r`
/// (0-based) sits at radius `first_radius + r * step`. With the default
/// [`new`](Self::new) the first ring is at `step` (i.e. `(r + 1) * step`); use
/// [`linear`](Self::linear) or [`with_first_radius`](Self::with_first_radius) to
/// set a different inner offset (e.g. start at 50 m, then +100 m per ring).
///
/// Spoke `s` (0-based) sits at bearing `bearing_offset + s * 360 /
/// points_per_ring` degrees, measured clockwise from North. Iteration is
/// **ring-major**: every spoke of ring 0, then ring 1, and so on.
///
/// `points_per_ring` may be any positive count — powers/multiples of two are the
/// common case, but odd counts work identically. A grid with `ring_count` or
/// `points_per_ring` of zero is simply empty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialGrid {
    center: GeoPoint,
    ring_count: u32,
    first_radius_m: f64,
    step_m: f64,
    points_per_ring: u32,
    bearing_offset_deg: f64,
}

impl RadialGrid {
    /// Build a grid of `ring_count` rings spaced `step` apart, each with
    /// `points_per_ring` equally-spaced points, the first spoke pointing North.
    /// The innermost ring sits at `step` (i.e. ring `r` at `(r + 1) * step`).
    pub fn new(
        center: GeoPoint,
        ring_count: u32,
        step: LengthMeters,
        points_per_ring: u32,
    ) -> Self {
        Self::linear(center, ring_count, step, step, points_per_ring)
    }

    /// Build a grid whose innermost ring sits at `first_radius` and whose
    /// successive rings are `step` apart: ring `r` at `first_radius + r * step`.
    /// This is the "initial constant + linear" schedule (e.g. `first_radius`
    /// 50 m, `step` 100 m → rings at 50, 150, 250, … m).
    pub fn linear(
        center: GeoPoint,
        ring_count: u32,
        first_radius: LengthMeters,
        step: LengthMeters,
        points_per_ring: u32,
    ) -> Self {
        Self {
            center,
            ring_count,
            first_radius_m: first_radius.as_meters(),
            step_m: step.as_meters(),
            points_per_ring,
            bearing_offset_deg: 0.0,
        }
    }

    /// Override the innermost ring's radius, keeping the same `step` between
    /// rings (ring `r` at `first_radius + r * step`).
    pub fn with_first_radius(mut self, first_radius: LengthMeters) -> Self {
        self.first_radius_m = first_radius.as_meters();
        self
    }

    /// Rotate the first spoke `degrees` clockwise off North (the others follow).
    pub fn with_bearing_offset(mut self, degrees: f64) -> Self {
        self.bearing_offset_deg = degrees;
        self
    }

    /// The center the rings are anchored at.
    pub fn center(&self) -> GeoPoint {
        self.center
    }

    /// Number of rings.
    pub fn ring_count(&self) -> u32 {
        self.ring_count
    }

    /// Number of equally-spaced points on each ring.
    pub fn points_per_ring(&self) -> u32 {
        self.points_per_ring
    }

    /// Radius of the innermost ring (ring 0).
    pub fn first_radius(&self) -> LengthMeters {
        LengthMeters::new(self.first_radius_m)
    }

    /// Radial spacing between successive rings.
    pub fn step(&self) -> LengthMeters {
        LengthMeters::new(self.step_m)
    }

    /// Total number of sample points (`ring_count * points_per_ring`).
    pub fn len(&self) -> usize {
        self.ring_count as usize * self.points_per_ring as usize
    }

    /// Whether the grid has no sample points.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Radius of ring `r` (0-based): `first_radius + r * step`.
    pub fn radius(&self, ring: u32) -> LengthMeters {
        LengthMeters::new(self.first_radius_m + f64::from(ring) * self.step_m)
    }

    /// Compass bearing of spoke `s`, clockwise from North. A grid with zero
    /// `points_per_ring` is empty and has no spokes, so the bearing is just the
    /// offset (avoiding a divide-by-zero).
    pub fn bearing(&self, spoke: u32) -> Bearing {
        if self.points_per_ring == 0 {
            return Bearing::from_degrees(self.bearing_offset_deg);
        }
        Bearing::from_degrees(
            self.bearing_offset_deg + f64::from(spoke) * 360.0 / f64::from(self.points_per_ring),
        )
    }

    /// The geographic point at `(ring, spoke)`. Indices beyond the configured
    /// counts are not rejected — the polar formula is simply evaluated at the
    /// implied radius/bearing.
    pub fn point(&self, ring: u32, spoke: u32) -> GeoPoint {
        destination_m(
            self.center,
            self.bearing(spoke),
            self.radius(ring).as_meters(),
        )
    }

    /// The full [`RingPoint`] (point plus polar metadata) at `(ring, spoke)`.
    pub fn ring_point(&self, ring: u32, spoke: u32) -> RingPoint {
        RingPoint {
            ring,
            spoke,
            radius: self.radius(ring),
            bearing: self.bearing(spoke),
            point: self.point(ring, spoke),
        }
    }

    /// The points of a single ring, ordered by spoke.
    pub fn ring(&self, ring: u32) -> Vec<GeoPoint> {
        (0..self.points_per_ring)
            .map(|spoke| self.point(ring, spoke))
            .collect()
    }

    /// Iterator over each ring's points, from innermost to outermost.
    pub fn rings(&self) -> impl Iterator<Item = Vec<GeoPoint>> + '_ {
        (0..self.ring_count).map(move |ring| self.ring(ring))
    }

    /// Flat, ring-major iterator over every sample point with its metadata.
    pub fn points(&self) -> impl Iterator<Item = RingPoint> + '_ {
        (0..self.ring_count).flat_map(move |ring| {
            (0..self.points_per_ring).map(move |spoke| self.ring_point(ring, spoke))
        })
    }

    /// Flat, ring-major collection of just the geographic points — handy for
    /// snapping a whole grid to a graph at once.
    pub fn to_points(&self) -> Vec<GeoPoint> {
        self.points().map(|rp| rp.point).collect()
    }
}
