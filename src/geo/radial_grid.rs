//! Radial sampling: concentric rings of equally-spaced points around a center.
//!
//! A [`RadialGrid`] describes `ring_count` concentric rings spaced `step` apart,
//! each carrying `points_per_ring` points spread evenly around the ring. It is a
//! cheap `Copy` descriptor — points are computed on demand, so constructing a
//! grid is free and you only pay for the points you query.
//!
//! ```
//! use euclid::LengthMeters;
//! use hausdorff::geo::{GeoPoint, RadialGrid};
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
    /// Ring index, 0-based; ring 0 is the innermost, one `step` from the center.
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

    /// Compass bearing of spoke `s`, clockwise from North.
    pub fn bearing(&self, spoke: u32) -> Bearing {
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

#[cfg(test)]
mod tests {
    use super::super::{GeoPoint, haversine_m};
    use super::{RadialGrid, RingPoint};
    use euclid::LengthMeters;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn counts_and_dimensions() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(52.52, 13.405),
            15,
            LengthMeters::new(100.0),
            32,
        );
        assert_eq!(g.ring_count(), 15);
        assert_eq!(g.points_per_ring(), 32);
        assert_eq!(g.len(), 480);
        assert!(!g.is_empty());
        assert_eq!(g.to_points().len(), 480);
        assert_eq!(g.points().count(), 480);
        assert_eq!(g.rings().count(), 15);
    }

    #[test]
    fn radius_and_bearing_schedule() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            3,
            LengthMeters::new(100.0),
            4,
        );
        assert!(approx(g.radius(0).as_meters(), 100.0, 1e-9));
        assert!(approx(g.radius(2).as_meters(), 300.0, 1e-9));
        assert!(approx(g.bearing(0).degrees(), 0.0, 1e-9));
        assert!(approx(g.bearing(1).degrees(), 90.0, 1e-9));
        assert!(approx(g.bearing(2).degrees(), 180.0, 1e-9));
        assert!(approx(g.bearing(3).degrees(), 270.0, 1e-9));
    }

    #[test]
    fn ring_major_ordering() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            2,
            LengthMeters::new(100.0),
            4,
        );
        let pts: Vec<RingPoint> = g.points().collect();
        assert_eq!((pts[0].ring, pts[0].spoke), (0, 0));
        assert_eq!((pts[3].ring, pts[3].spoke), (0, 3));
        assert_eq!((pts[4].ring, pts[4].spoke), (1, 0));
        // `ring(i)` matches the corresponding flat slice.
        let ring1 = g.ring(1);
        assert_eq!(ring1.len(), 4);
        assert_eq!(ring1[0], pts[4].point);
    }

    #[test]
    fn cardinal_directions_point_the_right_way() {
        let center = GeoPoint::from_lat_lon(0.0, 0.0);
        let g = RadialGrid::new(center, 1, LengthMeters::new(1_000.0), 4);
        let north = g.point(0, 0);
        let east = g.point(0, 1);
        let south = g.point(0, 2);
        let west = g.point(0, 3);
        assert!(north.lat > center.lat && approx(north.lon, 0.0, 1e-6));
        assert!(east.lon > center.lon && approx(east.lat, 0.0, 1e-6));
        assert!(south.lat < center.lat && approx(south.lon, 0.0, 1e-6));
        assert!(west.lon < center.lon && approx(west.lat, 0.0, 1e-6));
    }

    #[test]
    fn points_sit_at_the_ring_radius() {
        let center = GeoPoint::from_lat_lon(52.52, 13.405);
        let g = RadialGrid::new(center, 15, LengthMeters::new(100.0), 32);
        for rp in g.points() {
            let d = haversine_m(center, rp.point);
            assert!(
                approx(d, rp.radius.as_meters(), 0.5),
                "ring {} spoke {}: {d} m vs nominal {} m",
                rp.ring,
                rp.spoke,
                rp.radius.as_meters(),
            );
        }
    }

    #[test]
    fn linear_with_initial_constant_offsets_the_inner_ring() {
        // first ring at 50 m, then +100 m per ring -> 50, 150, 250 m.
        let g = RadialGrid::linear(
            GeoPoint::from_lat_lon(0.0, 0.0),
            3,
            LengthMeters::new(50.0),
            LengthMeters::new(100.0),
            8,
        );
        assert!(approx(g.first_radius().as_meters(), 50.0, 1e-9));
        assert!(approx(g.radius(0).as_meters(), 50.0, 1e-9));
        assert!(approx(g.radius(1).as_meters(), 150.0, 1e-9));
        assert!(approx(g.radius(2).as_meters(), 250.0, 1e-9));
    }

    #[test]
    fn with_first_radius_matches_linear_constructor() {
        let center = GeoPoint::from_lat_lon(52.52, 13.405);
        let a = RadialGrid::new(center, 4, LengthMeters::new(100.0), 8)
            .with_first_radius(LengthMeters::new(50.0));
        let b = RadialGrid::linear(
            center,
            4,
            LengthMeters::new(50.0),
            LengthMeters::new(100.0),
            8,
        );
        assert_eq!(a, b);
    }

    #[test]
    fn new_keeps_first_ring_at_step() {
        // The simple `new` schedule: ring r at (r + 1) * step.
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            3,
            LengthMeters::new(100.0),
            8,
        );
        assert!(approx(g.radius(0).as_meters(), 100.0, 1e-9));
        assert!(approx(g.radius(2).as_meters(), 300.0, 1e-9));
    }

    #[test]
    fn bearing_offset_rotates_the_first_spoke() {
        let center = GeoPoint::from_lat_lon(0.0, 0.0);
        let g = RadialGrid::new(center, 1, LengthMeters::new(1_000.0), 4).with_bearing_offset(90.0);
        assert!(approx(g.bearing(0).degrees(), 90.0, 1e-9));
        let p = g.point(0, 0);
        assert!(p.lon > center.lon && approx(p.lat, 0.0, 1e-6)); // first spoke now points East
    }

    #[test]
    fn odd_points_per_ring_is_supported() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            2,
            LengthMeters::new(100.0),
            5,
        );
        assert_eq!(g.len(), 10);
        assert!(approx(g.bearing(1).degrees(), 72.0, 1e-9));
    }

    #[test]
    fn empty_grid_yields_no_points() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            0,
            LengthMeters::new(100.0),
            8,
        );
        assert!(g.is_empty());
        assert_eq!(g.points().count(), 0);
        assert_eq!(g.to_points().len(), 0);
    }
}
