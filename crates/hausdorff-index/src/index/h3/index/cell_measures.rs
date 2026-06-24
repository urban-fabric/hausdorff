//! [`H3Cell`] geometric and tabular measures: the cell boundary
//! (`boundary_lat_lng`), exact spherical area (`area*`), great-circle distance
//! between centers (`great_circle_distance*`), and the resolution-indexed H3
//! reference tables (`max_grid_disk_size`, `average_hexagon_*`, `num_cells`).

use hausdorff_base::earth::{AreaSquareMeters, LengthMeters};

use super::*;

impl H3Cell {
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

    /// Great-circle (haversine) distance between this cell's center and
    /// `other`'s, in radians.
    pub fn great_circle_distance_rads(self, other: Self) -> f64 {
        distance_rads(self.center_lat_lng_rad(), other.center_lat_lng_rad())
    }

    /// Great-circle (haversine) distance between this cell's center and
    /// `other`'s, in kilometers.
    pub fn great_circle_distance_km(self, other: Self) -> f64 {
        self.great_circle_distance_rads(other) * EARTH_RADIUS_KM
    }

    /// Great-circle (haversine) distance between this cell's center and
    /// `other`'s, in meters.
    pub fn great_circle_distance_m(self, other: Self) -> f64 {
        self.great_circle_distance_km(other) * 1000.0
    }

    /// Great-circle (haversine) distance between this cell's center and
    /// `other`'s as a typed [`LengthMeters`] — prefer this over
    /// [`great_circle_distance_m`](Self::great_circle_distance_m) so the
    /// Earth-distance unit travels with the value (and isn't mistaken for the
    /// radian measure [`great_circle_distance_rads`](Self::great_circle_distance_rads)).
    pub fn great_circle_distance(self, other: Self) -> LengthMeters {
        LengthMeters::new(self.great_circle_distance_m(other))
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
}
