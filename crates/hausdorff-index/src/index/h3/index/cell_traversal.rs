//! [`H3Cell`] grid traversal: disk/ring enumeration (`grid_disk`,
//! `grid_disk_distances`, `grid_ring`), grid distance and adjacency
//! (`grid_distance`, `is_neighbor_with`), local-IJ conversion (`to_local_ij`,
//! `from_local_ij`), and shortest grid paths (`grid_path_cells*`).

use super::*;

impl H3Cell {
    /// All cells within grid distance `k` (the filled disk, including `self`).
    pub fn grid_disk(self, k: u32) -> Vec<Self> {
        grid::grid_disk(self.0, k).into_iter().map(Self).collect()
    }

    /// All cells within grid distance `k`, paired with distance from `self`.
    pub fn grid_disk_distances(self, k: u32) -> Vec<(Self, u32)> {
        grid::grid_disk_distances(self.0, k)
            .into_iter()
            .map(|(cell, distance)| (Self(cell), distance))
            .collect()
    }

    /// Cells at exactly grid distance `k` from `self`.
    pub fn grid_ring(self, k: u32) -> Vec<Self> {
        self.grid_disk_distances(k)
            .into_iter()
            .filter_map(|(cell, distance)| (distance == k).then_some(cell))
            .collect()
    }

    /// Grid distance to another cell, or `None` if undefined (different
    /// resolution, too far, or across a pentagon).
    pub fn grid_distance(self, other: Self) -> Option<i32> {
        grid::grid_distance(self.0, other.0)
    }

    /// This cell's IJ coordinate in `origin`'s local coordinate system.
    pub fn to_local_ij(self, origin: Self) -> Option<H3LocalIJ> {
        let ij = coord::CoordIJ::from(grid::local_ijk(self.0, origin.0)?);
        Some(H3LocalIJ {
            anchor: origin,
            i: ij.i,
            j: ij.j,
        })
    }

    /// The cell at `ij` in its anchor's local coordinate system — the inverse
    /// of [`to_local_ij`](Self::to_local_ij). Returns `None` when the
    /// coordinate is outside the anchor's local frame (too far away or across
    /// a pentagon's deleted subsequence).
    pub fn from_local_ij(ij: H3LocalIJ) -> Option<Self> {
        ij.to_cell()
    }

    /// True iff `other` is one grid step away from this cell.
    pub fn is_neighbor_with(self, other: Self) -> bool {
        self.grid_distance(other) == Some(1)
    }

    /// Number of cells in a shortest grid path from this cell to `other`,
    /// including both endpoints.
    pub fn grid_path_cells_size(self, other: Self) -> Option<u32> {
        self.grid_distance(other)
            .and_then(|distance| u32::try_from(distance + 1).ok())
    }

    /// A shortest grid path from this cell to `other`, including both
    /// endpoints. Returns `None` when H3 local-IJ distance is undefined.
    pub fn grid_path_cells(self, other: Self) -> Option<Vec<Self>> {
        let mut remaining = self.grid_distance(other)?;
        let mut current = self;
        let mut path = Vec::with_capacity(usize::try_from(remaining + 1).ok()?);
        path.push(current);

        while remaining > 0 {
            let mut neighbors = current.grid_ring(1);
            neighbors.sort_unstable();
            current = neighbors
                .into_iter()
                .find(|neighbor| neighbor.grid_distance(other) == Some(remaining - 1))?;
            path.push(current);
            remaining -= 1;
        }

        Some(path)
    }
}
