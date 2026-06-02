//! Uniform lon/lat grid index over a point set, for radius (`ST_DWithin`) and
//! k-nearest-neighbour joins. Points are bucketed into `cell_deg`-degree cells;
//! a radius probe scans the block of cells the search radius can reach — sized
//! from the radius and `cell_deg`, and widened in longitude toward the poles —
//! then refines by great-circle ([`crate::earth::haversine_m`]) distance. So
//! `cell_deg` is only a bucketing granularity (it trades cells-scanned-per-probe
//! against memory) and need not match the radius. This is a **flat** lon/lat
//! grid: it does not wrap the antimeridian and degrades within one cell of a
//! pole. [`PointGrid::build`] rejects a non-finite or non-positive `cell_deg`.
//!
//! Buckets are a **dense** CSR grid when the occupied cell box is compact, or a
//! **sparse** hash map when the points are scattered, so memory stays ~`O(points)`
//! even for globe-spanning sets at a fine `cell_deg`. Generic over a per-point
//! payload `P` (e.g. a source row id), so the caller maps results back to its rows.

use std::collections::HashMap;

use crate::earth::{GeoPoint, haversine_m};

/// Go dense while `occupied_cells <= max(points * FACTOR, FLOOR)`; past that the
/// grid is too sparse to be worth its footprint, so fall back to the hash map.
const DENSE_CELLS_PER_POINT: i128 = 2;
const DENSE_MIN_CELLS: i128 = 4096;

/// Grid cell `(row, col)` of a point, flooring toward negative infinity so it is
/// correct for any coordinate sign.
#[inline]
fn cell(lat: f64, lon: f64, cell_deg: f64) -> (i64, i64) {
    (
        (lat / cell_deg).floor() as i64,
        (lon / cell_deg).floor() as i64,
    )
}

/// Great-circle distance in metres between two lon/lat points.
#[inline]
fn distance_m(a_lat: f64, a_lon: f64, b_lat: f64, b_lon: f64) -> f64 {
    haversine_m(
        GeoPoint::from_lat_lon(a_lat, a_lon),
        GeoPoint::from_lat_lon(b_lat, b_lon),
    )
}

struct GridPoint<P> {
    lat: f64,
    lon: f64,
    payload: P,
}

/// Bucketing of the indexed points by grid cell: a dense flat grid (direct
/// addressing, no hashing) when the occupied cell box is compact, or a sparse
/// hash fallback otherwise.
enum Buckets {
    Dense(DenseGrid),
    Sparse(HashMap<(i64, i64), Vec<u32>>),
}

/// Dense grid in compressed-sparse-row form: `items[starts[s]..starts[s + 1]]`
/// are the point indices in cell `s`. Two contiguous allocations (no per-cell heap
/// node), so probing a cell is a bounds check plus a slice.
struct DenseGrid {
    min_row: i64,
    min_col: i64,
    rows: i64,
    cols: i64,
    starts: Vec<u32>,
    items: Vec<u32>,
}

impl Buckets {
    fn build<P>(points: &[GridPoint<P>], cell_deg: f64) -> Self {
        let Some(first) = points.first() else {
            return Buckets::Sparse(HashMap::new());
        };
        let (mut min_row, mut min_col) = cell(first.lat, first.lon, cell_deg);
        let (mut max_row, mut max_col) = (min_row, min_col);
        for p in &points[1..] {
            let (r, c) = cell(p.lat, p.lon, cell_deg);
            min_row = min_row.min(r);
            max_row = max_row.max(r);
            min_col = min_col.min(c);
            max_col = max_col.max(c);
        }

        // i128 spans avoid overflow for pathological `cell_deg`/coordinate ranges.
        let row_span = (max_row as i128) - (min_row as i128) + 1;
        let col_span = (max_col as i128) - (min_col as i128) + 1;
        let cell_count = row_span.saturating_mul(col_span);
        let budget = ((points.len() as i128) * DENSE_CELLS_PER_POINT).max(DENSE_MIN_CELLS);

        if cell_count > budget {
            let mut map: HashMap<(i64, i64), Vec<u32>> = HashMap::with_capacity(points.len());
            for (point_idx, p) in points.iter().enumerate() {
                map.entry(cell(p.lat, p.lon, cell_deg))
                    .or_default()
                    .push(point_idx as u32);
            }
            return Buckets::Sparse(map);
        }

        // Dense path: counting sort of the points into a CSR grid. `cell_count`
        // is within `budget` here, so the casts to `i64`/`usize` cannot overflow.
        let (rows, cols) = (row_span as i64, col_span as i64);
        let ncells = (rows * cols) as usize;
        let slot = |r: i64, c: i64| ((r - min_row) * cols + (c - min_col)) as usize;

        let mut starts = vec![0u32; ncells + 1];
        for p in points {
            let (r, c) = cell(p.lat, p.lon, cell_deg);
            starts[slot(r, c) + 1] += 1;
        }
        for s in 0..ncells {
            starts[s + 1] += starts[s];
        }
        let mut items = vec![0u32; points.len()];
        let mut cursor = starts.clone();
        for (point_idx, p) in points.iter().enumerate() {
            let (r, c) = cell(p.lat, p.lon, cell_deg);
            let s = slot(r, c);
            items[cursor[s] as usize] = point_idx as u32;
            cursor[s] += 1;
        }

        Buckets::Dense(DenseGrid {
            min_row,
            min_col,
            rows,
            cols,
            starts,
            items,
        })
    }

    /// The point indices bucketed into cell `(row, col)`, or an empty slice when
    /// the cell is empty or outside the indexed extent.
    #[inline]
    fn get(&self, row: i64, col: i64) -> &[u32] {
        match self {
            Buckets::Dense(g) => {
                let (dr, dc) = (row - g.min_row, col - g.min_col);
                if dr < 0 || dc < 0 || dr >= g.rows || dc >= g.cols {
                    return &[];
                }
                let s = (dr * g.cols + dc) as usize;
                &g.items[g.starts[s] as usize..g.starts[s + 1] as usize]
            }
            Buckets::Sparse(map) => map.get(&(row, col)).map_or(&[], Vec::as_slice),
        }
    }
}

/// A reusable uniform-grid index over lon/lat points tagged with a payload —
/// build once, probe many (the broadcast analogue of [`within_join`]).
pub struct PointGrid<P> {
    cell_deg: f64,
    points: Vec<GridPoint<P>>,
    buckets: Buckets,
}

impl<P: Copy> PointGrid<P> {
    /// Index `(latitude, longitude, payload)` points into a `cell_deg`-degree
    /// grid. Returns `None` unless `cell_deg` is finite and positive (a zero,
    /// negative, or `NaN` cell size would produce nonsense buckets).
    #[must_use]
    pub fn build(points: impl IntoIterator<Item = (f64, f64, P)>, cell_deg: f64) -> Option<Self> {
        if !cell_deg.is_finite() || cell_deg <= 0.0 {
            return None;
        }
        let points: Vec<GridPoint<P>> = points
            .into_iter()
            .map(|(lat, lon, payload)| GridPoint { lat, lon, payload })
            .collect();
        let buckets = Buckets::build(&points, cell_deg);
        Some(Self {
            cell_deg,
            points,
            buckets,
        })
    }

    /// Number of indexed points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the index holds no points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// `(probe_payload, point_payload)` for every probe point within `meters`
    /// (great-circle) of an indexed point. Each probe scans the block of cells the
    /// radius can reach — `ceil(radius / cell)` rings in latitude, widened by
    /// `1/cos(lat)` in longitude — then refines by great-circle distance, so the
    /// result is complete for any positive `cell_deg` (a smaller one just scans
    /// more cells). Ring counts are capped at half the grid, so a pole-adjacent
    /// probe stays bounded but may be incomplete there (a flat grid can't
    /// represent a pole neighbourhood).
    #[must_use]
    pub fn within<Q: Copy>(
        &self,
        probe: impl IntoIterator<Item = (f64, f64, Q)>,
        meters: f64,
    ) -> Vec<(Q, P)> {
        const DEG_PER_M: f64 = 1.0 / 111_320.0;
        let radius_deg = meters * DEG_PER_M;
        // cell_deg is finite and > 0 (enforced by `build`), so these are finite.
        let max_rings = ((180.0 / self.cell_deg).ceil() as i64).max(1);
        let lat_rings = ((radius_deg / self.cell_deg).ceil() as i64).clamp(1, max_rings);
        let mut out = Vec::new();
        for (lat, lon, q) in probe {
            // A degree of longitude shrinks as cos(lat), so the radius spans more
            // longitude cells nearer the poles.
            let coslat = lat.to_radians().cos().abs();
            let lon_rings = if coslat < 1e-9 {
                max_rings
            } else {
                ((radius_deg / (coslat * self.cell_deg)).ceil() as i64).clamp(1, max_rings)
            };
            let (cy, cx) = cell(lat, lon, self.cell_deg);
            for dy in -lat_rings..=lat_rings {
                for dx in -lon_rings..=lon_rings {
                    for &i in self.buckets.get(cy + dy, cx + dx) {
                        let p = &self.points[i as usize];
                        if distance_m(lat, lon, p.lat, p.lon) <= meters {
                            out.push((q, p.payload));
                        }
                    }
                }
            }
        }
        out
    }
}

/// k-nearest-neighbour over a candidate point set: for each probe point, the `k`
/// closest of `points` by great-circle distance, as `(probe_payload,
/// point_payload, distance_m)`. A brute-force scan of the candidates (no grid) —
/// the in-memory refine kernel a caller feeds with a pre-filtered neighbourhood;
/// it is exact only when that neighbourhood already holds every true neighbour.
#[must_use]
pub fn k_nearest<P, Q>(
    probe: impl IntoIterator<Item = (f64, f64, Q)>,
    points: &[(f64, f64, P)],
    k: usize,
) -> Vec<(Q, P, f64)>
where
    P: Copy,
    Q: Copy,
{
    let mut out = Vec::new();
    let mut candidates: Vec<(f64, P)> = Vec::with_capacity(points.len());
    for (lat, lon, q) in probe {
        candidates.clear();
        for &(plat, plon, payload) in points {
            candidates.push((distance_m(lat, lon, plat, plon), payload));
        }
        let keep = k.min(candidates.len());
        if keep == 0 {
            continue;
        }
        // Keep the k smallest, then sort just that selection (late materialization:
        // never fully order rows that will not be emitted).
        if keep < candidates.len() {
            candidates.select_nth_unstable_by(keep, |a, b| a.0.total_cmp(&b.0));
        }
        candidates[..keep].sort_by(|a, b| a.0.total_cmp(&b.0));
        for &(d, payload) in candidates.iter().take(keep) {
            out.push((q, payload, d));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    fn pts(lat: &[f64], lon: &[f64]) -> Vec<(f64, f64, u32)> {
        lat.iter()
            .zip(lon)
            .enumerate()
            .map(|(i, (&y, &x))| (y, x, i as u32))
            .collect()
    }

    /// O(n·m) reference: every (probe, point) pair within `meters`.
    fn brute(
        probe: &[(f64, f64, u32)],
        points: &[(f64, f64, u32)],
        meters: f64,
    ) -> BTreeSet<(u32, u32)> {
        let mut pairs = BTreeSet::new();
        for &(ly, lx, li) in probe {
            for &(ry, rx, ri) in points {
                if distance_m(ly, lx, ry, rx) <= meters {
                    pairs.insert((li, ri));
                }
            }
        }
        pairs
    }

    fn within_pairs(
        grid: &PointGrid<u32>,
        probe: &[(f64, f64, u32)],
        meters: f64,
    ) -> BTreeSet<(u32, u32)> {
        grid.within(probe.iter().copied(), meters)
            .into_iter()
            .collect()
    }

    #[test]
    fn within_agrees_with_brute_force_on_both_bucket_paths() {
        let cell_deg = 1.0;
        let meters = 50_000.0;

        // Compact cluster spanning a few cells -> dense CSR grid.
        let r = pts(&[0.0, 1.2, 2.1, 1.1], &[0.0, 0.2, 0.1, 2.2]);
        let l = pts(&[0.05, 1.15], &[0.05, 0.15]);
        let grid = PointGrid::build(r.iter().copied(), cell_deg).expect("valid cell_deg");
        assert!(
            matches!(grid.buckets, Buckets::Dense(_)),
            "expected dense path"
        );
        assert_eq!(within_pairs(&grid, &l, meters), brute(&l, &r, meters));

        // Two clusters at opposite ends of the globe -> sparse hash fallback.
        let r = pts(&[0.0, 0.001, 70.0, 70.001], &[0.0, 0.001, 170.0, 170.001]);
        let l = pts(&[0.0005, 70.0005], &[0.0005, 170.0005]);
        let grid = PointGrid::build(r.iter().copied(), cell_deg).expect("valid cell_deg");
        assert!(
            matches!(grid.buckets, Buckets::Sparse(_)),
            "expected sparse path"
        );
        assert_eq!(within_pairs(&grid, &l, meters), brute(&l, &r, meters));
    }

    #[test]
    fn k_nearest_returns_k_closest_sorted() {
        let points = pts(&[0.0, 0.1, 0.2, 5.0], &[0.0, 0.0, 0.0, 0.0]);
        let probe = pts(&[0.0], &[0.0]);
        let knn = k_nearest(probe.iter().copied(), &points, 2);
        assert_eq!(knn.len(), 2);
        // Two closest to (0,0) are payloads 0 (d=0) then 1, in ascending distance.
        assert_eq!((knn[0].0, knn[0].1), (0, 0));
        assert_eq!((knn[1].0, knn[1].1), (0, 1));
        assert!(knn[0].2 <= knn[1].2);
    }

    #[test]
    fn within_is_complete_when_cell_deg_is_smaller_than_the_radius() {
        // cell_deg (0.1°) far smaller than the ~0.45° radius: matches sit several
        // cells away, so a fixed 3x3 probe would miss them — `within` must size
        // its neighbourhood from the radius.
        let cell_deg = 0.1;
        let meters = 50_000.0;
        let r = pts(&[0.0, 0.3, -0.25, 0.1], &[0.0, 0.2, 0.1, -0.3]);
        let l = pts(&[0.05], &[0.05]);
        let grid = PointGrid::build(r.iter().copied(), cell_deg).expect("valid cell_deg");
        assert_eq!(within_pairs(&grid, &l, meters), brute(&l, &r, meters));
    }

    #[test]
    fn build_rejects_invalid_cell_deg() {
        let p = pts(&[0.0], &[0.0]);
        assert!(PointGrid::build(p.iter().copied(), 0.0).is_none());
        assert!(PointGrid::build(p.iter().copied(), -1.0).is_none());
        assert!(PointGrid::build(p.iter().copied(), f64::NAN).is_none());
    }
}
