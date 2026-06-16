//! Hexagon grid traversal: `neighbor_rotations` (the pentagon-aware step) and a
//! safe BFS `grid_disk`. Ported from the H3 reference (h3o), using `u8`
//! direction digits (Center=0, K=1, J=2, JK=3, I=4, IK=5, IJ=6).

use std::collections::{HashSet, VecDeque};

use super::base_cells;
use super::coord::{CCW, CW, CoordIJK};
use super::faceijk::FaceIJK;
use super::index;

/// CCW ring traversal directions around {1,0,0} (J, JK, K, IK, I, IJ).
const RING_DIRECTIONS: [u8; 6] = [2, 3, 1, 5, 4, 6];

#[inline]
const fn is_class3(res: u8) -> bool {
    res & 1 == 1
}

/// New digit when traversing (current digit → direction → new digit), Class II.
#[rustfmt::skip]
const NEW_DIRECTION_II: [[u8; 7]; 7] = [
    [0,1,2,3,4,5,6],
    [1,4,3,6,5,2,0],
    [2,3,1,4,6,0,5],
    [3,6,4,5,0,1,2],
    [4,5,6,0,2,3,1],
    [5,2,0,1,3,6,4],
    [6,0,5,2,1,4,3],
];
/// New coarser-level move (current digit → direction → adjustment), Class II.
#[rustfmt::skip]
const NEW_ADJUSTMENT_II: [[u8; 7]; 7] = [
    [0,0,0,0,0,0,0],
    [0,1,0,1,0,5,0],
    [0,0,2,3,0,0,2],
    [0,1,3,3,0,0,0],
    [0,0,0,0,4,4,6],
    [0,5,0,0,4,5,0],
    [0,0,2,0,6,0,6],
];
/// New digit, Class III.
#[rustfmt::skip]
const NEW_DIRECTION_III: [[u8; 7]; 7] = [
    [0,1,2,3,4,5,6],
    [1,2,3,4,5,6,0],
    [2,3,4,5,6,0,1],
    [3,4,5,6,0,1,2],
    [4,5,6,0,1,2,3],
    [5,6,0,1,2,3,4],
    [6,0,1,2,3,4,5],
];
/// New coarser-level move, Class III.
#[rustfmt::skip]
const NEW_ADJUSTMENT_III: [[u8; 7]; 7] = [
    [0,0,0,0,0,0,0],
    [0,1,0,3,0,1,0],
    [0,0,2,2,0,0,6],
    [0,3,2,3,0,0,0],
    [0,0,0,0,4,5,4],
    [0,1,0,0,5,5,0],
    [0,0,6,0,4,0,6],
];

/// The cell neighboring `origin` in `direction`, plus the updated rotation
/// count. Returns `None` only when stepping into a pentagon's deleted K axis.
pub(super) fn neighbor_rotations(
    origin: u64,
    direction: u8,
    mut rotations: u8,
) -> Option<(u64, u8)> {
    let mut current = origin;
    let mut direction = index::rotate_direction(direction, rotations as usize, CCW);
    let origin_res = index::get_resolution(origin);

    let mut res = origin_res;
    let mut r = origin_res;
    while r >= 1 {
        let from = index::get_direction(current, r) as usize;
        let to = direction as usize;
        let next_direction = if is_class3(r) {
            current = index::set_direction(current, NEW_DIRECTION_II[from][to], r);
            NEW_ADJUSTMENT_II[from][to]
        } else {
            current = index::set_direction(current, NEW_DIRECTION_III[from][to], r);
            NEW_ADJUSTMENT_III[from][to]
        };
        if next_direction == 0 {
            break;
        }
        direction = next_direction;
        res = r - 1;
        r -= 1;
    }

    let old_base_cell = index::get_base_cell(origin);
    let new_rotations = if res == 0 {
        if let Some(bc) = base_cells::neighbor(old_base_cell, direction) {
            current = index::set_base_cell(current, bc);
            base_cells::neighbor_rotation(old_base_cell, direction)
        } else {
            let bc = base_cells::neighbor(old_base_cell, 5)?;
            current = index::set_base_cell(current, bc);
            current = index::rotate60(current, 1, CCW);
            rotations += 1;
            base_cells::neighbor_rotation(old_base_cell, 5)
        }
    } else {
        0
    };

    let new_base_cell = index::get_base_cell(current);
    if base_cells::is_pentagon(new_base_cell) {
        let mut already_adjusted = false;
        if index::first_axe(current) == Some(1) {
            if old_base_cell == new_base_cell {
                let old_leading = index::first_axe(origin);
                old_leading?;
                if old_leading == Some(3) {
                    current = index::rotate60(current, 1, CCW);
                    rotations += 1;
                } else if old_leading == Some(5) {
                    current = index::rotate60(current, 1, CW);
                    rotations += 5;
                }
            } else {
                if base_cells::is_cw_offset(new_base_cell, base_cells::home_face(old_base_cell)) {
                    current = index::rotate60(current, 1, CW);
                }
                already_adjusted = true;
            }
        }

        for _ in 0..new_rotations {
            current = index::pentagon_rotate60(current, CCW);
        }

        if old_base_cell != new_base_cell {
            let dir = index::first_axe(current);
            if base_cells::is_polar_pentagon(new_base_cell) {
                rotations += u8::from(old_base_cell != 118 && old_base_cell != 8 && dir != Some(3));
            } else {
                rotations += u8::from(!already_adjusted && dir == Some(5));
            }
        }
    } else if new_rotations != 0 {
        current = index::rotate60(current, new_rotations as usize, CCW);
    }

    Some((current, (rotations + new_rotations) % 6))
}

/// All cells within grid distance `k` of `origin` (the filled grid disk,
/// including `origin`). Safe BFS — robust around pentagons.
///
/// Proof-backed by `riemann` `Hausdorff/H3Grid.gridDistance_stepPath_le` (at the
/// ideal-`CoordIJK` level): a BFS path of `n` unit ring-steps stays within grid
/// distance `n`, so the BFS ring index `k` soundly bounds the grid distance and
/// `grid_disk k` ⊆ the radius-`k` disk. (The cell-id/pentagon traversal itself —
/// `neighbor_rotations` — is graph plumbing, out of scope.)
pub(super) fn grid_disk(origin: u64, k: u32) -> Vec<u64> {
    grid_disk_distances(origin, k)
        .into_iter()
        .map(|(cell, _distance)| cell)
        .collect()
}

/// All cells within grid distance `k` of `origin`, paired with their BFS grid
/// distance from `origin`.
pub(super) fn grid_disk_distances(origin: u64, k: u32) -> Vec<(u64, u32)> {
    let mut seen: HashSet<u64> = HashSet::new();
    let mut candidates: VecDeque<(u64, u32)> = VecDeque::new();
    let mut out = Vec::new();
    candidates.push_back((origin, 0));

    while let Some((cell, ring)) = candidates.pop_front() {
        if ring > k || seen.contains(&cell) {
            continue;
        }
        if ring < k {
            for &d in &RING_DIRECTIONS {
                if let Some((n, _)) = neighbor_rotations(cell, d, 0) {
                    candidates.push_back((n, ring + 1));
                }
            }
        }
        seen.insert(cell);
        out.push((cell, ring));
    }
    out
}

/// Origin leading digit -> index leading digit -> 60° CW rotations. `0xff` =
/// invalid (the deleted K axis).
#[rustfmt::skip]
const PENTAGON_ROTATIONS: [[u8; 7]; 7] = [
    [0,    0xff, 0,    0,    0,    0,    0],
    [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff],
    [0,    0xff, 0,    0,    0,    1,    0],
    [0,    0xff, 0,    0,    1,    1,    0],
    [0,    0xff, 0,    5,    0,    0,    0],
    [0,    0xff, 5,    5,    0,    0,    0],
    [0,    0xff, 0,    0,    0,    0,    0],
];

/// Direction pairs (origin_dir*7 + index_dir) that cannot be unfolded across a
/// pentagon.
const FAILED_DIRECTIONS: u64 = 0b0101000_1000100_0001100_1010000_0110000_0000000_0000000;

fn direction_unfoldable(origin_dir: u8, index_dir: u8) -> bool {
    let offset = origin_dir * 7 + index_dir;
    FAILED_DIRECTIONS & (1 << offset) != 0
}

/// The cell's IJK coordinate in `origin`'s local coordinate system (mirrors H3
/// `cellToLocalIjk`). `None` if the cells can't be unfolded into one frame.
pub(super) fn local_ijk(cell: u64, origin: u64) -> Option<CoordIJK> {
    let origin_res = index::get_resolution(origin);
    let origin_bc = index::get_base_cell(origin);
    let bc = index::get_base_cell(cell);
    let mut bits = cell;

    let mut dir = 0u8;
    let mut rev_dir;
    if origin_bc == bc {
        rev_dir = 0u8;
    } else {
        dir = base_cells::direction_to(origin_bc, bc)?;
        rev_dir = base_cells::direction_to(bc, origin_bc)?;
    }

    let origin_on_pent = base_cells::is_pentagon(origin_bc);
    let index_on_pent = base_cells::is_pentagon(bc);

    if dir != 0 {
        let base_cell_rotations = base_cells::neighbor_rotation(origin_bc, dir) as usize;
        if index_on_pent {
            for _ in 0..base_cell_rotations {
                rev_dir = if rev_dir == 5 {
                    index::rotate_direction(rev_dir, 2, CW)
                } else {
                    index::rotate_direction(rev_dir, 1, CW)
                };
                bits = index::pentagon_rotate60(bits, CW);
            }
        } else {
            rev_dir = index::rotate_direction(rev_dir, base_cell_rotations, CW);
            bits = index::rotate60(bits, base_cell_rotations, CW);
        }
    }

    let mut ijk = FaceIJK::local_ijk_coord(bits, origin_res as usize);

    if dir != 0 {
        let (pentagon_rotations, direction_rotations) = if origin_on_pent {
            let leading = index::first_axe(origin).unwrap_or(0);
            if direction_unfoldable(leading, dir) {
                return None;
            }
            let rot = PENTAGON_ROTATIONS[leading as usize][dir as usize];
            (rot, rot)
        } else if index_on_pent {
            let leading = index::first_axe(bits).unwrap_or(0);
            if direction_unfoldable(leading, rev_dir) {
                return None;
            }
            (PENTAGON_ROTATIONS[rev_dir as usize][leading as usize], 0u8)
        } else {
            (0u8, 0u8)
        };
        if pentagon_rotations == 0xff || direction_rotations == 0xff {
            return None;
        }
        for _ in 0..pentagon_rotations {
            ijk = ijk.rotate60::<CW>();
        }

        let mut offset = CoordIJK::new(0, 0, 0).neighbor(dir);
        for res in (1..=origin_res).rev() {
            offset = if is_class3(res) {
                offset.down_aperture7::<CCW>()
            } else {
                offset.down_aperture7::<CW>()
            };
        }
        for _ in 0..direction_rotations {
            offset = offset.rotate60::<CW>();
        }
        ijk = (ijk + offset).normalize();
    } else if origin_on_pent && index_on_pent {
        let o_leading = index::first_axe(origin).unwrap_or(0);
        let i_leading = index::first_axe(bits).unwrap_or(0);
        if direction_unfoldable(o_leading, i_leading) {
            return None;
        }
        let rot = PENTAGON_ROTATIONS[o_leading as usize][i_leading as usize];
        if rot == 0xff {
            return None;
        }
        for _ in 0..rot {
            ijk = ijk.rotate60::<CW>();
        }
    }

    Some(ijk)
}

/// Reverse base cell direction -> leading index digit -> 60° CCW rotations,
/// undoing the warp `PENTAGON_ROTATIONS` introduces when the origin is on a
/// pentagon. `0xff` = invalid (the deleted K axis).
#[rustfmt::skip]
const PENTAGON_ROTATIONS_REVERSE: [[u8; 7]; 7] = [
    [0,    0,    0,    0,    0,    0,    0],
    [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff],
    [0,    1,    0,    0,    0,    0,    0],
    [0,    1,    0,    0,    0,    1,    0],
    [0,    5,    0,    0,    0,    0,    0],
    [0,    5,    0,    5,    0,    0,    0],
    [0,    0,    0,    0,    0,    0,    0],
];

/// As `PENTAGON_ROTATIONS_REVERSE` for an index on a non-polar pentagon with a
/// non-pentagon origin.
#[rustfmt::skip]
const PENTAGON_ROTATIONS_REVERSE_NONPOLAR: [[u8; 7]; 7] = [
    [0,    0,    0,    0,    0,    0,    0],
    [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff],
    [0,    1,    0,    0,    0,    0,    0],
    [0,    1,    0,    0,    0,    1,    0],
    [0,    5,    0,    0,    0,    0,    0],
    [0,    1,    0,    5,    1,    1,    0],
    [0,    0,    0,    0,    0,    0,    0],
];

/// As `PENTAGON_ROTATIONS_REVERSE` for an index on a polar pentagon with a
/// non-pentagon origin.
#[rustfmt::skip]
const PENTAGON_ROTATIONS_REVERSE_POLAR: [[u8; 7]; 7] = [
    [0,    0,    0,    0,    0,    0,    0],
    [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff],
    [0,    1,    1,    1,    1,    1,    1],
    [0,    1,    0,    0,    0,    1,    0],
    [0,    1,    0,    0,    1,    1,    1],
    [0,    1,    0,    5,    1,    1,    0],
    [0,    1,    1,    0,    1,    1,    1],
];

/// Sets the cell index direction bits from a local IJK coordinate, finest
/// resolution up, returning the leftover base-cell-level IJK. `None` when the
/// coordinate is outside the representable grid (i32 overflow or a non-unit
/// per-resolution step).
fn checked_directions_bits_from_ijk(
    mut ijk: CoordIJK,
    bits: &mut u64,
    res: u8,
) -> Option<CoordIJK> {
    for r in (1..=res).rev() {
        let last_ijk = ijk;
        let last_center = if is_class3(r) {
            ijk = ijk.checked_up_aperture7::<CCW>()?;
            ijk.down_aperture7::<CCW>()
        } else {
            ijk = ijk.checked_up_aperture7::<CW>()?;
            ijk.down_aperture7::<CW>()
        };
        let digit = (last_ijk - last_center).normalize().checked_digit()?;
        *bits = index::set_direction(*bits, digit, r);
    }
    Some(ijk)
}

/// The cell at IJK coordinate `coord` in `origin`'s local coordinate system —
/// the inverse of [`local_ijk`] (mirrors H3 `localIjkToCell`). `None` when the
/// coordinate is out of range or crosses a pentagon's deleted subsequence.
pub(super) fn cell_from_local_ijk(origin: u64, coord: CoordIJK) -> Option<u64> {
    let res = index::get_resolution(origin);
    let origin_bc = index::get_base_cell(origin);
    let origin_on_pent = base_cells::is_pentagon(origin_bc);

    let mut bits = index::set_resolution(index::DEFAULT_CELL_INDEX, res);

    if res == 0 {
        let dir = coord.checked_digit()?;
        let bc = base_cells::neighbor(origin_bc, dir)?;
        return Some(index::set_base_cell(bits, bc));
    }

    let ijk = checked_directions_bits_from_ijk(coord, &mut bits, res)?;

    let mut dir = ijk.checked_digit()?;
    let mut base_cell = base_cells::neighbor(origin_bc, dir);
    let index_on_pent = base_cell.is_some_and(base_cells::is_pentagon);

    if dir != 0 {
        let mut pentagon_rotations = 0u8;
        if origin_on_pent {
            let leading = index::first_axe(origin).unwrap_or(0);
            pentagon_rotations = PENTAGON_ROTATIONS_REVERSE[leading as usize][dir as usize];
            if pentagon_rotations == 0xff {
                return None;
            }
            dir = index::rotate_direction(dir, pentagon_rotations as usize, CCW);
            if dir == 1 {
                return None;
            }
            base_cell = Some(base_cells::neighbor(origin_bc, dir)?);
        }
        let fixed_base_cell = base_cell?;
        let base_cell_rotations = base_cells::neighbor_rotation(origin_bc, dir);

        if index_on_pent {
            let rev_dir = base_cells::direction_to(fixed_base_cell, origin_bc)?;
            bits = index::rotate60(bits, base_cell_rotations as usize, CCW);

            let leading = index::first_axe(bits).unwrap_or(0);
            let pentagon_rotations = if base_cells::is_polar_pentagon(fixed_base_cell) {
                PENTAGON_ROTATIONS_REVERSE_POLAR[rev_dir as usize][leading as usize]
            } else {
                PENTAGON_ROTATIONS_REVERSE_NONPOLAR[rev_dir as usize][leading as usize]
            };
            if pentagon_rotations == 0xff {
                return None;
            }
            for _ in 0..pentagon_rotations {
                bits = index::pentagon_rotate60(bits, CCW);
            }
        } else {
            bits = index::rotate60(
                bits,
                (pentagon_rotations + base_cell_rotations) as usize,
                CCW,
            );
        }
    } else if origin_on_pent && index_on_pent {
        let origin_leading = index::first_axe(origin).unwrap_or(0);
        let index_leading = index::first_axe(bits).unwrap_or(0);
        let rotations = PENTAGON_ROTATIONS_REVERSE[origin_leading as usize][index_leading as usize];
        if rotations == 0xff {
            return None;
        }
        bits = index::rotate60(bits, rotations as usize, CCW);
    }

    if index_on_pent && index::first_axe(bits) == Some(1) {
        return None;
    }

    Some(index::set_base_cell(bits, base_cell?))
}

/// Grid distance between two cells (same resolution), or `None` if they can't
/// be placed in a common local frame (different resolution / far apart / across
/// a pentagon).
pub(super) fn grid_distance(a: u64, b: u64) -> Option<i32> {
    if index::get_resolution(a) != index::get_resolution(b) {
        return None;
    }
    let src = local_ijk(a, a)?;
    let dst = local_ijk(b, a)?;
    Some(src.distance(dst))
}

/// Direction digit from `origin` to `destination`, or `None` if the cells are
/// not grid neighbors.
pub(super) fn direction_for_neighbor(origin: u64, destination: u64) -> Option<u8> {
    let start = if index::H3Cell::try_new(origin)?.is_pentagon() {
        2
    } else {
        1
    };
    (start..=6).find(|&direction| {
        neighbor_rotations(origin, direction, 0).map(|(neighbor, _)| neighbor) == Some(destination)
    })
}
