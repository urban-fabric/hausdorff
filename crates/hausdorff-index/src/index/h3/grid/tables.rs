/// CCW ring traversal directions around {1,0,0} (J, JK, K, IK, I, IJ).
pub(super) const RING_DIRECTIONS: [u8; 6] = [2, 3, 1, 5, 4, 6];

/// New digit when traversing (current digit → direction → new digit), Class II.
#[rustfmt::skip]
pub(super) const NEW_DIRECTION_II: [[u8; 7]; 7] = [
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
pub(super) const NEW_ADJUSTMENT_II: [[u8; 7]; 7] = [
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
pub(super) const NEW_DIRECTION_III: [[u8; 7]; 7] = [
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
pub(super) const NEW_ADJUSTMENT_III: [[u8; 7]; 7] = [
    [0,0,0,0,0,0,0],
    [0,1,0,3,0,1,0],
    [0,0,2,2,0,0,6],
    [0,3,2,3,0,0,0],
    [0,0,0,0,4,5,4],
    [0,1,0,0,5,5,0],
    [0,0,6,0,4,0,6],
];

/// Origin leading digit -> index leading digit -> 60° CW rotations. `0xff` =
/// invalid (the deleted K axis).
#[rustfmt::skip]
pub(super) const PENTAGON_ROTATIONS: [[u8; 7]; 7] = [
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

pub(super) fn direction_unfoldable(origin_dir: u8, index_dir: u8) -> bool {
    let offset = origin_dir * 7 + index_dir;
    FAILED_DIRECTIONS & (1 << offset) != 0
}

/// Reverse base cell direction -> leading index digit -> 60° CCW rotations,
/// undoing the warp `PENTAGON_ROTATIONS` introduces when the origin is on a
/// pentagon. `0xff` = invalid (the deleted K axis).
#[rustfmt::skip]
pub(super) const PENTAGON_ROTATIONS_REVERSE: [[u8; 7]; 7] = [
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
pub(super) const PENTAGON_ROTATIONS_REVERSE_NONPOLAR: [[u8; 7]; 7] = [
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
pub(super) const PENTAGON_ROTATIONS_REVERSE_POLAR: [[u8; 7]; 7] = [
    [0,    0,    0,    0,    0,    0,    0],
    [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff],
    [0,    1,    1,    1,    1,    1,    1],
    [0,    1,    0,    0,    0,    1,    0],
    [0,    1,    0,    0,    1,    1,    1],
    [0,    1,    0,    5,    1,    1,    0],
    [0,    1,    1,    0,    1,    1,    1],
];
