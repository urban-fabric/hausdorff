//! The H3 index u64 bit codec.
//!
//! An H3 index packs, from the high bits down: a reserved high bit (0), a 4-bit
//! mode (1 = cell), 3 reserved bits (0 for cells), a 4-bit resolution, a 7-bit
//! base cell, then 15 three-bit digits (one per resolution level). Digits
//! beyond the index's resolution are set to the sentinel `7`.
//!
//! Layout (bit offsets):
//! ```text
//!  63        high bit (must be 0)
//!  59..=62   mode (cell = 1)
//!  56..=58   reserved (0 for cells)
//!  52..=55   resolution (0..=15)
//!  45..=51   base cell (0..=121)
//!  0..=44    15 digits x 3 bits  (digit r at bits (15-r)*3 .. +2)
//! ```
//!
//! These accessors are shared across the whole `h3` module — the FaceIJK
//! projection, grid traversal, and child enumeration all manipulate raw bits —
//! so [`super`] re-exports the subset reached via the `index::` path.

const MODE_OFFSET: u64 = 59;
const RESERVED_OFFSET: u64 = 56;
const RESOLUTION_OFFSET: u64 = 52;
const BASE_CELL_OFFSET: u64 = 45;

pub(in crate::index::h3) const CELL_MODE: u64 = 1;
pub(in crate::index::h3) const DIRECTED_EDGE_MODE: u64 = 2;
pub(in crate::index::h3) const VERTEX_MODE: u64 = 4;
pub const MAX_RESOLUTION: u8 = 15;
/// Sentinel digit value stored in unused (beyond-resolution) digit slots.
pub(in crate::index::h3) const INVALID_DIGIT: u8 = 7;
/// The K-axis digit — the "deleted subsequence" leading digit for pentagons.
pub(in crate::index::h3) const K_AXIS_DIGIT: u8 = 1;

#[inline]
pub(in crate::index::h3) const fn get_mode(bits: u64) -> u8 {
    ((bits >> MODE_OFFSET) & 0xf) as u8
}

#[inline]
pub(in crate::index::h3) const fn get_resolution(bits: u64) -> u8 {
    ((bits >> RESOLUTION_OFFSET) & 0xf) as u8
}

#[inline]
pub(in crate::index::h3) const fn get_base_cell(bits: u64) -> u8 {
    ((bits >> BASE_CELL_OFFSET) & 0x7f) as u8
}

#[inline]
pub(in crate::index::h3) const fn get_reserved_bits(bits: u64) -> u8 {
    ((bits >> RESERVED_OFFSET) & 0x7) as u8
}

/// Returns the 3-bit digit at resolution `res` (1..=15).
#[inline]
pub(in crate::index::h3) const fn get_digit(bits: u64, res: u8) -> u8 {
    debug_assert!(res >= 1 && res <= MAX_RESOLUTION);
    ((bits >> ((MAX_RESOLUTION - res) as u64 * 3)) & 0x7) as u8
}

/// A cell index with mode=cell, resolution 0, base cell 0, and all 15 digits
/// set to the sentinel `7` — the starting point for building an index.
pub(in crate::index::h3) const DEFAULT_CELL_INDEX: u64 =
    ((CELL_MODE) << MODE_OFFSET) | ((1u64 << (3 * MAX_RESOLUTION as u64)) - 1);

#[inline]
const fn digit_shift(res: u8) -> u64 {
    (MAX_RESOLUTION - res) as u64 * 3
}

/// `get_digit`, named for the projection code.
#[inline]
pub(in crate::index::h3) const fn get_direction(bits: u64, res: u8) -> u8 {
    get_digit(bits, res)
}

#[inline]
pub(in crate::index::h3) const fn set_resolution(bits: u64, res: u8) -> u64 {
    (bits & !(0xf << RESOLUTION_OFFSET)) | ((res as u64) << RESOLUTION_OFFSET)
}

#[inline]
pub(in crate::index::h3) const fn set_mode(bits: u64, mode: u64) -> u64 {
    (bits & !(0xf << MODE_OFFSET)) | (mode << MODE_OFFSET)
}

#[inline]
pub(in crate::index::h3) const fn set_reserved_bits(bits: u64, value: u8) -> u64 {
    (bits & !(0x7 << RESERVED_OFFSET)) | ((value as u64) << RESERVED_OFFSET)
}

#[inline]
pub(in crate::index::h3) const fn set_base_cell(bits: u64, base_cell: u8) -> u64 {
    (bits & !(0x7f << BASE_CELL_OFFSET)) | ((base_cell as u64) << BASE_CELL_OFFSET)
}

#[inline]
pub(in crate::index::h3) const fn set_direction(bits: u64, digit: u8, res: u8) -> u64 {
    let shift = digit_shift(res);
    (bits & !(0x7 << shift)) | ((digit as u64) << shift)
}

#[inline]
pub(in crate::index::h3) const fn clr_direction(bits: u64, res: u8) -> u64 {
    bits & !(0x7 << digit_shift(res))
}

/// Rotates a single direction digit `count` times by 60°.
#[inline]
pub(in crate::index::h3) const fn rotate_direction(digit: u8, count: usize, ccw: bool) -> u8 {
    rotate_digit(digit, count, ccw)
}

/// Sets all unused (beyond-`res`) digit slots to the sentinel `7`.
#[inline]
pub(in crate::index::h3) const fn set_unused(bits: u64, res: u8) -> u64 {
    bits | ((1u64 << digit_shift(res)) - 1)
}

/// The first non-center digit (1..=6) reading from resolution 1, or `None`.
pub(in crate::index::h3) const fn first_axe(bits: u64) -> Option<u8> {
    let res = get_resolution(bits);
    let mut r = 1u8;
    while r <= res {
        let d = get_digit(bits, r);
        if d != 0 {
            return Some(d);
        }
        r += 1;
    }
    None
}

/// A single 60° digit rotation (`rotate60_once`).
#[inline]
const fn rotate_digit_once(digit: u8, ccw: bool) -> u8 {
    match digit {
        1 => {
            if ccw {
                5
            } else {
                3
            }
        } // K -> IK / JK
        2 => {
            if ccw {
                3
            } else {
                6
            }
        } // J -> JK / IJ
        3 => {
            if ccw {
                1
            } else {
                2
            }
        } // JK -> K / J
        4 => {
            if ccw {
                6
            } else {
                5
            }
        } // I -> IJ / IK
        5 => {
            if ccw {
                4
            } else {
                1
            }
        } // IK -> I / K
        6 => {
            if ccw {
                2
            } else {
                4
            }
        } // IJ -> J / I
        _ => 0, // Center
    }
}

/// `count` 60° digit rotations.
const fn rotate_digit(digit: u8, count: usize, ccw: bool) -> u8 {
    const CCW_SEQ: [u8; 6] = [1, 5, 4, 6, 2, 3];
    const CW_SEQ: [u8; 6] = [1, 3, 2, 6, 4, 5];
    if count == 0 || digit == 0 {
        return digit;
    }
    let offset = match digit {
        1 => 0,
        2 => {
            if ccw {
                4
            } else {
                2
            }
        }
        3 => {
            if ccw {
                5
            } else {
                1
            }
        }
        4 => {
            if ccw {
                2
            } else {
                4
            }
        }
        5 => {
            if ccw {
                1
            } else {
                5
            }
        }
        6 => 3,
        _ => return digit,
    };
    let index = (count + offset) % 6;
    if ccw { CCW_SEQ[index] } else { CW_SEQ[index] }
}

/// Rotates a whole index `count` times by 60°.
pub(in crate::index::h3) fn rotate60(mut bits: u64, count: usize, ccw: bool) -> u64 {
    let res = get_resolution(bits);
    let mut r = 1u8;
    while r <= res {
        let d = get_direction(bits, r);
        let d = if count == 1 {
            rotate_digit_once(d, ccw)
        } else {
            rotate_digit(d, count, ccw)
        };
        bits = set_direction(bits, d, r);
        r += 1;
    }
    bits
}

/// Rotates a pentagon index by 60°, accounting for the missing K-axis sequence.
pub(in crate::index::h3) fn pentagon_rotate60(mut bits: u64, ccw: bool) -> u64 {
    let res = get_resolution(bits);
    if res == 0 {
        return bits;
    }
    // Direction that would rotate to the (deleted) K axis.
    let target = if ccw { 3 } else { 5 }; // JK.axe() / IK.axe()
    let double = first_axe(bits) == Some(target);
    let mut r = 1u8;
    while r <= res {
        let d = get_direction(bits, r);
        let d = if double {
            rotate_digit(d, 2, ccw)
        } else {
            rotate_digit_once(d, ccw)
        };
        bits = set_direction(bits, d, r);
        r += 1;
    }
    bits
}

/// True iff any used digit (1..=res) is non-center.
pub(in crate::index::h3) const fn has_nonzero_digit(bits: u64) -> bool {
    let res = get_resolution(bits);
    let mut r = 1u8;
    while r <= res {
        if get_digit(bits, r) != 0 {
            return true;
        }
        r += 1;
    }
    false
}
