//! Index validity rules — the `is_valid_*` predicates backing
//! [`H3Cell::try_new`](super::H3Cell::try_new) and its directed-edge / vertex
//! counterparts. Each mirrors the corresponding H3 reference check.

use super::bits::{
    CELL_MODE, DIRECTED_EDGE_MODE, INVALID_DIGIT, K_AXIS_DIGIT, MAX_RESOLUTION, VERTEX_MODE,
    get_base_cell, get_digit, get_mode, get_reserved_bits, get_resolution, set_mode,
    set_reserved_bits,
};
use super::{H3Cell, H3Vertex};
use crate::index::h3::base_cells;

/// Mirrors H3 `isValidCell`: high bit clear, cell mode, reserved bits clear,
/// base cell in range, used digits 0..=6 (no `7`), the pentagon deleted-K-axis
/// rule, and unused digits all `7`.
pub(super) fn is_valid_cell(bits: u64) -> bool {
    if (bits >> 63) & 1 != 0 {
        return false;
    }
    if get_mode(bits) as u64 != CELL_MODE {
        return false;
    }
    if get_reserved_bits(bits) != 0 {
        return false;
    }
    let base_cell = get_base_cell(bits);
    if base_cell >= base_cells::COUNT {
        return false;
    }
    let res = get_resolution(bits);

    let is_pent_base = base_cells::is_pentagon(base_cell);
    let mut found_first_nonzero = false;
    let mut r = 1u8;
    while r <= res {
        let d = get_digit(bits, r);
        if !found_first_nonzero && d != 0 {
            found_first_nonzero = true;
            // Pentagons have no K-axis child (deleted subsequence).
            if is_pent_base && d == K_AXIS_DIGIT {
                return false;
            }
        }
        if d == INVALID_DIGIT {
            return false;
        }
        r += 1;
    }
    // Unused digits must be the sentinel 7.
    while r <= MAX_RESOLUTION {
        if get_digit(bits, r) != INVALID_DIGIT {
            return false;
        }
        r += 1;
    }
    true
}

pub(super) fn is_valid_directed_edge(bits: u64) -> bool {
    if get_mode(bits) as u64 != DIRECTED_EDGE_MODE {
        return false;
    }
    let direction = get_reserved_bits(bits);
    if !(1..=6).contains(&direction) {
        return false;
    }

    let origin = set_reserved_bits(set_mode(bits, CELL_MODE), 0);
    let Some(origin) = H3Cell::try_new(origin) else {
        return false;
    };

    !(origin.is_pentagon() && direction == K_AXIS_DIGIT)
}

pub(super) fn is_valid_vertex(bits: u64) -> bool {
    if get_mode(bits) as u64 != VERTEX_MODE {
        return false;
    }

    let owner_bits = set_reserved_bits(set_mode(bits, CELL_MODE), 0);
    let Some(owner) = H3Cell::try_new(owner_bits) else {
        return false;
    };
    let vertex = get_reserved_bits(bits);
    owner.vertex(vertex).map(H3Vertex::as_u64) == Some(bits)
}
