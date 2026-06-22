use std::sync::LazyLock;

pub(super) const LOOKUP_BITS: u32 = 4;
pub(super) const SWAP_MASK: i64 = 0x01;
pub(super) const INVERT_MASK: i64 = 0x02;
const LOOKUP_LEN: usize = 1 << (2 * LOOKUP_BITS + 2);

const POS_TO_IJ: [[i64; 4]; 4] = [[0, 1, 3, 2], [0, 2, 3, 1], [3, 2, 0, 1], [3, 1, 0, 2]];
const POS_TO_ORIENTATION: [i64; 4] = [SWAP_MASK, 0, 0, INVERT_MASK | SWAP_MASK];

pub(super) struct LookupTables {
    pub(super) pos: [i64; LOOKUP_LEN],
    pub(super) ij: [i64; LOOKUP_LEN],
}

/// `lookup_pos` maps `(i, j, orientation)` to the Hilbert curve position;
/// `lookup_ij` is the inverse. Built once, on first use.
pub(super) static LOOKUP: LazyLock<LookupTables> = LazyLock::new(|| {
    let mut t = LookupTables {
        pos: [0; LOOKUP_LEN],
        ij: [0; LOOKUP_LEN],
    };
    init_lookup_cell(0, 0, 0, 0, 0, 0, &mut t);
    init_lookup_cell(0, 0, 0, SWAP_MASK, 0, SWAP_MASK, &mut t);
    init_lookup_cell(0, 0, 0, INVERT_MASK, 0, INVERT_MASK, &mut t);
    init_lookup_cell(
        0,
        0,
        0,
        SWAP_MASK | INVERT_MASK,
        0,
        SWAP_MASK | INVERT_MASK,
        &mut t,
    );
    t
});

#[allow(clippy::too_many_arguments)]
fn init_lookup_cell(
    level: u32,
    i: i64,
    j: i64,
    orig_orientation: i64,
    pos: i64,
    orientation: i64,
    t: &mut LookupTables,
) {
    if level == LOOKUP_BITS {
        let ij = (i << LOOKUP_BITS) + j;
        t.pos[((ij << 2) + orig_orientation) as usize] = (pos << 2) + orientation;
        t.ij[((pos << 2) + orig_orientation) as usize] = (ij << 2) + orientation;
        return;
    }
    let level = level + 1;
    let i = i << 1;
    let j = j << 1;
    let pos = pos << 2;
    let r = POS_TO_IJ[orientation as usize];
    init_lookup_cell(
        level,
        i + (r[0] >> 1),
        j + (r[0] & 1),
        orig_orientation,
        pos,
        orientation ^ POS_TO_ORIENTATION[0],
        t,
    );
    init_lookup_cell(
        level,
        i + (r[1] >> 1),
        j + (r[1] & 1),
        orig_orientation,
        pos + 1,
        orientation ^ POS_TO_ORIENTATION[1],
        t,
    );
    init_lookup_cell(
        level,
        i + (r[2] >> 1),
        j + (r[2] & 1),
        orig_orientation,
        pos + 2,
        orientation ^ POS_TO_ORIENTATION[2],
        t,
    );
    init_lookup_cell(
        level,
        i + (r[3] >> 1),
        j + (r[3] & 1),
        orig_orientation,
        pos + 3,
        orientation ^ POS_TO_ORIENTATION[3],
        t,
    );
}
