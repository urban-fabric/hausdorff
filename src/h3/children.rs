//! Children enumeration (aperture-7 descent, with the pentagon deleted-K-axis
//! skip). Ported from the H3 reference (h3o).

use super::index;

#[rustfmt::skip]
const HEXAGON_CHILDREN_COUNTS: [u64; 16] = [
    1, 7, 49, 343, 2401, 16_807, 117_649, 823_543, 5_764_801, 40_353_607,
    282_475_249, 1_977_326_743, 13_841_287_201, 96_889_010_407, 678_223_072_849,
    4_747_561_509_943,
];

#[rustfmt::skip]
const PENTAGON_CHILDREN_COUNTS: [u64; 16] = [
    1, 6, 41, 286, 2001, 14_006, 98_041, 686_286, 4_804_001, 33_628_006,
    235_396_041, 1_647_772_286, 11_534_406_001, 80_740_842_006, 565_185_894_041,
    3_956_301_258_286,
];

pub(super) fn children_count(parent_res: u8, target_res: u8, is_pentagon: bool) -> u64 {
    if parent_res > target_res {
        return 0;
    }
    if parent_res == target_res {
        return 1;
    }
    let n = (target_res - parent_res) as usize;
    if is_pentagon {
        PENTAGON_CHILDREN_COUNTS[n]
    } else {
        HEXAGON_CHILDREN_COUNTS[n]
    }
}

pub(super) fn hexagon_children_count(resolution_delta: u8) -> u64 {
    HEXAGON_CHILDREN_COUNTS[usize::from(resolution_delta)]
}

pub(super) fn pentagon_children_count(resolution_delta: u8) -> u64 {
    PENTAGON_CHILDREN_COUNTS[usize::from(resolution_delta)]
}

fn starting_state(bits: u64, parent_res: u8, target_res: u8) -> u64 {
    let range = target_res.saturating_sub(parent_res) as u64;
    let mut scratchpad = bits;
    if range != 0 {
        let mask = (1u64 << (range * 3)) - 1;
        let offset = (index::MAX_RESOLUTION - target_res) as u64 * 3;
        scratchpad &= !(mask << offset);
        scratchpad = index::set_resolution(scratchpad, target_res);
    }
    scratchpad
}

/// Iterator over a cell's children at a target resolution.
pub(super) struct Children {
    parent_res: u8,
    target_res: u8,
    scratchpad: u64,
    /// Resolution at/below which the first K-axis child is skipped (pentagons),
    /// or -1 for hexagons.
    skip_count: i16,
    count: u64,
}

impl Children {
    pub(super) fn new(bits: u64, parent_res: u8, is_pentagon: bool, target_res: u8) -> Self {
        Self {
            parent_res,
            target_res,
            scratchpad: starting_state(bits, parent_res, target_res),
            skip_count: if is_pentagon {
                i16::from(target_res)
            } else {
                -1
            },
            count: children_count(parent_res, target_res, is_pentagon),
        }
    }

    fn next_direction(&mut self, res: u8) -> u8 {
        let one = 1u64 << ((index::MAX_RESOLUTION - res) as u64 * 3);
        self.scratchpad = self.scratchpad.wrapping_add(one);
        index::get_direction(self.scratchpad, res)
    }
}

impl Iterator for Children {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        if self.count == 0 {
            return None;
        }
        let index = self.scratchpad;
        self.count -= 1;

        if self.count != 0 {
            for resolution in (self.parent_res..=self.target_res).rev() {
                let direction = self.next_direction(resolution);

                // First K-axis child of each resolution is skipped for pentagons.
                if self.skip_count == i16::from(resolution) && direction == 1 {
                    self.next_direction(resolution);
                    self.skip_count -= 1;
                }

                // Digit exhausted (rolled past IJ to the invalid 7): carry up.
                if direction > 6 {
                    self.scratchpad = index::clr_direction(self.scratchpad, resolution);
                    continue;
                }
                break;
            }
        }
        Some(index)
    }
}
