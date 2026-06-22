//! FaceIJK: a face number + IJK coords on that face's grid, and the conversions
//! to/from a cell index (with icosahedron-edge "overage" handling). Ported from
//! the H3 reference (h3o); resolution is a `usize`, Class III is `res & 1 == 1`.

use super::base_cells;
use super::coord::{CCW, CW, CoordIJK, Vec2d, Vec3d};
use super::face::{self, IJ, JK, KI};
use super::index;

mod boundary;
mod tables;

use tables::FACE_IJK_BASE_CELLS;

static ZERO: CoordIJK = CoordIJK::new(0, 0, 0);

/// Overage classification when adjusting to the correct icosahedron face.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub(super) enum Overage {
    /// On the original face.
    None,
    /// On a face edge (only on substrate grids).
    FaceEdge,
    /// Interior of a new face.
    NewFace,
}

/// A base cell and its 60° CCW rotation count.
#[derive(Clone, Copy)]
pub(super) struct Rotation {
    pub(super) base_cell: u8,
    pub(super) count: u8,
}

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub(super) struct FaceIJK {
    pub(super) face: usize,
    pub(super) coord: CoordIJK,
}

#[inline]
const fn is_class3(res: usize) -> bool {
    res & 1 == 1
}

impl FaceIJK {
    pub(super) const fn new(face: usize, coord: CoordIJK) -> Self {
        Self { face, coord }
    }

    /// Home FaceIJK for a base cell.
    fn from_base_cell(base_cell: u8) -> Self {
        Self::new(
            base_cells::home_face(base_cell),
            base_cells::home_coord(base_cell),
        )
    }

    /// Encodes an n-vector to the FaceIJK of the containing cell at `res`.
    pub(super) fn from_vec3d(value: Vec3d, res: usize) -> Self {
        let (face, distance) = value.closest_face();
        let coord = Vec2d::from_vec3d(value, res, face, distance).into();
        Self::new(face, coord)
    }

    /// As [`from_vec3d`](Self::from_vec3d) but via the faster gnomonic radial
    /// factor ([`Vec2d::from_vec3d_fast`](super::coord::Vec2d::from_vec3d_fast)); equal
    /// except at the last ULPs, so a result can shift to an adjacent cell only for a
    /// near-boundary point.
    pub(super) fn from_vec3d_fast(value: Vec3d, res: usize) -> Self {
        let (face, distance) = value.closest_face();
        let coord = Vec2d::from_vec3d_fast(value, res, face, distance).into();
        Self::new(face, coord)
    }

    /// 60° CCW rotation count into the base cell's coordinate system here.
    fn base_cell_rotation(&self) -> Rotation {
        let i = self.coord.i() as usize;
        let j = self.coord.j() as usize;
        let k = self.coord.k() as usize;
        FACE_IJK_BASE_CELLS[self.face][i][j][k]
    }

    /// Converts a FaceIJK at `res` to a raw cell index.
    pub(super) fn to_cell_bits(mut self, res: usize) -> u64 {
        let mut bits = index::set_resolution(index::DEFAULT_CELL_INDEX, res as u8);

        if res == 0 {
            let rotation = self.base_cell_rotation();
            return index::set_base_cell(bits, rotation.base_cell);
        }

        self.coord = directions_bits_from_ijk(self.coord, &mut bits, res);

        let rotation = self.base_cell_rotation();
        bits = index::set_base_cell(bits, rotation.base_cell);

        if base_cells::is_pentagon(rotation.base_cell) {
            if index::first_axe(bits) == Some(1) {
                if base_cells::is_cw_offset(rotation.base_cell, self.face) {
                    bits = index::rotate60(bits, 1, CW);
                } else {
                    bits = index::rotate60(bits, 1, CCW);
                }
            }
            for _ in 0..rotation.count {
                bits = index::pentagon_rotate60(bits, CCW);
            }
        } else {
            bits = index::rotate60(bits, rotation.count as usize, CCW);
        }

        bits
    }

    /// FaceIJK for a cell index (home base cell + descend digits), plus whether
    /// an overage is possible.
    fn from_bits(bits: u64, res: usize, base_cell: u8) -> (Self, bool) {
        let mut fijk = Self::from_base_cell(base_cell);
        let possible_overage =
            base_cells::is_pentagon(base_cell) || (res != 0 || fijk.coord != ZERO);

        for r in 1..=res {
            if is_class3(r) {
                fijk.coord = fijk.coord.down_aperture7::<CCW>();
            } else {
                fijk.coord = fijk.coord.down_aperture7::<CW>();
            }
            let digit = index::get_direction(bits, r as u8);
            fijk.coord = fijk.coord.neighbor(digit);
        }
        (fijk, possible_overage)
    }

    /// Adjusts this FaceIJK in place so it is relative to the correct face.
    fn adjust_overage_class2(
        &mut self,
        class2_res: usize,
        is_pent4: bool,
        substrate: bool,
    ) -> Overage {
        let factor = if substrate { 3 } else { 1 };
        let face = self.face;
        let dimension = self.coord.i() + self.coord.j() + self.coord.k();
        let max_dim = MAX_DIM_BY_CII_RES[class2_res] * factor;

        if substrate && dimension == max_dim {
            return Overage::FaceEdge;
        }
        if dimension > max_dim {
            let orientation = if self.coord.k() > 0 {
                if self.coord.j() > 0 {
                    face::NEIGHBORS[face][JK]
                } else {
                    if is_pent4 {
                        let origin = CoordIJK::new(max_dim, 0, 0);
                        let tmp = (self.coord - origin).rotate60::<CW>();
                        self.coord = tmp + origin;
                    }
                    face::NEIGHBORS[face][KI]
                }
            } else {
                face::NEIGHBORS[face][IJ]
            };
            self.face = orientation.face as usize;

            for _ in 0..orientation.ccw_rot60 {
                self.coord = self.coord.rotate60::<CCW>();
            }
            let mut trans = orientation.translate;
            let unit_scale = UNIT_SCALE_BY_CII_RES[class2_res] * factor;
            trans *= unit_scale;
            self.coord = (self.coord + trans).normalize();

            if substrate && self.coord.i() + self.coord.j() + self.coord.k() == max_dim {
                return Overage::FaceEdge;
            }
            return Overage::NewFace;
        }
        Overage::None
    }

    /// FaceIJK for a cell index, with overage applied (the center-on-face form).
    pub(super) fn from_cell(mut bits: u64, base_cell: u8, res: usize) -> Self {
        if base_cells::is_pentagon(base_cell) && index::first_axe(bits) == Some(5) {
            bits = index::rotate60(bits, 1, CW);
        }

        let (mut fijk, overage) = Self::from_bits(bits, res, base_cell);
        if !overage {
            return fijk;
        }

        let original_coord = fijk.coord;
        let class3 = is_class3(res);
        let class2_res = if class3 {
            fijk.coord = fijk.coord.down_aperture7::<CW>();
            res + 1
        } else {
            res
        };

        let is_pent4 = base_cells::is_pentagon(base_cell) && index::first_axe(bits) == Some(4);

        if fijk.adjust_overage_class2(class2_res, is_pent4, false) != Overage::None {
            if base_cells::is_pentagon(base_cell) {
                while fijk.adjust_overage_class2(class2_res, false, false) != Overage::None {}
            }
            if class3 {
                fijk.coord = fijk.coord.up_aperture7::<CW>();
            }
        } else if class3 {
            fijk.coord = original_coord;
        }

        fijk
    }

    /// The cell's digit sequence resolved into a neutral base-cell coordinate
    /// space (rooted at base cell 2, whose home IJK is the origin). The face is
    /// irrelevant here; only the IJK is used (for local-IJK / grid distance).
    pub(super) fn local_ijk_coord(bits: u64, res: usize) -> CoordIJK {
        const BASE_CELL_FOR_LOCAL: u8 = 2;
        Self::from_bits(bits, res, BASE_CELL_FOR_LOCAL).0.coord
    }

    /// The cell center as geographic coordinates in radians `(lat, lng)`.
    pub(super) fn to_center_lat_lng_rad(self, res: usize) -> (f64, f64) {
        let v = Vec3d::from_vec2d(Vec2d::from(self.coord), self.face, res, false);
        v.to_lat_lng_rad()
    }

    #[allow(clippy::needless_range_loop)]
    pub(super) fn base_cell_rotation_count(base_cell: u8, face: usize) -> Option<u8> {
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let rotation = FACE_IJK_BASE_CELLS[face][i][j][k];
                    if rotation.base_cell == base_cell {
                        return Some(rotation.count);
                    }
                }
            }
        }
        None
    }
}

/// Sets the cell index direction bits from the FaceIJK, finest resolution up,
/// returning the IJK adjusted into the base cell's coordinate system.
fn directions_bits_from_ijk(mut ijk: CoordIJK, bits: &mut u64, res: usize) -> CoordIJK {
    for r in (1..=res).rev() {
        let last_ijk = ijk;
        let last_center = if is_class3(r) {
            ijk = ijk.up_aperture7::<CCW>();
            ijk.down_aperture7::<CCW>()
        } else {
            ijk = ijk.up_aperture7::<CW>();
            ijk.down_aperture7::<CW>()
        };
        let diff = (last_ijk - last_center).normalize();
        *bits = index::set_direction(*bits, diff.to_digit(), r as u8);
    }
    ijk
}

#[rustfmt::skip]
static MAX_DIM_BY_CII_RES: [i32; 17] = [
    2, -1, 14, -1, 98, -1, 686, -1, 4802, -1, 33614, -1, 235298, -1, 1647086, -1, 11529602,
];

#[rustfmt::skip]
static UNIT_SCALE_BY_CII_RES: [i32; 17] = [
    1, -1, 7, -1, 49, -1, 343, -1, 2401, -1, 16807, -1, 117649, -1, 823543, -1, 5764801,
];
