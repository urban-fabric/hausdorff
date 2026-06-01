//! FaceIJK: a face number + IJK coords on that face's grid, and the conversions
//! to/from a cell index (with icosahedron-edge "overage" handling). Ported from
//! the H3 reference (h3o); resolution is a `usize`, Class III is `res & 1 == 1`.

use super::base_cells;
use super::coord::{CCW, CW, CoordIJK, SQRT3_2, Vec2d, Vec3d};
use super::face::{self, IJ, JK, KI};
use super::index;

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
            // Force rotation out of the missing K-axis sub-sequence.
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
        // Adjust for the pentagonal missing sequence (sub-sequence 5).
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

    /// Substrate-grid vertices of this (center) cell. Returns the adjusted
    /// (extended) resolution used for the vertices. `verts` is filled with `n`
    /// entries (6 hexagon, 5 pentagon).
    fn vertices(&mut self, res: usize, n: usize, verts: &mut [FaceIJK; 6]) -> usize {
        // Vertices CCW from the i-axis on a 33r (Class II) / 33r7r (Class III)
        // substrate grid.
        const VERTS_CII: [CoordIJK; 6] = [
            CoordIJK::new(2, 1, 0),
            CoordIJK::new(1, 2, 0),
            CoordIJK::new(0, 2, 1),
            CoordIJK::new(0, 1, 2),
            CoordIJK::new(1, 0, 2),
            CoordIJK::new(2, 0, 1),
        ];
        const VERTS_CIII: [CoordIJK; 6] = [
            CoordIJK::new(5, 4, 0),
            CoordIJK::new(1, 5, 0),
            CoordIJK::new(0, 5, 4),
            CoordIJK::new(0, 1, 5),
            CoordIJK::new(4, 0, 5),
            CoordIJK::new(5, 0, 1),
        ];

        self.coord = self.coord.down_aperture3::<CCW>();
        self.coord = self.coord.down_aperture3::<CW>();

        let (table, adjusted_res) = if is_class3(res) {
            self.coord = self.coord.down_aperture7::<CW>();
            (&VERTS_CIII, res + 1)
        } else {
            (&VERTS_CII, res)
        };

        for i in 0..n {
            verts[i].face = self.face;
            verts[i].coord = (self.coord + table[i]).normalize();
        }
        adjusted_res
    }

    /// Cell boundary as geographic `(lat, lng)` radian points (corners plus
    /// Class-III icosahedron-edge-crossing vertices). Mirrors H3 `cellToBoundary`.
    fn boundary_lat_lng_rad(self, res: usize, is_pentagon: bool) -> Vec<(f64, f64)> {
        let n = if is_pentagon { 5 } else { 6 };
        self.boundary_segment_lat_lng_rad(res, is_pentagon, 0, n as u8)
    }

    /// Boundary segment as geographic `(lat, lng)` radian points, starting at a
    /// topological vertex and returning `length` topological vertices. Segment
    /// edges crossing an icosahedron edge include an extra distortion point.
    fn boundary_segment_lat_lng_rad(
        self,
        res: usize,
        is_pentagon: bool,
        start: u8,
        length: u8,
    ) -> Vec<(f64, f64)> {
        let n = if is_pentagon { 5 } else { 6 };
        debug_assert!(usize::from(start) < n);
        debug_assert!(usize::from(length) <= n);

        let mut center = self;
        let mut verts = [FaceIJK::default(); 6];
        let adjusted_res = center.vertices(res, n, &mut verts);
        let additional_iteration = u8::from(usize::from(length) == n);

        let mut out: Vec<(f64, f64)> = Vec::with_capacity(n + 2);
        let max_dim = f64::from(MAX_DIM_BY_CII_RES[adjusted_res]);
        let edge_v0 = Vec2d::new(3.0 * max_dim, 0.0);
        let edge_v1 = Vec2d::new(-1.5 * max_dim, 3.0 * SQRT3_2 * max_dim);
        let edge_v2 = Vec2d::new(-1.5 * max_dim, -3.0 * SQRT3_2 * max_dim);

        if is_pentagon {
            // Pentagons: each vertex is overage-adjusted; Class III edges cross
            // icosahedron edges and need an interpolated vertex.
            let mut last_fijk = FaceIJK::default();
            for vert in start..(start + length + additional_iteration) {
                let mut fijk = verts[usize::from(vert) % n];
                while fijk.adjust_overage_class2(adjusted_res, false, true) == Overage::NewFace {}

                if is_class3(res) && vert > start {
                    let mut tmp = fijk;
                    let orig2d0 = Vec2d::from(last_fijk.coord);
                    let cur_to_last = face::adjacent_face_dir(tmp.face as u8, last_fijk.face as u8);
                    let orient = &face::NEIGHBORS[tmp.face][cur_to_last as usize];
                    tmp.face = orient.face as usize;
                    for _ in 0..orient.ccw_rot60 {
                        tmp.coord = tmp.coord.rotate60::<CCW>();
                    }
                    let mut trans = orient.translate;
                    trans *= UNIT_SCALE_BY_CII_RES[adjusted_res] * 3;
                    tmp.coord = (tmp.coord + trans).normalize();
                    let orig2d1 = Vec2d::from(tmp.coord);

                    let (e0, e1) =
                        match face::adjacent_face_dir(tmp.face as u8, fijk.face as u8) as usize {
                            IJ => (edge_v0, edge_v1),
                            JK => (edge_v1, edge_v2),
                            KI => (edge_v2, edge_v0),
                            _ => (edge_v0, edge_v1),
                        };
                    let inter = Vec2d::intersection((orig2d0, orig2d1), (e0, e1));
                    out.push(
                        Vec3d::from_vec2d(inter, tmp.face, adjusted_res, true).to_lat_lng_rad(),
                    );
                }
                if usize::from(vert) < usize::from(start) + n {
                    out.push(
                        Vec3d::from_vec2d(Vec2d::from(fijk.coord), fijk.face, adjusted_res, true)
                            .to_lat_lng_rad(),
                    );
                }
                last_fijk = fijk;
            }
        } else {
            let mut last_face = usize::MAX;
            let mut last_overage = Overage::None;
            for vert in start..(start + length + additional_iteration) {
                let v = usize::from(vert) % n;
                let mut fijk = verts[v];
                let overage = fijk.adjust_overage_class2(adjusted_res, false, true);

                if is_class3(res)
                    && vert > start
                    && fijk.face != last_face
                    && last_overage != Overage::FaceEdge
                {
                    let last_v = (v + 5) % n;
                    let orig2d0 = Vec2d::from(verts[last_v].coord);
                    let orig2d1 = Vec2d::from(verts[v].coord);
                    let face2 = if last_face == center.face {
                        fijk.face
                    } else {
                        last_face
                    };
                    let (e0, e1) =
                        match face::adjacent_face_dir(center.face as u8, face2 as u8) as usize {
                            IJ => (edge_v0, edge_v1),
                            JK => (edge_v1, edge_v2),
                            KI => (edge_v2, edge_v0),
                            _ => (edge_v0, edge_v1),
                        };
                    let inter = Vec2d::intersection((orig2d0, orig2d1), (e0, e1));
                    if orig2d0 != inter && orig2d1 != inter {
                        out.push(
                            Vec3d::from_vec2d(inter, center.face, adjusted_res, true)
                                .to_lat_lng_rad(),
                        );
                    }
                }
                if usize::from(vert) < usize::from(start) + n {
                    out.push(
                        Vec3d::from_vec2d(Vec2d::from(fijk.coord), fijk.face, adjusted_res, true)
                            .to_lat_lng_rad(),
                    );
                }
                last_face = fijk.face;
                last_overage = overage;
            }
        }
        out
    }

    /// Full cell-boundary points (radians) for a cell index.
    pub(super) fn cell_boundary_lat_lng_rad(
        bits: u64,
        base_cell: u8,
        res: usize,
        is_pentagon: bool,
    ) -> Vec<(f64, f64)> {
        Self::from_cell(bits, base_cell, res).boundary_lat_lng_rad(res, is_pentagon)
    }

    pub(super) fn cell_boundary_segment_lat_lng_rad(
        bits: u64,
        base_cell: u8,
        res: usize,
        is_pentagon: bool,
        start: u8,
        length: u8,
    ) -> Vec<(f64, f64)> {
        Self::from_cell(bits, base_cell, res).boundary_segment_lat_lng_rad(
            res,
            is_pentagon,
            start,
            length,
        )
    }

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

// ---- overage scale tables ---------------------------------------------------

#[rustfmt::skip]
static MAX_DIM_BY_CII_RES: [i32; 17] = [
    2, -1, 14, -1, 98, -1, 686, -1, 4802, -1, 33614, -1, 235298, -1, 1647086, -1, 11529602,
];

#[rustfmt::skip]
static UNIT_SCALE_BY_CII_RES: [i32; 17] = [
    1, -1, 7, -1, 49, -1, 343, -1, 2401, -1, 16807, -1, 117649, -1, 823543, -1, 5764801,
];

// ---- res-0 base cell lookup per face ---------------------------------------

const fn r(base_cell: u8, count: u8) -> Rotation {
    Rotation { base_cell, count }
}

#[rustfmt::skip]
static FACE_IJK_BASE_CELLS: [[[[Rotation; 3]; 3]; 3]; 20] = [
    [
        [[r(16,0),r(18,0),r(24,0)], [r(33,0),r(30,0),r(32,3)], [r(49,1),r(48,3),r(50,3)]],
        [[r(8,0),r(5,5),r(10,5)], [r(22,0),r(16,0),r(18,0)], [r(41,1),r(33,0),r(30,0)]],
        [[r(4,0),r(0,5),r(2,5)], [r(15,1),r(8,0),r(5,5)], [r(31,1),r(22,0),r(16,0)]],
    ],
    [
        [[r(2,0),r(6,0),r(14,0)], [r(10,0),r(11,0),r(17,3)], [r(24,1),r(23,3),r(25,3)]],
        [[r(0,0),r(1,5),r(9,5)], [r(5,0),r(2,0),r(6,0)], [r(18,1),r(10,0),r(11,0)]],
        [[r(4,1),r(3,5),r(7,5)], [r(8,1),r(0,0),r(1,5)], [r(16,1),r(5,0),r(2,0)]],
    ],
    [
        [[r(7,0),r(21,0),r(38,0)], [r(9,0),r(19,0),r(34,3)], [r(14,1),r(20,3),r(36,3)]],
        [[r(3,0),r(13,5),r(29,5)], [r(1,0),r(7,0),r(21,0)], [r(6,1),r(9,0),r(19,0)]],
        [[r(4,2),r(12,5),r(26,5)], [r(0,1),r(3,0),r(13,5)], [r(2,1),r(1,0),r(7,0)]],
    ],
    [
        [[r(26,0),r(42,0),r(58,0)], [r(29,0),r(43,0),r(62,3)], [r(38,1),r(47,3),r(64,3)]],
        [[r(12,0),r(28,5),r(44,5)], [r(13,0),r(26,0),r(42,0)], [r(21,1),r(29,0),r(43,0)]],
        [[r(4,3),r(15,5),r(31,5)], [r(3,1),r(12,0),r(28,5)], [r(7,1),r(13,0),r(26,0)]],
    ],
    [
        [[r(31,0),r(41,0),r(49,0)], [r(44,0),r(53,0),r(61,3)], [r(58,1),r(65,3),r(75,3)]],
        [[r(15,0),r(22,5),r(33,5)], [r(28,0),r(31,0),r(41,0)], [r(42,1),r(44,0),r(53,0)]],
        [[r(4,4),r(8,5),r(16,5)], [r(12,1),r(15,0),r(22,5)], [r(26,1),r(28,0),r(31,0)]],
    ],
    [
        [[r(50,0),r(48,0),r(49,3)], [r(32,0),r(30,3),r(33,3)], [r(24,3),r(18,3),r(16,3)]],
        [[r(70,0),r(67,0),r(66,3)], [r(52,3),r(50,0),r(48,0)], [r(37,3),r(32,0),r(30,3)]],
        [[r(83,0),r(87,3),r(85,3)], [r(74,3),r(70,0),r(67,0)], [r(57,1),r(52,3),r(50,0)]],
    ],
    [
        [[r(25,0),r(23,0),r(24,3)], [r(17,0),r(11,3),r(10,3)], [r(14,3),r(6,3),r(2,3)]],
        [[r(45,0),r(39,0),r(37,3)], [r(35,3),r(25,0),r(23,0)], [r(27,3),r(17,0),r(11,3)]],
        [[r(63,0),r(59,3),r(57,3)], [r(56,3),r(45,0),r(39,0)], [r(46,3),r(35,3),r(25,0)]],
    ],
    [
        [[r(36,0),r(20,0),r(14,3)], [r(34,0),r(19,3),r(9,3)], [r(38,3),r(21,3),r(7,3)]],
        [[r(55,0),r(40,0),r(27,3)], [r(54,3),r(36,0),r(20,0)], [r(51,3),r(34,0),r(19,3)]],
        [[r(72,0),r(60,3),r(46,3)], [r(73,3),r(55,0),r(40,0)], [r(71,3),r(54,3),r(36,0)]],
    ],
    [
        [[r(64,0),r(47,0),r(38,3)], [r(62,0),r(43,3),r(29,3)], [r(58,3),r(42,3),r(26,3)]],
        [[r(84,0),r(69,0),r(51,3)], [r(82,3),r(64,0),r(47,0)], [r(76,3),r(62,0),r(43,3)]],
        [[r(97,0),r(89,3),r(71,3)], [r(98,3),r(84,0),r(69,0)], [r(96,3),r(82,3),r(64,0)]],
    ],
    [
        [[r(75,0),r(65,0),r(58,3)], [r(61,0),r(53,3),r(44,3)], [r(49,3),r(41,3),r(31,3)]],
        [[r(94,0),r(86,0),r(76,3)], [r(81,3),r(75,0),r(65,0)], [r(66,3),r(61,0),r(53,3)]],
        [[r(107,0),r(104,3),r(96,3)], [r(101,3),r(94,0),r(86,0)], [r(85,3),r(81,3),r(75,0)]],
    ],
    [
        [[r(57,0),r(59,0),r(63,3)], [r(74,0),r(78,3),r(79,3)], [r(83,3),r(92,3),r(95,3)]],
        [[r(37,0),r(39,3),r(45,3)], [r(52,0),r(57,0),r(59,0)], [r(70,3),r(74,0),r(78,3)]],
        [[r(24,0),r(23,3),r(25,3)], [r(32,3),r(37,0),r(39,3)], [r(50,3),r(52,0),r(57,0)]],
    ],
    [
        [[r(46,0),r(60,0),r(72,3)], [r(56,0),r(68,3),r(80,3)], [r(63,3),r(77,3),r(90,3)]],
        [[r(27,0),r(40,3),r(55,3)], [r(35,0),r(46,0),r(60,0)], [r(45,3),r(56,0),r(68,3)]],
        [[r(14,0),r(20,3),r(36,3)], [r(17,3),r(27,0),r(40,3)], [r(25,3),r(35,0),r(46,0)]],
    ],
    [
        [[r(71,0),r(89,0),r(97,3)], [r(73,0),r(91,3),r(103,3)], [r(72,3),r(88,3),r(105,3)]],
        [[r(51,0),r(69,3),r(84,3)], [r(54,0),r(71,0),r(89,0)], [r(55,3),r(73,0),r(91,3)]],
        [[r(38,0),r(47,3),r(64,3)], [r(34,3),r(51,0),r(69,3)], [r(36,3),r(54,0),r(71,0)]],
    ],
    [
        [[r(96,0),r(104,0),r(107,3)], [r(98,0),r(110,3),r(115,3)], [r(97,3),r(111,3),r(119,3)]],
        [[r(76,0),r(86,3),r(94,3)], [r(82,0),r(96,0),r(104,0)], [r(84,3),r(98,0),r(110,3)]],
        [[r(58,0),r(65,3),r(75,3)], [r(62,3),r(76,0),r(86,3)], [r(64,3),r(82,0),r(96,0)]],
    ],
    [
        [[r(85,0),r(87,0),r(83,3)], [r(101,0),r(102,3),r(100,3)], [r(107,3),r(112,3),r(114,3)]],
        [[r(66,0),r(67,3),r(70,3)], [r(81,0),r(85,0),r(87,0)], [r(94,3),r(101,0),r(102,3)]],
        [[r(49,0),r(48,3),r(50,3)], [r(61,3),r(66,0),r(67,3)], [r(75,3),r(81,0),r(85,0)]],
    ],
    [
        [[r(95,0),r(92,0),r(83,0)], [r(79,0),r(78,0),r(74,3)], [r(63,1),r(59,3),r(57,3)]],
        [[r(109,0),r(108,0),r(100,5)], [r(93,1),r(95,0),r(92,0)], [r(77,1),r(79,0),r(78,0)]],
        [[r(117,4),r(118,5),r(114,5)], [r(106,1),r(109,0),r(108,0)], [r(90,1),r(93,1),r(95,0)]],
    ],
    [
        [[r(90,0),r(77,0),r(63,0)], [r(80,0),r(68,0),r(56,3)], [r(72,1),r(60,3),r(46,3)]],
        [[r(106,0),r(93,0),r(79,5)], [r(99,1),r(90,0),r(77,0)], [r(88,1),r(80,0),r(68,0)]],
        [[r(117,3),r(109,5),r(95,5)], [r(113,1),r(106,0),r(93,0)], [r(105,1),r(99,1),r(90,0)]],
    ],
    [
        [[r(105,0),r(88,0),r(72,0)], [r(103,0),r(91,0),r(73,3)], [r(97,1),r(89,3),r(71,3)]],
        [[r(113,0),r(99,0),r(80,5)], [r(116,1),r(105,0),r(88,0)], [r(111,1),r(103,0),r(91,0)]],
        [[r(117,2),r(106,5),r(90,5)], [r(121,1),r(113,0),r(99,0)], [r(119,1),r(116,1),r(105,0)]],
    ],
    [
        [[r(119,0),r(111,0),r(97,0)], [r(115,0),r(110,0),r(98,3)], [r(107,1),r(104,3),r(96,3)]],
        [[r(121,0),r(116,0),r(103,5)], [r(120,1),r(119,0),r(111,0)], [r(112,1),r(115,0),r(110,0)]],
        [[r(117,1),r(113,5),r(105,5)], [r(118,1),r(121,0),r(116,0)], [r(114,1),r(120,1),r(119,0)]],
    ],
    [
        [[r(114,0),r(112,0),r(107,0)], [r(100,0),r(102,0),r(101,3)], [r(83,1),r(87,3),r(85,3)]],
        [[r(118,0),r(120,0),r(115,5)], [r(108,1),r(114,0),r(112,0)], [r(92,1),r(100,0),r(102,0)]],
        [[r(117,0),r(121,5),r(119,5)], [r(109,1),r(118,0),r(120,0)], [r(95,1),r(108,1),r(114,0)]],
    ],
];
