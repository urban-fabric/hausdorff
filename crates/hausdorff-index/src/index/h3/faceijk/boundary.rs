use super::super::coord::{CCW, CW, CoordIJK, SQRT3_2, Vec2d, Vec3d};
use super::super::face::{self, IJ, JK, KI};
use super::{FaceIJK, MAX_DIM_BY_CII_RES, Overage, UNIT_SCALE_BY_CII_RES, is_class3};

impl FaceIJK {
    /// Substrate-grid vertices of this (center) cell. Returns the adjusted
    /// (extended) resolution used for the vertices. `verts` is filled with `n`
    /// entries (6 hexagon, 5 pentagon).
    fn vertices(&mut self, res: usize, n: usize, verts: &mut [FaceIJK; 6]) -> usize {
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
    pub(in crate::index::h3) fn cell_boundary_lat_lng_rad(
        bits: u64,
        base_cell: u8,
        res: usize,
        is_pentagon: bool,
    ) -> Vec<(f64, f64)> {
        Self::from_cell(bits, base_cell, res).boundary_lat_lng_rad(res, is_pentagon)
    }

    pub(in crate::index::h3) fn cell_boundary_segment_lat_lng_rad(
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
}
