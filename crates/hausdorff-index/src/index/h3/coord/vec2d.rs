use super::super::face;
use super::{
    AP7_ROT_RADS, CoordIJK, EPSILON, INV_RES0_U_GNOMONIC, RSIN60, SQRT3_2, SQRT7_POWERS, Vec3d,
    is_class3,
};

impl From<CoordIJK> for Vec2d {
    #[inline]
    fn from(v: CoordIJK) -> Self {
        let i = f64::from(v.i - v.k);
        let j = f64::from(v.j - v.k);
        Self::new(0.5_f64.mul_add(-j, i), j * SQRT3_2)
    }
}

#[derive(Debug, Clone, Copy)]
pub(in crate::index::h3) struct Vec2d {
    pub(super) x: f64,
    pub(super) y: f64,
}

impl PartialEq for Vec2d {
    fn eq(&self, other: &Self) -> bool {
        let e = f64::from(f32::EPSILON);
        (self.x - other.x).abs() <= e && (self.y - other.y).abs() <= e
    }
}

impl Vec2d {
    pub(in crate::index::h3) const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub(super) fn magnitude(self) -> f64 {
        self.x.hypot(self.y)
    }

    /// Encodes an n-vector to icosahedral 2D hex coordinates relative to the
    /// given face center. `distance` is the squared euclidean distance to it.
    pub(in crate::index::h3) fn from_vec3d(
        value: Vec3d,
        res: usize,
        face: usize,
        distance: f64,
    ) -> Self {
        let r = {
            let r = distance.mul_add(-0.5, 1.0).acos();
            if r < EPSILON {
                return Self::new(0.0, 0.0);
            }
            (r.tan() * INV_RES0_U_GNOMONIC) * SQRT7_POWERS[res]
        };
        let theta = {
            let mut theta = face::AXES_AZ_RADS_CII[face] - face::CENTER_POINT[face].azimuth(&value);
            if is_class3(res) {
                theta -= AP7_ROT_RADS;
            }
            theta
        };
        Self::new(r * theta.cos(), r * theta.sin())
    }

    /// As [`from_vec3d`](Self::from_vec3d) but computing the gnomonic radial
    /// factor via the identity `tan(acos(x)) = sqrt(1 - x²)/x`, trading the exact
    /// path's `acos` + `tan` for one `sqrt` and a divide. It is mathematically
    /// equal; only the last ULPs can differ (`1 - x²` is evaluated as the
    /// cancellation-free `distance·(1 - distance/4)`, since `x ≈ 1` in the common
    /// small-angle case), so a cell assignment can shift to an adjacent — still
    /// valid — cell only for a point essentially on a boundary.
    pub(in crate::index::h3) fn from_vec3d_fast(
        value: Vec3d,
        res: usize,
        face: usize,
        distance: f64,
    ) -> Self {
        let r = {
            let x = distance.mul_add(-0.5, 1.0);
            let one_minus_x2 = distance * 0.25_f64.mul_add(-distance, 1.0);
            let tan_r = one_minus_x2.sqrt() / x;
            if tan_r < EPSILON {
                return Self::new(0.0, 0.0);
            }
            (tan_r * INV_RES0_U_GNOMONIC) * SQRT7_POWERS[res]
        };
        let theta = {
            let mut theta = face::AXES_AZ_RADS_CII[face] - face::CENTER_POINT[face].azimuth(&value);
            if is_class3(res) {
                theta -= AP7_ROT_RADS;
            }
            theta
        };
        Self::new(r * theta.cos(), r * theta.sin())
    }

    pub(in crate::index::h3) fn intersection(line1: (Self, Self), line2: (Self, Self)) -> Self {
        let s1 = Self::new(line1.1.x - line1.0.x, line1.1.y - line1.0.y);
        let s2 = Self::new(line2.1.x - line2.0.x, line2.1.y - line2.0.y);
        let t =
            s2.x.mul_add(line1.0.y - line2.0.y, -s2.y * (line1.0.x - line2.0.x))
                / (-s2.x).mul_add(s1.y, s1.x * s2.y);
        Self::new(t.mul_add(s1.x, line1.0.x), t.mul_add(s1.y, line1.0.y))
    }
}

impl From<Vec2d> for CoordIJK {
    #[inline]
    fn from(value: Vec2d) -> Self {
        let k = 0;
        let a1 = value.x.abs();
        let a2 = value.y.abs();

        let x2 = a2 * RSIN60;
        let x1 = a1 + x2 / 2.0;

        let m1 = x1 as i32;
        let m2 = x2 as i32;

        let r1 = x1 - f64::from(m1);
        let r2 = x2 - f64::from(m2);

        let (mut i, mut j) = if r1 < 0.5 {
            if r1 < 1.0 / 3.0 {
                let i = m1;
                let j = m2 + i32::from(r2 >= f64::midpoint(1.0, r1));
                (i, j)
            } else {
                let i = m1 + i32::from((1.0 - r1) <= r2 && r2 < (2.0 * r1));
                let j = m2 + i32::from(r2 >= (1.0 - r1));
                (i, j)
            }
        } else if r1 < 2.0 / 3.0 {
            let j = m2 + i32::from(r2 >= (1.0 - r1));
            let i = m1 + i32::from(2.0_f64.mul_add(r1, -1.0) >= r2 || r2 >= (1.0 - r1));
            (i, j)
        } else {
            let i = m1 + 1;
            let j = m2 + i32::from(r2 >= (r1 / 2.0));
            (i, j)
        };

        if value.x < 0.0 {
            let offset = j % 2;
            let axis_i = i32::midpoint(j, offset);
            let diff = i - axis_i;
            i -= 2 * diff + offset;
        }
        if value.y < 0.0 {
            i -= (2 * j + 1) / 2;
            j = -j;
        }

        Self::new(i, j, k).normalize()
    }
}
