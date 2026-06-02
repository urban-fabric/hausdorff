//! H3 internal coordinate systems: IJK / IJ / cube hex coords, the local 2D hex
//! plane (`Vec2d`), and unit-sphere n-vectors (`Vec3d`) with the icosahedral
//! gnomonic projection. Ported from the H3 reference (h3o); resolution is a
//! plain `usize` index and Class III is `res & 1 == 1`.

use std::ops::{Add, MulAssign, Sub};

use super::face;

/// CCW rotation marker for the const-generic rotation helpers.
pub(super) const CCW: bool = true;
/// CW rotation marker for the const-generic rotation helpers.
pub(super) const CW: bool = false;

const EPSILON: f64 = 0.000_000_000_000_000_1;
const RES0_U_GNOMONIC: f64 = 0.381_966_011_250_105;
const INV_RES0_U_GNOMONIC: f64 = 2.618_033_988_749_896;
const AP7_ROT_RADS: f64 = 0.333_473_172_251_832_1;
pub(super) const SQRT3_2: f64 = 0.866_025_403_784_438_6;
const ONE_SEVENTH: f64 = 1.0 / 7.0;
const ONE_THIRD: f64 = 1.0 / 3.0;
const RSIN60: f64 = 1.154_700_538_379_251_5;
const TWO_PI: f64 = std::f64::consts::TAU;

#[rustfmt::skip]
const SQRT7_POWERS: [f64; 17] = [
    1.0, 2.6457513110645907, 7., 18.520259177452136, 49.00000000000001,
    129.64181424216497, 343.0000000000001, 907.4926996951549, 2401.000000000001,
    6352.448897866085, 16807.000000000007, 44467.1422850626, 117649.00000000007,
    311269.9959954382, 823543.0000000006, 2178889.971968068, 5764801.0,
];

#[rustfmt::skip]
const INV_SQRT7_POWERS: [f64; 17] = [
    1., 0.3779644730092272, 0.14285714285714285, 0.05399492471560388,
    0.02040816326530612, 0.007713560673657697, 0.002915451895043731,
    0.0011019372390939565, 0.0004164931278633901, 0.00015741960558485093,
    0.00005949901826619858, 0.00002248851508355013, 0.000008499859752314082,
    0.000003212645011935733, 0.0000012142656789020115, 0.0000004589492874193903,
    0.00000017346652555743034,
];

#[inline]
const fn is_class3(res: usize) -> bool {
    res & 1 == 1
}

fn to_positive_angle(mut angle: f64) -> f64 {
    if angle < 0.0 {
        angle += TWO_PI;
    } else if angle >= TWO_PI {
        angle -= TWO_PI;
    }
    angle
}

// ---- IJ ---------------------------------------------------------------------

/// IJ hex coordinates (each axis 120° apart).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct CoordIJ {
    pub(super) i: i32,
    pub(super) j: i32,
}

impl CoordIJ {
    pub(super) const fn new(i: i32, j: i32) -> Self {
        Self { i, j }
    }
}

// ---- IJK --------------------------------------------------------------------

/// IJK hex coordinates (each axis 120° apart).
#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub(super) struct CoordIJK {
    i: i32,
    j: i32,
    k: i32,
}

impl CoordIJK {
    #[inline]
    pub(super) const fn new(i: i32, j: i32, k: i32) -> Self {
        Self { i, j, k }
    }

    pub(super) const fn i(self) -> i32 {
        self.i
    }
    pub(super) const fn j(self) -> i32 {
        self.j
    }
    pub(super) const fn k(self) -> i32 {
        self.k
    }

    /// Normalizes to the smallest non-negative components (at most two nonzero).
    #[inline]
    pub(super) fn normalize(mut self) -> Self {
        let min = self.i.min(self.j.min(self.k));
        self.i -= min;
        self.j -= min;
        self.k -= min;
        self
    }

    /// Grid distance between two IJK coordinates.
    pub(super) fn distance(self, other: Self) -> i32 {
        let d = (self - other).normalize();
        d.i.abs().max(d.j.abs().max(d.k.abs()))
    }

    /// Unit vector → direction digit (0..=6). Caller guarantees a unit vector.
    pub(super) fn to_digit(self) -> u8 {
        let v = self.normalize();
        ((v.i << 2) | (v.j << 1) | v.k) as u8
    }

    /// Indexing parent in an aperture-7 grid.
    #[inline]
    pub(super) fn up_aperture7<const CCW_: bool>(self) -> Self {
        let ij = CoordIJ::from(self);
        let (i, j) = if CCW_ {
            (
                f64::from(3 * ij.i - ij.j) * ONE_SEVENTH,
                f64::from(ij.i + 2 * ij.j) * ONE_SEVENTH,
            )
        } else {
            (
                f64::from(2 * ij.i + ij.j) * ONE_SEVENTH,
                f64::from(3 * ij.j - ij.i) * ONE_SEVENTH,
            )
        };
        Self::new(i.round() as i32, j.round() as i32, 0).normalize()
    }

    /// Center of this hex at the next finer aperture-7 resolution.
    #[inline]
    pub(super) fn down_aperture7<const CCW_: bool>(self) -> Self {
        let (mut iv, mut jv, mut kv) = if CCW_ {
            (Self::new(3, 0, 1), Self::new(1, 3, 0), Self::new(0, 1, 3))
        } else {
            (Self::new(3, 1, 0), Self::new(0, 3, 1), Self::new(1, 0, 3))
        };
        iv *= self.i;
        jv *= self.j;
        kv *= self.k;
        (iv + jv + kv).normalize()
    }

    /// Center of this hex at the next finer aperture-3 resolution.
    #[inline]
    pub(super) fn down_aperture3<const CCW_: bool>(self) -> Self {
        let (mut iv, mut jv, mut kv) = if CCW_ {
            (Self::new(2, 0, 1), Self::new(1, 2, 0), Self::new(0, 1, 2))
        } else {
            (Self::new(2, 1, 0), Self::new(0, 2, 1), Self::new(1, 0, 2))
        };
        iv *= self.i;
        jv *= self.j;
        kv *= self.k;
        (iv + jv + kv).normalize()
    }

    /// Neighbor in the given direction digit (0..=6).
    #[inline]
    pub(super) fn neighbor(self, digit: u8) -> Self {
        let c = Self::new(
            i32::from((digit >> 2) & 1),
            i32::from((digit >> 1) & 1),
            i32::from(digit & 1),
        );
        (self + c).normalize()
    }

    /// 60° rotation.
    #[inline]
    pub(super) fn rotate60<const CCW_: bool>(self) -> Self {
        let (mut iv, mut jv, mut kv) = if CCW_ {
            (Self::new(1, 1, 0), Self::new(0, 1, 1), Self::new(1, 0, 1))
        } else {
            (Self::new(1, 0, 1), Self::new(1, 1, 0), Self::new(0, 1, 1))
        };
        iv *= self.i;
        jv *= self.j;
        kv *= self.k;
        (iv + jv + kv).normalize()
    }
}

impl Add for CoordIJK {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.i + o.i, self.j + o.j, self.k + o.k)
    }
}
impl Sub for CoordIJK {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.i - o.i, self.j - o.j, self.k - o.k)
    }
}
impl MulAssign<i32> for CoordIJK {
    fn mul_assign(&mut self, rhs: i32) {
        self.i *= rhs;
        self.j *= rhs;
        self.k *= rhs;
    }
}

impl From<CoordIJK> for CoordIJ {
    #[inline]
    fn from(v: CoordIJK) -> Self {
        Self::new(v.i - v.k, v.j - v.k)
    }
}

impl From<CoordIJK> for Vec2d {
    // Center point in 2D cartesian coordinates of a hex.
    #[inline]
    fn from(v: CoordIJK) -> Self {
        let i = f64::from(v.i - v.k);
        let j = f64::from(v.j - v.k);
        Self::new(0.5_f64.mul_add(-j, i), j * SQRT3_2)
    }
}

// ---- Vec2d (local hex plane) ------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub(super) struct Vec2d {
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
    pub(super) const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub(super) fn magnitude(self) -> f64 {
        self.x.hypot(self.y)
    }

    /// Encodes an n-vector to icosahedral 2D hex coordinates relative to the
    /// given face center. `distance` is the squared euclidean distance to it.
    pub(super) fn from_vec3d(value: Vec3d, res: usize, face: usize, distance: f64) -> Self {
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

    pub(super) fn intersection(line1: (Self, Self), line2: (Self, Self)) -> Self {
        let s1 = Self::new(line1.1.x - line1.0.x, line1.1.y - line1.0.y);
        let s2 = Self::new(line2.1.x - line2.0.x, line2.1.y - line2.0.y);
        let t =
            s2.x.mul_add(line1.0.y - line2.0.y, -s2.y * (line1.0.x - line2.0.x))
                / (-s2.x).mul_add(s1.y, s1.x * s2.y);
        Self::new(t.mul_add(s1.x, line1.0.x), t.mul_add(s1.y, line1.0.y))
    }
}

impl From<Vec2d> for CoordIJK {
    // Containing hex in IJK for a 2D cartesian coordinate vector (DGGRID).
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

// ---- Vec3d (unit-sphere n-vector) -------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Vec3d {
    pub(super) x: f64,
    pub(super) y: f64,
    pub(super) z: f64,
}

const NORTH_POLE: Vec3d = Vec3d {
    x: 0.0,
    y: 0.0,
    z: 1.0,
};

impl Vec3d {
    pub(super) const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub(super) fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    /// n-vector of a cell center, given its 2D hex coords on `face`.
    /// `res` is the (extended) resolution index; `is_substrate` selects the
    /// substrate grid used for boundary vertices.
    pub(super) fn from_vec2d(value: Vec2d, face: usize, res: usize, is_substrate: bool) -> Self {
        let r = {
            let mut r = value.magnitude();
            if r < EPSILON {
                return face::CENTER_POINT[face];
            }
            r *= INV_SQRT7_POWERS[res];
            if is_substrate {
                r *= ONE_THIRD;
            }
            (r * RES0_U_GNOMONIC).atan()
        };
        let theta = {
            let mut theta = value.y.atan2(value.x);
            if !is_substrate && is_class3(res) {
                theta = to_positive_angle(theta + AP7_ROT_RADS);
            }
            to_positive_angle(face::AXES_AZ_RADS_CII[face] - theta)
        };
        let center = face::CENTER_POINT[face];
        if r < EPSILON {
            return center;
        }
        let (north, east) = center.tangent_basis();
        let dir = linear_combination(theta.cos(), &north, theta.sin(), &east);
        let mut res_v = linear_combination(r.cos(), &center, r.sin(), &dir);
        res_v.normalize();
        res_v
    }

    /// Closest icosahedron face and the squared euclidean distance to it.
    pub(super) fn closest_face(self) -> (usize, f64) {
        const MAX_DIST: f64 = 5.0;
        let mut best = (0usize, MAX_DIST);
        for (i, center) in face::CENTER_POINT.iter().enumerate() {
            let d = self.distance_squared(center);
            if d < best.1 {
                best = (i, d);
            }
        }
        best
    }

    /// Azimuth (radians) from `self` to `other`.
    pub(super) fn azimuth(&self, other: &Self) -> f64 {
        let (north, east) = self.tangent_basis();
        let mut proj = linear_combination(1.0, other, -other.dot(self), self);
        proj.normalize();
        proj.dot(&east).atan2(proj.dot(&north))
    }

    fn tangent_basis(&self) -> (Self, Self) {
        let mut north = linear_combination(1.0, &NORTH_POLE, -NORTH_POLE.dot(self), self);
        north.normalize();
        let east = north.cross(self);
        (north, east)
    }

    fn distance_squared(&self, other: &Self) -> f64 {
        linear_combination(1.0, self, -1.0, other).norm_squared()
    }

    fn normalize(&mut self) {
        let norm = self.norm();
        if norm > 0.0 {
            let s = 1.0 / norm;
            self.x *= s;
            self.y *= s;
            self.z *= s;
        } else {
            self.x = 0.0;
            self.y = 0.0;
            self.z = 0.0;
        }
    }

    fn cross(&self, o: &Self) -> Self {
        Self::new(
            self.y.mul_add(o.z, -(self.z * o.y)),
            self.z.mul_add(o.x, -(self.x * o.z)),
            self.x.mul_add(o.y, -(self.y * o.x)),
        )
    }

    fn dot(&self, o: &Self) -> f64 {
        self.x.mul_add(o.x, self.y.mul_add(o.y, self.z * o.z))
    }

    fn norm_squared(&self) -> f64 {
        self.dot(self)
    }
}

#[inline]
fn linear_combination(a: f64, v1: &Vec3d, b: f64, v2: &Vec3d) -> Vec3d {
    Vec3d::new(
        a.mul_add(v1.x, b * v2.x),
        a.mul_add(v1.y, b * v2.y),
        a.mul_add(v1.z, b * v2.z),
    )
}

// ---- LatLng <-> Vec3d -------------------------------------------------------

impl Vec3d {
    /// n-vector from geographic coordinates in radians.
    pub(super) fn from_lat_lng_rad(lat: f64, lng: f64) -> Self {
        let r = lat.cos();
        Self::new(lng.cos() * r, lng.sin() * r, lat.sin())
    }

    /// Geographic coordinates (radians) of this n-vector: `(lat, lng)`.
    /// Matches h3o's `LatLng::from(Vec3d)` exactly: `asin(z)` / `atan2(y, x)`.
    pub(super) fn to_lat_lng_rad(self) -> (f64, f64) {
        (self.z.asin(), self.y.atan2(self.x))
    }
}
