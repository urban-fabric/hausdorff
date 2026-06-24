//! H3 internal coordinate systems: IJK / IJ / cube hex coords, the local 2D hex
//! plane (`Vec2d`), and unit-sphere n-vectors (`Vec3d`) with the icosahedral
//! gnomonic projection. Ported from the H3 reference (h3o); resolution is a
//! plain `usize` index and Class III is `res & 1 == 1`.

use std::ops::{Add, MulAssign, Sub};

mod vec2d;
mod vec3d;

pub(super) use vec2d::Vec2d;
pub(super) use vec3d::Vec3d;

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

/// IJK hex coordinates (each axis 120° apart).
///
/// Proof-backed by `riemann` `Hausdorff/H3Coord` (`CoordIJK` over ℤ): `normalize`
/// is the canonical smallest-non-negative name (`normalize_{nonneg,idem}`,
/// `cmin_normalize`) and is translation-invariant (`normalize_translate`);
/// `distance` is the component-spread **pseudometric** on the hex equivalence
/// (`gridDistance_{eq_spread,self,nonneg,comm,triangle,eq_zero_iff}`); `rotate60`
/// CW/CCW are mutual inverses and compose to the 120° axis 3-cycle
/// (`rotate60_{cw_ccw,ccw_cw}`, `rotate60{ccw,cw}_sq`); `down_aperture{7,3}` are
/// well-defined on the hex equivalence (`downAp*_normalize`); `neighbor` moves
/// one cell (`gridDistance_neighbor_le_one`); and `to_digit` ↔ the digit decode
/// is a bijection on `0..6` (`toDigit_mem`, `unpackDigit_isUnit`,
/// `unpackDigit_toDigit`, `toDigit_unpackDigit`, `neighborByDigit_toDigit`). The
/// `IJ` packing is faithful (`normalize_ofIJ_toIJ`).
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

    /// As [`normalize`](Self::normalize), returning `None` on i32 overflow.
    pub(super) fn checked_normalize(mut self) -> Option<Self> {
        let min = self.i.min(self.j.min(self.k));
        self.i = self.i.checked_sub(min)?;
        self.j = self.j.checked_sub(min)?;
        self.k = self.k.checked_sub(min)?;
        Some(self)
    }

    /// As [`to_digit`](Self::to_digit), returning `None` when the coordinate is
    /// not a unit vector.
    pub(super) fn checked_digit(self) -> Option<u8> {
        let v = self.normalize();
        ((v.i | v.j | v.k) & !1 == 0).then_some(((v.i << 2) | (v.j << 1) | v.k) as u8)
    }

    /// As [`up_aperture7`](Self::up_aperture7), returning `None` on i32
    /// overflow (arbitrary local IJ input can exceed the grid range).
    pub(super) fn checked_up_aperture7<const CCW_: bool>(self) -> Option<Self> {
        let ij = CoordIJ::from(self);
        let (i, j) = if CCW_ {
            (
                f64::from(ij.i.checked_mul(3)?.checked_sub(ij.j)?) * ONE_SEVENTH,
                f64::from(ij.j.checked_mul(2)?.checked_add(ij.i)?) * ONE_SEVENTH,
            )
        } else {
            (
                f64::from(ij.i.checked_mul(2)?.checked_add(ij.j)?) * ONE_SEVENTH,
                f64::from(ij.j.checked_mul(3)?.checked_sub(ij.i)?) * ONE_SEVENTH,
            )
        };
        Self::new(i.round() as i32, j.round() as i32, 0).checked_normalize()
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
