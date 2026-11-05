//! `S2CellId` — a 64-bit S2 cell identifier (golang/geo's `cellid.go`).
//!
//! The most-significant 3 bits encode the cube face (0..5); the remaining 61
//! bits encode the cell's position along the face's Hilbert curve, which
//! implicitly carries both the cell's location and its subdivision level (the
//! level is fixed by the position of the lowest set bit). Sequentially
//! increasing ids trace a continuous space-filling curve over the sphere.

use serde::{Deserialize, Serialize};

use super::coords::{self, MAX_LEVEL, MAX_SIZE, NUM_FACES};
use crate::space::s2::S2Point;
use hausdorff_base::earth::GeoPoint;

mod lookup;
mod neighbors;
mod token;

use lookup::{INVERT_MASK, LOOKUP, LOOKUP_BITS, SWAP_MASK};

/// Bits used to encode the face number.
const FACE_BITS: u32 = 3;
/// Total number of position bits (`2*MAX_LEVEL + 1`).
const POS_BITS: u32 = 2 * MAX_LEVEL + 1;

/// A 64-bit S2 cell identifier. The zero value and `u64::MAX` are invalid.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct S2CellId(u64);

impl S2CellId {
    /// Wraps a raw 64-bit value. Use [`S2CellId::is_valid`] to check validity.
    #[must_use]
    pub const fn from_u64(value: u64) -> Self {
        Self(value)
    }

    /// The raw 64-bit value.
    #[must_use]
    pub const fn as_u64(self) -> u64 {
        self.0
    }

    /// The face cell (level 0) for a given cube face in `[0, 5]`.
    #[must_use]
    pub fn from_face(face: u8) -> Self {
        Self((u64::from(face) << POS_BITS) + lsb_for_level(0))
    }

    /// A cell from a face, 61-bit Hilbert position, and level; the position is
    /// truncated to the cell center at `level`.
    #[must_use]
    pub fn from_face_pos_level(face: u8, pos: u64, level: u32) -> Self {
        Self(((u64::from(face) << POS_BITS) + pos) | 1).parent(level)
    }

    /// The leaf cell containing the given WGS84 coordinate.
    #[must_use]
    pub fn from_geo_point(p: GeoPoint) -> Self {
        Self::from_point(S2Point::from_geo_point(p))
    }

    /// The cell at `level` containing the given WGS84 coordinate.
    #[must_use]
    pub fn from_geo_point_at_level(p: GeoPoint, level: u32) -> Self {
        Self::from_geo_point(p).parent(level)
    }

    /// The leaf cell containing the given unit-sphere point.
    #[must_use]
    pub fn from_point(p: S2Point) -> Self {
        let (f, u, v) = coords::xyz_to_face_uv(p.vector());
        let i = coords::st_to_ij(coords::uv_to_st(u));
        let j = coords::st_to_ij(coords::uv_to_st(v));
        Self::from_face_ij(i64::from(f), i, j)
    }

    /// The cube face (0..5).
    #[must_use]
    pub const fn face(self) -> u8 {
        (self.0 >> POS_BITS) as u8
    }

    /// The Hilbert curve position within the face, in `[0, 2^PosBits - 1]`.
    #[must_use]
    pub const fn pos(self) -> u64 {
        self.0 & (u64::MAX >> FACE_BITS)
    }

    /// The subdivision level in `[0, MAX_LEVEL]`. Only meaningful for a valid
    /// cell.
    #[must_use]
    pub const fn level(self) -> u32 {
        MAX_LEVEL - (self.0.trailing_zeros() >> 1)
    }

    /// Whether this is a leaf cell (deepest level).
    #[must_use]
    pub const fn is_leaf(self) -> bool {
        self.0 & 1 != 0
    }

    /// Whether this is a top-level (face) cell.
    #[must_use]
    pub const fn is_face(self) -> bool {
        self.0 & (lsb_for_level(0) - 1) == 0
    }

    /// Whether this id represents a valid cell.
    #[must_use]
    pub const fn is_valid(self) -> bool {
        (self.face() as u64) < NUM_FACES as u64 && (self.lsb() & 0x1555_5555_5555_5555 != 0)
    }

    const fn lsb(self) -> u64 {
        self.0 & self.0.wrapping_neg()
    }

    /// The child position (0..3) of this cell's ancestor at `level`, relative
    /// to that ancestor's parent. `level` must be in `1..=MAX_LEVEL`.
    #[must_use]
    pub const fn child_position(self, level: u32) -> u8 {
        ((self.0 >> (2 * (MAX_LEVEL - level) + 1)) & 3) as u8
    }

    /// The ancestor cell at `level`, which must be no greater than this cell's
    /// level.
    #[must_use]
    pub const fn parent(self, level: u32) -> Self {
        let lsb = lsb_for_level(level);
        Self((self.0 & lsb.wrapping_neg()) | lsb)
    }

    /// The immediate parent. Assumes this is not a face cell.
    #[must_use]
    pub const fn immediate_parent(self) -> Self {
        let nlsb = self.lsb() << 2;
        Self((self.0 & nlsb.wrapping_neg()) | nlsb)
    }

    /// The four immediate children in Hilbert curve order. For a leaf cell the
    /// result is not meaningful.
    #[must_use]
    pub fn children(self) -> [Self; 4] {
        let mut lsb = self.lsb();
        let c0 = self.0 - lsb + (lsb >> 2);
        lsb >>= 1;
        let c1 = c0 + lsb;
        let c2 = c1 + lsb;
        let c3 = c2 + lsb;
        [Self(c0), Self(c1), Self(c2), Self(c3)]
    }

    /// The first child in a Hilbert-order traversal of this cell's children.
    #[must_use]
    pub const fn child_begin(self) -> Self {
        let ol = self.lsb();
        Self(self.0 - ol + (ol >> 2))
    }

    /// One past the last child in a Hilbert-order traversal (may be invalid).
    #[must_use]
    pub const fn child_end(self) -> Self {
        let ol = self.lsb();
        Self(self.0 + ol + (ol >> 2))
    }

    /// The next cell along the Hilbert curve at this cell's level.
    #[must_use]
    pub const fn next(self) -> Self {
        Self(self.0.wrapping_add(self.lsb() << 1))
    }

    /// The previous cell along the Hilbert curve at this cell's level.
    #[must_use]
    pub const fn prev(self) -> Self {
        Self(self.0.wrapping_sub(self.lsb() << 1))
    }

    /// The minimum (leaf) cell id contained within this cell.
    #[must_use]
    pub const fn range_min(self) -> Self {
        Self(self.0 - (self.lsb() - 1))
    }

    /// The maximum (leaf) cell id contained within this cell.
    #[must_use]
    pub const fn range_max(self) -> Self {
        Self(self.0 + (self.lsb() - 1))
    }

    /// Whether this cell contains `other` (i.e. `other` is this cell or a
    /// descendant).
    ///
    /// Proof-backed by `riemann` `Hausdorff/S2CellOrder`: with the half-open
    /// id-range `[range_min, range_max]`, `contains` is reflexive
    /// (`contains_refl`) and implies [`intersects`](Self::intersects)
    /// (`contains_imp_intersects`); true range-containment is a partial order
    /// (`rangeSubset_{refl,trans,antisymm}`) of which this id-test is the sound
    /// image (`contains_of_rangeSubset`).
    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        self.range_min().0 <= other.0 && other.0 <= self.range_max().0
    }

    /// Whether this cell's subtree intersects `other`'s subtree.
    ///
    /// Proof-backed by `riemann` `Hausdorff/S2CellOrder`: range-overlap is
    /// reflexive and *symmetric* (`intersects_refl`, `intersects_comm`), and the
    /// negative form `not_intersects_of_hi_lt_lo` certifies the pruning — disjoint
    /// id-ranges share no cell, so the query can skip this subtree.
    #[must_use]
    pub const fn intersects(self, other: Self) -> bool {
        other.range_min().0 <= self.range_max().0 && other.range_max().0 >= self.range_min().0
    }

    /// The level of the lowest common ancestor of the two cells, or `None` if
    /// they share no common ancestor (different faces).
    #[must_use]
    pub fn common_ancestor_level(self, other: Self) -> Option<u32> {
        let diff = (self.0 ^ other.0).max(self.lsb()).max(other.lsb());
        let msb_pos = 63 - diff.leading_zeros();
        if msb_pos > 60 {
            return None;
        }
        Some((60 - msb_pos) >> 1)
    }

    /// The cell center as a unit-sphere [`S2Point`].
    #[must_use]
    pub fn center_point(self) -> S2Point {
        S2Point::from_vector(self.raw_point().normalize())
    }

    /// The cell center as a WGS84 coordinate.
    #[must_use]
    pub fn center(self) -> GeoPoint {
        S2Point::from_vector(self.raw_point()).to_geo_point()
    }

    fn raw_point(self) -> crate::space::s2::vector::R3Vector {
        let (f, si, ti) = self.face_si_ti();
        coords::face_uv_to_xyz(
            f as u8,
            coords::st_to_uv((0.5 / MAX_SIZE as f64) * f64::from(si)),
            coords::st_to_uv((0.5 / MAX_SIZE as f64) * f64::from(ti)),
        )
    }

    fn face_si_ti(self) -> (i64, u32, u32) {
        let (f, i, j, _) = self.face_ij_orientation();
        let delta = if self.is_leaf() {
            1
        } else if (i ^ ((self.0 >> 2) as i64)) & 1 != 0 {
            2
        } else {
            0
        };
        (f, (2 * i + delta) as u32, (2 * j + delta) as u32)
    }

    /// The cell's four corner vertices as unit-sphere [`S2Point`]s, ordered
    /// counter-clockwise in the face's `(u, v)` frame (`(i,j)`, `(i+s,j)`,
    /// `(i+s,j+s)`, `(i,j+s)`). The corners are the cell's `(i, j)` lattice
    /// square at this level, projected through the quadratic cube→sphere
    /// transform — the geometry a region coverer tests cells against.
    #[must_use]
    pub fn vertices(self) -> [S2Point; 4] {
        let (f, i, j, _) = self.face_ij_orientation();
        let size = size_ij(self.level());
        let (i, j) = (i & -size, j & -size);
        let corner = |ci: i64, cj: i64| {
            let u = coords::st_to_uv(ci as f64 / MAX_SIZE as f64);
            let v = coords::st_to_uv(cj as f64 / MAX_SIZE as f64);
            S2Point::from_vector(coords::face_uv_to_xyz(f as u8, u, v).normalize())
        };
        [
            corner(i, j),
            corner(i + size, j),
            corner(i + size, j + size),
            corner(i, j + size),
        ]
    }

    /// The cell's `(u, v)` rectangle on its cube face, as
    /// `(face, (u_lo, u_hi), (v_lo, v_hi))` — the exact footprint the lat/lng
    /// bound of a cell is computed from.
    /// `pub(crate)`: `shape::rect` builds its cell bounds from this.
    pub(crate) fn face_uv_bounds(self) -> (u8, (f64, f64), (f64, f64)) {
        let (f, i, j, _) = self.face_ij_orientation();
        let size = size_ij(self.level());
        let (i, j) = (i & -size, j & -size);
        let bound = |lo: i64| {
            (
                coords::st_to_uv(lo as f64 / MAX_SIZE as f64),
                coords::st_to_uv((lo + size) as f64 / MAX_SIZE as f64),
            )
        };
        (f as u8, bound(i), bound(j))
    }

    fn face_ij_orientation(self) -> (i64, i64, i64, i64) {
        let lookup = &*LOOKUP;
        let f = i64::from(self.face());
        let mut orientation = f & SWAP_MASK;
        let mut i: i64 = 0;
        let mut j: i64 = 0;
        let mut nbits = MAX_LEVEL as i64 - 7 * LOOKUP_BITS as i64;

        for k in (0..=7i64).rev() {
            let shift = (k as u32) * 2 * LOOKUP_BITS + 1;
            let chunk = ((self.0 >> shift) & ((1u64 << (2 * nbits)) - 1)) as i64;
            orientation += chunk << 2;
            orientation = lookup.ij[orientation as usize];
            i += (orientation >> (LOOKUP_BITS + 2)) << (k * LOOKUP_BITS as i64);
            j += ((orientation >> 2) & ((1 << LOOKUP_BITS) - 1)) << (k * LOOKUP_BITS as i64);
            orientation &= SWAP_MASK | INVERT_MASK;
            nbits = LOOKUP_BITS as i64;
        }

        if self.lsb() & 0x1111_1111_1111_1110 != 0 {
            orientation ^= SWAP_MASK;
        }
        (f, i, j, orientation)
    }

    fn from_face_ij(f: i64, i: i64, j: i64) -> Self {
        let lookup = &*LOOKUP;
        let mut n: u64 = (f as u64) << (POS_BITS - 1);
        let mut bits = f & SWAP_MASK;
        for k in (0..=7i64).rev() {
            let mask = (1 << LOOKUP_BITS) - 1;
            bits += ((i >> (k * LOOKUP_BITS as i64)) & mask) << (LOOKUP_BITS + 2);
            bits += ((j >> (k * LOOKUP_BITS as i64)) & mask) << 2;
            bits = lookup.pos[bits as usize];
            n |= ((bits >> 2) as u64) << ((k as u64) * 2 * u64::from(LOOKUP_BITS));
            bits &= SWAP_MASK | INVERT_MASK;
        }
        Self(n * 2 + 1)
    }

    fn from_face_ij_wrap(f: i64, i: i64, j: i64) -> Self {
        let i = coords::clamp_i64(i, -1, MAX_SIZE);
        let j = coords::clamp_i64(j, -1, MAX_SIZE);

        let scale = 1.0 / MAX_SIZE as f64;
        let limit = f64::from_bits(1.0_f64.to_bits() + 1);
        let u = (scale * ((i << 1) + 1 - MAX_SIZE) as f64)
            .max(-limit)
            .min(limit);
        let v = (scale * ((j << 1) + 1 - MAX_SIZE) as f64)
            .max(-limit)
            .min(limit);

        let (f, u, v) = coords::xyz_to_face_uv(coords::face_uv_to_xyz(f as u8, u, v));
        Self::from_face_ij(
            i64::from(f),
            coords::st_to_ij(0.5 * (u + 1.0)),
            coords::st_to_ij(0.5 * (v + 1.0)),
        )
    }

    fn from_face_ij_same(f: i64, i: i64, j: i64, same_face: bool) -> Self {
        if same_face {
            Self::from_face_ij(f, i, j)
        } else {
            Self::from_face_ij_wrap(f, i, j)
        }
    }
}

/// The lowest-numbered set bit for a cell at `level`.
const fn lsb_for_level(level: u32) -> u64 {
    1u64 << (2 * (MAX_LEVEL - level))
}

/// Edge length, in leaf-cell `(i,j)` units, of a cell at `level`.
const fn size_ij(level: u32) -> i64 {
    1 << (MAX_LEVEL - level)
}
