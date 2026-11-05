//! S2 cube-sphere coordinate transforms (golang/geo's `stuv.go`).
//!
//! The chain converting a lat/lng to a cell id walks several coordinate
//! systems: direction vector `(x,y,z)` → cube face + `(u,v)` → cell-space
//! `(s,t)` → leaf-cell indices `(i,j)`. S2 offers linear/tangent/quadratic
//! projections between `(s,t)` and `(u,v)`; like golang/geo we use the
//! **quadratic** transform.

use crate::space::s2::vector::{Axis, R3Vector};

/// Number of subdivision levels needed to specify a leaf cell.
pub(super) const MAX_LEVEL: u32 = 30;

/// Number of cube faces.
pub(super) const NUM_FACES: u8 = 6;

/// Number of leaf cells along an edge of a face (`2^MAX_LEVEL`); also one past
/// the maximum valid leaf `i`/`j` index.
pub(super) const MAX_SIZE: i64 = 1 << MAX_LEVEL;

/// Clamp `v` to `[lo, hi]` (integer analogue of S2's `clamp`).
pub(super) fn clamp_i64(v: i64, lo: i64, hi: i64) -> i64 {
    v.max(lo).min(hi)
}

/// The quadratic `s`/`t` → `u`/`v` transform (makes projected cells more
/// uniform in size).
///
/// Proof-backed by `riemann` `Hausdorff/S2Coords.uvToST_stToUV`: `uv_to_st` is
/// the **exact** left inverse of `st_to_uv` on `[0, 1]` (each branch's radicand
/// collapses to a perfect square), so any float non-idempotence is pure rounding;
/// with the corner anchors `stToUV_{zero,half,one}` (`-1`/`0`/`1`).
pub(super) fn st_to_uv(s: f64) -> f64 {
    if s >= 0.5 {
        (1.0 / 3.0) * (4.0 * s * s - 1.0)
    } else {
        (1.0 / 3.0) * (1.0 - 4.0 * (1.0 - s) * (1.0 - s))
    }
}

/// Inverse of [`st_to_uv`]. Not exactly idempotent due to floating-point error.
pub(super) fn uv_to_st(u: f64) -> f64 {
    if u >= 0.0 {
        0.5 * (1.0 + 3.0 * u).sqrt()
    } else {
        1.0 - 0.5 * (1.0 - 3.0 * u).sqrt()
    }
}

/// Converts an `s`/`t` value to a leaf-cell `i`/`j` index in `[0, MAX_SIZE-1]`.
pub(super) fn st_to_ij(s: f64) -> i64 {
    clamp_i64((MAX_SIZE as f64 * s).floor() as i64, 0, MAX_SIZE - 1)
}

/// The cube face (0..5) containing direction `r`. Boundary points resolve
/// deterministically.
pub(super) fn face(r: R3Vector) -> u8 {
    match r.largest_component() {
        Axis::X => {
            if r.x < 0.0 {
                3
            } else {
                0
            }
        }
        Axis::Y => {
            if r.y < 0.0 {
                4
            } else {
                1
            }
        }
        Axis::Z => {
            if r.z < 0.0 {
                5
            } else {
                2
            }
        }
    }
}

/// Given a face for which `r` projects validly (positive dot product with the
/// face normal), returns the `(u, v)` coordinates (which may fall outside
/// `[-1,1]`).
pub(super) fn valid_face_xyz_to_uv(face: u8, r: R3Vector) -> (f64, f64) {
    match face {
        0 => (r.y / r.x, r.z / r.x),
        1 => (-r.x / r.y, r.z / r.y),
        2 => (-r.x / r.z, -r.y / r.z),
        3 => (r.z / r.x, r.y / r.x),
        4 => (r.z / r.y, -r.x / r.y),
        _ => (-r.y / r.z, -r.x / r.z),
    }
}

/// Converts a direction vector (not necessarily unit length) to `(face, u, v)`.
pub(super) fn xyz_to_face_uv(r: R3Vector) -> (u8, f64, f64) {
    let f = face(r);
    let (u, v) = valid_face_xyz_to_uv(f, r);
    (f, u, v)
}

/// Converts `(face, u, v)` to an unnormalized direction vector.
/// `pub(crate)`: `shape::rect` inverts the projection too.
pub(crate) fn face_uv_to_xyz(face: u8, u: f64, v: f64) -> R3Vector {
    match face {
        0 => R3Vector::new(1.0, u, v),
        1 => R3Vector::new(-u, 1.0, v),
        2 => R3Vector::new(-u, -v, 1.0),
        3 => R3Vector::new(-1.0, -v, -u),
        4 => R3Vector::new(v, -1.0, -u),
        _ => R3Vector::new(v, u, -1.0),
    }
}
