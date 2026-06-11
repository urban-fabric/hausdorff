//! `s2` orientation predicates on the unit sphere, ported from golang/geo's
//! `predicates.go`. Computed from vector triple products — no projection.

use super::S2Point;

/// True iff `a, b, c` wind counter-clockwise, i.e. `c` lies to the left of the
/// directed great-circle edge `a → b`.
///
/// Computed as the signed volume of the parallelepiped, `(c × a) · b > 0` — the
/// order golang/geo's `Sign` uses so that `sign(a,b,c)` and `sign(c,b,a)` are
/// never both true.
#[must_use]
pub fn sign(a: S2Point, b: S2Point, c: S2Point) -> bool {
    c.vector().cross(a.vector()).dot(b.vector()) > 0.0
}

/// Orientation of `a, b, c`: `+1` counter-clockwise, `−1` clockwise, `0`
/// degenerate (two points equal, or exactly collinear in floating point).
///
/// This is the fast triage form; unlike golang/geo's full `RobustSign` it does
/// not fall back to exact arithmetic for points that are collinear to within
/// floating-point error (it returns `0` instead). Rotation-invariant and
/// sign-reversing under argument swaps for all non-degenerate inputs.
#[must_use]
pub fn robust_sign(a: S2Point, b: S2Point, c: S2Point) -> i32 {
    let det = c.vector().cross(a.vector()).dot(b.vector());
    if det > 0.0 {
        1
    } else if det < 0.0 {
        -1
    } else {
        0
    }
}
