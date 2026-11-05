//! Shared low-level helpers for the planar algorithms: ring normalization, the
//! shoelace kernel, and a point-in-ring test.

use super::Coord;

/// A ring's distinct vertices: drops a trailing coordinate equal to the first
/// (the closing point of a closed ring) so shoelace/centroid sums don't
/// double-count the seam.
pub(crate) fn ring_vertices(ring: &[Coord]) -> &[Coord] {
    if ring.len() >= 2 && ring[0] == ring[ring.len() - 1] {
        &ring[..ring.len() - 1]
    } else {
        ring
    }
}

/// Twice the signed shoelace area of a ring. Sign is positive for a
/// counter-clockwise ring in x-right / y-up axes. Works on closed or open rings;
/// degenerate rings (< 3 distinct vertices) give `0`.
///
/// **Vertices are shifted to a local origin first, and that is not an
/// optimisation.** The shoelace is translation-invariant in exact arithmetic but
/// badly conditioned in `f64`: the raw form sums terms of size `|coord|²` to reach
/// an answer of size `area`, so a small shape far from the origin cancels almost
/// everything it computes. Once `|coord|² > 2⁵³` the products are no longer exact
/// integers and a unit-sized contribution falls below their ULP, which puts the
/// threshold at about `1e8` — well inside real projected coordinates, and
/// millimetre-based grids reach it easily.
///
/// Untranslated, a unit square reported area `0` at an offset of `1e8`, `16384` at
/// `1e10` and `134217728` at `1e12`, and [`Winding::is_ccw`](super::Winding) — which
/// reads this function's *sign* — answered `false` for a counter-clockwise ring at
/// every one of those. Subtracting `v[0]` makes each term the size of the shape
/// rather than of its position, which is exact for the differences involved.
pub(crate) fn signed_area2(ring: &[Coord]) -> f64 {
    let v = ring_vertices(ring);
    let n = v.len();
    if n < 3 {
        return 0.0;
    }
    let o = v[0];
    let mut s = 0.0;
    for i in 0..n {
        let a = v[i];
        let b = v[(i + 1) % n];
        s += (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
    }
    s
}

/// Twice the signed area of triangle `a,b,c` — the 2-D cross product of `b-a`
/// and `c-a`. `> 0` if `a→b→c` turns left (counter-clockwise), `< 0` right,
/// `0` collinear.
pub(crate) fn orient(a: Coord, b: Coord, c: Coord) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// Even-odd ray-cast test: is `p` inside `ring`? Reports interior membership;
/// results for points exactly on the boundary are not guaranteed either way.
///
/// Proof-backed by `riemann` `ExactArith/PointInPolygon`: each edge's crossing
/// test is exactly a rightward-ray segment crossing (`crossesRay_iff_segment`,
/// for `p` off the endpoint heights — hence the boundary caveat), the verdict is
/// independent of edge direction (`crossesRay_symm`) and of the start vertex /
/// ring rotation (`rayParity_perm`, since the `inside` toggle is a permutation-
/// invariant `xor` fold). The remaining parity-equals-interior step is the
/// Jordan-curve theorem (out of scope).
pub(crate) fn point_in_ring(p: Coord, ring: &[Coord]) -> bool {
    let v = ring_vertices(ring);
    let n = v.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (pi, pj) = (v[i], v[j]);
        let crosses = (pi.y > p.y) != (pj.y > p.y)
            && p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x;
        if crosses {
            inside = !inside;
        }
        j = i;
    }
    inside
}
