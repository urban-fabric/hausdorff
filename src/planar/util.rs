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
pub(crate) fn signed_area2(ring: &[Coord]) -> f64 {
    let v = ring_vertices(ring);
    let n = v.len();
    if n < 3 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..n {
        let a = v[i];
        let b = v[(i + 1) % n];
        s += a.x * b.y - b.x * a.y;
    }
    s
}

/// Even-odd ray-cast test: is `p` inside `ring`? Reports interior membership;
/// results for points exactly on the boundary are not guaranteed either way.
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
