//! Turning the sweep's result rings into polygons: the contour-assembly phase
//! of [`bool_ops`](super). Operates purely on coordinate rings (no sweep
//! events) — it classifies rings by containment depth into shells and holes and
//! normalises their winding.

use super::super::util::point_in_ring;
use super::super::{Coord, LineString, Polygon};
use super::ring_vertices;

/// Turn raw closed rings into polygons. The Martinez output can wind either way,
/// so classify by *containment depth* rather than signed area: a ring contained
/// in an even number of others is a shell, odd is a hole of its immediate
/// container. Orientation is then normalised (shells CCW, holes CW).
pub(super) fn assemble(rings: Vec<Vec<Coord>>) -> Vec<Polygon> {
    let n = rings.len();
    if n == 0 {
        return Vec::new();
    }
    let reps: Vec<Coord> = rings.iter().map(|r| interior_point(r)).collect();
    let depth: Vec<usize> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i && point_in_ring(reps[i], &rings[j]))
                .count()
        })
        .collect();

    let mut polys: Vec<Polygon> = Vec::new();
    for i in 0..n {
        if !depth[i].is_multiple_of(2) {
            continue;
        }
        let holes: Vec<LineString> = (0..n)
            .filter(|&h| {
                depth[h] == depth[i] + 1 && immediate_parent(&rings, &reps, &depth, h) == Some(i)
            })
            .map(|h| LineString::new(oriented(rings[h].clone(), false)))
            .collect();
        polys.push(Polygon::new(
            LineString::new(oriented(rings[i].clone(), true)),
            holes,
        ));
    }
    polys
}

/// The innermost ring that contains hole `h` at one level shallower.
fn immediate_parent(
    rings: &[Vec<Coord>],
    reps: &[Coord],
    depth: &[usize],
    h: usize,
) -> Option<usize> {
    let mut best: Option<usize> = None;
    for j in 0..rings.len() {
        if j != h && depth[j] + 1 == depth[h] && point_in_ring(reps[h], &rings[j]) {
            best = match best {
                Some(b) if ring_area2(&rings[b]).abs() <= ring_area2(&rings[j]).abs() => Some(b),
                _ => Some(j),
            };
        }
    }
    best
}

/// Reverse the ring if its winding doesn't match the requested orientation
/// (`ccw = true` for shells, `false` for holes).
fn oriented(mut ring: Vec<Coord>, ccw: bool) -> Vec<Coord> {
    if (ring_area2(&ring) > 0.0) != ccw {
        ring.reverse();
    }
    ring
}

/// A point strictly inside the ring: step inward from the bottom-most (hence
/// convex) vertex along its edge bisector.
fn interior_point(ring: &[Coord]) -> Coord {
    let v = ring_vertices(ring);
    let n = v.len();
    if n < 3 {
        return v.first().copied().unwrap_or(Coord::new(0.0, 0.0));
    }
    let mut k = 0;
    for i in 1..n {
        if v[i].y < v[k].y || (v[i].y == v[k].y && v[i].x < v[k].x) {
            k = i;
        }
    }
    let cur = v[k];
    let e1 = v[(k + n - 1) % n] - cur;
    let e2 = v[(k + 1) % n] - cur;
    let (u1, l1) = unit(e1);
    let (u2, l2) = unit(e2);
    let mut bis = Coord::new(u1.x + u2.x, u1.y + u2.y);
    let blen = (bis.x * bis.x + bis.y * bis.y).sqrt();
    bis = if blen < 1e-12 {
        Coord::new(-u1.y, u1.x)
    } else {
        Coord::new(bis.x / blen, bis.y / blen)
    };
    let step = l1.min(l2) * 1e-3;
    let p = Coord::new(cur.x + bis.x * step, cur.y + bis.y * step);
    if point_in_ring(p, ring) {
        p
    } else {
        Coord::new(cur.x - bis.x * step, cur.y - bis.y * step)
    }
}

fn unit(c: Coord) -> (Coord, f64) {
    let l = (c.x * c.x + c.y * c.y).sqrt();
    if l < 1e-12 {
        (Coord::new(0.0, 0.0), 0.0)
    } else {
        (Coord::new(c.x / l, c.y / l), l)
    }
}

fn ring_area2(ring: &[Coord]) -> f64 {
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
