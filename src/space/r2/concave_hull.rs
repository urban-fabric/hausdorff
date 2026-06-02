//! `r2::concave_hull` — a concave hull via the chi-shape "edge digging"
//! heuristic: start from the convex hull and repeatedly replace a boundary edge
//! with a detour through a nearby interior point, as long as that does not create
//! a self-intersection.
//!
//! `concavity` follows geo's convention: **lower is more concave** (more
//! detail); large values approach the convex hull. This is a heuristic hull (not
//! geo's `concaveman`/R-tree port), but it is dependency-free and produces a
//! valid simple polygon that covers all the input points.

use super::convex_hull::convex_hull_coords;
use super::util::orient;
use super::{Coord, LineString, MultiPoint, Polygon};

/// Build a concave hull at the given `concavity`.
pub trait ConcaveHull {
    fn concave_hull(&self, concavity: f64) -> Polygon;
}

impl ConcaveHull for MultiPoint {
    fn concave_hull(&self, concavity: f64) -> Polygon {
        let coords: Vec<Coord> = self.0.iter().map(|p| p.0).collect();
        Polygon::new(
            LineString::new(concave_hull_coords(&coords, concavity)),
            vec![],
        )
    }
}

/// Closed (first vertex repeated at the end) concave-hull ring of a point set.
pub(crate) fn concave_hull_coords(points: &[Coord], concavity: f64) -> Vec<Coord> {
    let mut pts: Vec<Coord> = points.iter().copied().filter(|c| c.is_finite()).collect();
    pts.sort_by(|a, b| {
        a.x.partial_cmp(&b.x)
            .unwrap()
            .then(a.y.partial_cmp(&b.y).unwrap())
    });
    pts.dedup();
    if pts.len() < 4 {
        return convex_hull_coords(&pts);
    }

    // open boundary from the convex hull
    let mut boundary = convex_hull_coords(&pts);
    if boundary.len() >= 2 && boundary[0] == *boundary.last().unwrap() {
        boundary.pop();
    }
    let on_hull: std::collections::HashSet<(u64, u64)> = boundary.iter().map(|c| key(*c)).collect();
    let mut interior: Vec<Coord> = pts
        .iter()
        .copied()
        .filter(|c| !on_hull.contains(&key(*c)))
        .collect();

    let concavity = concavity.max(1.0);
    let mut changed = true;
    while changed && !interior.is_empty() {
        changed = false;
        let mut i = 0;
        while i < boundary.len() {
            let n = boundary.len();
            let a = boundary[i];
            let b = boundary[(i + 1) % n];
            let edge_len = dist(a, b);

            if let Some((idx, d2)) = nearest_to_edge(&interior, a, b) {
                let d = d2.sqrt();
                // dig when the edge is long relative to the detour distance
                if d > 0.0 && edge_len / d >= concavity {
                    let pt = interior[idx];
                    if !insert_self_intersects(&boundary, i, pt) {
                        boundary.insert(i + 1, pt);
                        interior.remove(idx);
                        changed = true;
                        continue; // re-examine the new edge (a → pt)
                    }
                }
            }
            i += 1;
        }
    }

    boundary.push(boundary[0]); // close
    boundary
}

fn key(c: Coord) -> (u64, u64) {
    (c.x.to_bits(), c.y.to_bits())
}

fn dist(a: Coord, b: Coord) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt()
}

/// Squared distance from `p` to segment `a–b`.
fn point_seg_d2(p: Coord, a: Coord, b: Coord) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    let t = if len2 == 0.0 {
        0.0
    } else {
        (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0)
    };
    let (px, py) = (a.x + t * dx, a.y + t * dy);
    (p.x - px).powi(2) + (p.y - py).powi(2)
}

fn nearest_to_edge(interior: &[Coord], a: Coord, b: Coord) -> Option<(usize, f64)> {
    interior
        .iter()
        .enumerate()
        .map(|(i, p)| (i, point_seg_d2(*p, a, b)))
        .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
}

/// Proper crossing: segments share an interior point (no shared endpoints).
fn proper_cross(a: Coord, b: Coord, c: Coord, d: Coord) -> bool {
    let (d1, d2, d3, d4) = (
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    );
    d1 != 0.0
        && d2 != 0.0
        && d3 != 0.0
        && d4 != 0.0
        && (d1 > 0.0) != (d2 > 0.0)
        && (d3 > 0.0) != (d4 > 0.0)
}

/// Would inserting `pt` after boundary index `i` make the new edges cross any
/// other boundary edge?
fn insert_self_intersects(boundary: &[Coord], i: usize, pt: Coord) -> bool {
    let n = boundary.len();
    let a = boundary[i];
    let b = boundary[(i + 1) % n];
    for j in 0..n {
        if j == i {
            continue; // the edge being split
        }
        let c = boundary[j];
        let d = boundary[(j + 1) % n];
        if proper_cross(a, pt, c, d) || proper_cross(pt, b, c, d) {
            return true;
        }
    }
    false
}
