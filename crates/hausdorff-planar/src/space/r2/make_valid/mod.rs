//! `r2::make_valid` — repairing polygonal geometry under the odd-even fill rule.
//!
//! The rule, following Ledoux et al. (2014): a point is inside the repaired
//! geometry when a ray from it crosses the input's edges an odd number of
//! times, counting every edge of every ring separately. So a bowtie splits into
//! two components, a dangling spur cancels against itself, a hole that shares an
//! edge with its shell cancels there, and a ring nested two deep is an island
//! rather than a hole.
//!
//! The route to that region avoids the boolean sweep entirely, which matters:
//! repair is what a caller reaches for when geometry is *already* untrustworthy,
//! and it should not depend on the operation that most wants its output.
//!
//! Instead:
//!
//! 1. cut every edge against every other, on one shared vertex pool ([`node`]);
//! 2. keep the arrangement edges covered an **odd** number of times.
//!
//!    This is the whole of the fill rule. Crossing an edge flips parity once per
//!    covering input edge, so a piece covered an even number of times has the
//!    same parity on both sides and is interior to the region, not boundary. No
//!    ray casting, and no stepping to one side of an edge to ask which side is
//!    in — the two questions that make repair fragile on slivers;
//! 3. walk the survivors into rings ([`trace`]). Those edges bound a region, so
//!    every vertex has even degree and every walk closes;
//! 4. classify rings into shells and holes by containment depth.
//!
//! The result is checked for validity before it is returned. A repair that
//! cannot produce valid output reports failure rather than handing back
//! something that merely looks repaired.

mod node;
mod trace;

use super::util::point_in_ring;
use super::validation::Validation;
use super::{Coord, LineString, MultiPolygon, Polygon};

/// Repair polygonal geometry into a valid [`MultiPolygon`].
pub trait MakeValid {
    /// The repaired geometry under the odd-even rule.
    ///
    /// # Errors
    /// Returns [`MakeValidError`] when the input carries non-finite
    /// coordinates, or when the repair cannot be assembled into valid output.
    fn make_valid(&self) -> Result<MultiPolygon, MakeValidError>;
}

/// Why a repair did not produce valid geometry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MakeValidError {
    /// A coordinate was NaN or infinite; nothing meaningful can be built.
    NonFiniteCoordinate,
    /// The arrangement could not be walked into closed rings.
    Unassembled,
    /// Rings were produced but do not form valid geometry.
    StillInvalid,
}

impl std::fmt::Display for MakeValidError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let message = match self {
            MakeValidError::NonFiniteCoordinate => "geometry has a non-finite coordinate",
            MakeValidError::Unassembled => "repaired edges did not close into rings",
            MakeValidError::StillInvalid => "repaired rings are not valid geometry",
        };
        f.write_str(message)
    }
}

impl std::error::Error for MakeValidError {}

impl MakeValid for Polygon {
    fn make_valid(&self) -> Result<MultiPolygon, MakeValidError> {
        repair(std::slice::from_ref(self))
    }
}

impl MakeValid for MultiPolygon {
    fn make_valid(&self) -> Result<MultiPolygon, MakeValidError> {
        repair(&self.0)
    }
}

fn repair(polygons: &[Polygon]) -> Result<MultiPolygon, MakeValidError> {
    let mut segments: Vec<(Coord, Coord)> = Vec::new();
    let mut extent: f64 = 1.0;
    for polygon in polygons {
        for ring in std::iter::once(polygon.exterior()).chain(polygon.interiors()) {
            for pair in ring.0.windows(2) {
                let (a, b) = (pair[0], pair[1]);
                for c in [a, b] {
                    if !c.x.is_finite() || !c.y.is_finite() {
                        return Err(MakeValidError::NonFiniteCoordinate);
                    }
                    extent = extent.max(c.x.abs()).max(c.y.abs());
                }
                if a != b {
                    segments.push((a, b));
                }
            }
        }
    }
    if segments.is_empty() {
        return Ok(MultiPolygon::new(Vec::new()));
    }

    let mut nodes = node::Nodes::new(extent);
    let pieces = node::arrange(&segments, &mut nodes);

    // Odd cover is the fill rule: a piece under an even number of input edges
    // has matching parity either side of it, so it is interior, not boundary.
    let mut cover: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    for &(a, b) in &pieces {
        *cover.entry((a.min(b), a.max(b))).or_insert(0) += 1;
    }
    let mut boundary: Vec<(usize, usize)> = cover
        .into_iter()
        .filter(|&(_, count)| !count.is_multiple_of(2))
        .map(|(edge, _)| edge)
        .collect();
    // A hash map hands these back in a different order every process, and the
    // trace must not depend on that.
    boundary.sort_unstable();
    if boundary.is_empty() {
        return Ok(MultiPolygon::new(Vec::new()));
    }

    let rings = trace::rings(&boundary, &nodes).ok_or(MakeValidError::Unassembled)?;
    let repaired = assemble(rings);
    if repaired.is_valid() {
        Ok(repaired)
    } else {
        Err(MakeValidError::StillInvalid)
    }
}

/// Classify rings into shells and holes by containment depth — even depth is a
/// shell, odd is a hole of its immediate container — and normalise winding.
///
/// Depth rather than signed area because the walk's winding records the order
/// edges happened to be visited, not whether a ring encloses material.
fn assemble(rings: Vec<Vec<Coord>>) -> MultiPolygon {
    let n = rings.len();
    let closed: Vec<Vec<Coord>> = rings
        .iter()
        .map(|r| {
            let mut c = r.clone();
            c.push(c[0]);
            c
        })
        .collect();
    let reps: Vec<Coord> = rings.iter().map(|r| representative(r)).collect();
    let depth: Vec<usize> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i && point_in_ring(reps[i], &closed[j]))
                .count()
        })
        .collect();

    let mut polygons = Vec::new();
    for i in 0..n {
        if !depth[i].is_multiple_of(2) {
            continue;
        }
        let holes: Vec<LineString> = (0..n)
            .filter(|&h| depth[h] == depth[i] + 1 && parent(&closed, &reps, &depth, h) == Some(i))
            .map(|h| LineString::new(oriented(closed[h].clone(), false)))
            .collect();
        polygons.push(Polygon::new(
            LineString::new(oriented(closed[i].clone(), true)),
            holes,
        ));
    }
    MultiPolygon::new(polygons)
}

/// The smallest ring one level shallower that contains `h`.
fn parent(rings: &[Vec<Coord>], reps: &[Coord], depth: &[usize], h: usize) -> Option<usize> {
    let mut best: Option<usize> = None;
    for j in 0..rings.len() {
        if j != h && depth[j] + 1 == depth[h] && point_in_ring(reps[h], &rings[j]) {
            best = match best {
                Some(b) if area2(&rings[b]).abs() <= area2(&rings[j]).abs() => Some(b),
                _ => Some(j),
            };
        }
    }
    best
}

fn oriented(mut ring: Vec<Coord>, ccw: bool) -> Vec<Coord> {
    if (area2(&ring) > 0.0) != ccw {
        ring.reverse();
    }
    ring
}

fn area2(ring: &[Coord]) -> f64 {
    let n = ring.len().saturating_sub(1);
    let mut s = 0.0;
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        s += a.x * b.y - b.x * a.y;
    }
    s
}

/// A point strictly inside `ring`: step in from the lowest vertex, which is
/// necessarily convex, along the bisector of its two edges.
fn representative(ring: &[Coord]) -> Coord {
    let n = ring.len();
    let mut k = 0;
    for i in 1..n {
        if ring[i].y < ring[k].y || (ring[i].y == ring[k].y && ring[i].x < ring[k].x) {
            k = i;
        }
    }
    let cur = ring[k];
    let unit = |c: Coord| {
        let l = (c.x * c.x + c.y * c.y).sqrt();
        if l < f64::MIN_POSITIVE {
            (Coord::new(0.0, 0.0), 0.0)
        } else {
            (Coord::new(c.x / l, c.y / l), l)
        }
    };
    let (u1, l1) = unit(ring[(k + n - 1) % n] - cur);
    let (u2, l2) = unit(ring[(k + 1) % n] - cur);
    let mut bis = Coord::new(u1.x + u2.x, u1.y + u2.y);
    let blen = (bis.x * bis.x + bis.y * bis.y).sqrt();
    bis = if blen < 1e-12 {
        Coord::new(-u1.y, u1.x)
    } else {
        Coord::new(bis.x / blen, bis.y / blen)
    };
    let step = l1.min(l2) * 1e-3;
    let mut closed = ring.to_vec();
    closed.push(closed[0]);
    let p = Coord::new(cur.x + bis.x * step, cur.y + bis.y * step);
    if point_in_ring(p, &closed) {
        p
    } else {
        Coord::new(cur.x - bis.x * step, cur.y - bis.y * step)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests;
