//! `s2::Region` — the common interface over spherical regions ([`Cap`], [`Loop`],
//! [`Polygon`]), mirroring golang/geo's `Region`. A region tests point
//! containment, reports a bounding [`Cap`], and answers the cell-level
//! intersect/contains queries that [`RegionCoverer`] uses to approximate the
//! region by a set of S2 cells.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use hausdorff_base::space::r3::Vector;
use hausdorff_base::space::s1::Angle;

use super::{Cap, CellUnion, Crossing, Loop, Polygon, Rect, S2CellId, S2Point, crossing_sign};

/// The deepest S2 subdivision level (a leaf cell).
const MAX_LEVEL: u32 = 30;

/// A region on the unit sphere.
pub trait Region {
    /// True iff the region contains `p`.
    fn contains_point(&self, p: S2Point) -> bool;

    /// A spherical cap that bounds the region (not necessarily the tightest).
    fn bounding_cap(&self) -> Cap;

    /// A lat/lng rectangle that bounds the region (not necessarily the tightest).
    fn bounding_rect(&self) -> Rect {
        self.bounding_cap().rect_bound()
    }

    /// True iff the region may intersect `cell`. **Conservative:** an
    /// implementation must return `true` whenever the cell and region share any
    /// point, so the coverer never prunes a cell it should keep (over-covering at
    /// worst). The [`Loop`]/[`Polygon`] impls achieve this by counting any
    /// boundary contact — a crossing *or* a touch — plus either-way vertex
    /// nesting; the [`Cap`] impl is approximate (see its note).
    fn intersects_cell(&self, cell: S2CellId) -> bool;

    /// True iff the region fully contains `cell`. **Strict:** a `false` result
    /// only means "not known to be fully inside", so a coverer keeps subdividing
    /// rather than wrongly emitting a coarse cell that pokes outside the region.
    fn contains_cell(&self, cell: S2CellId) -> bool;
}

/// The cell as a counter-clockwise spherical [`Loop`] (interior to the left). The
/// cube→sphere projection can wind a cell either way depending on the face, so we
/// orient by the cell center: a CCW loop contains its own center.
pub(super) fn cell_loop(cell: S2CellId) -> Loop {
    let verts = cell.vertices();
    let l = Loop::new(verts.to_vec());
    if l.contains_point(cell.center_point()) {
        l
    } else {
        let mut rev = verts.to_vec();
        rev.reverse();
        Loop::new(rev)
    }
}

/// The directed edges of a loop `(vᵢ, vᵢ₊₁)`, wrapping at the end.
fn loop_edges(l: &Loop) -> impl Iterator<Item = (S2Point, S2Point)> + '_ {
    let v = l.vertices();
    let n = v.len();
    (0..n).map(move |i| (v[i], v[(i + 1) % n]))
}

/// True iff any edge of `a` meets any edge of `b` — a proper interior crossing or
/// a shared/touching vertex ([`Crossing::Maybe`]). Used to deny "fully contains"
/// the moment boundaries touch, so the coverer keeps subdividing (the safe way).
fn boundaries_meet(a: &Loop, b: &Loop) -> bool {
    if a.num_vertices() < 2 || b.num_vertices() < 2 {
        return false;
    }
    for (a0, a1) in loop_edges(a) {
        for (b0, b1) in loop_edges(b) {
            if crossing_sign(a0, a1, b0, b1) != Crossing::DoNotCross {
                return true;
            }
        }
    }
    false
}

/// True iff `inner` lies entirely within `outer`: every `inner` vertex is inside
/// `outer` and their boundaries do not meet.
fn loop_contains_loop(outer: &Loop, inner: &Loop) -> bool {
    inner.vertices().iter().all(|&v| outer.contains_point(v)) && !boundaries_meet(outer, inner)
}

/// True iff loops `a` and `b` share any point — a proper crossing, a boundary
/// *touch* ([`Crossing::Maybe`]), or one loop nested in the other. Conservative:
/// counting touches means a coverer never prunes a cell that merely grazes the
/// region.
fn loops_overlap(a: &Loop, b: &Loop) -> bool {
    if a.num_vertices() < 3 || b.num_vertices() < 3 {
        return false;
    }
    for (a0, a1) in loop_edges(a) {
        for (b0, b1) in loop_edges(b) {
            if crossing_sign(a0, a1, b0, b1) != Crossing::DoNotCross {
                return true;
            }
        }
    }
    b.vertices().iter().any(|&v| a.contains_point(v))
        || a.vertices().iter().any(|&v| b.contains_point(v))
}

/// The bounding cap of a set of points: centered on their normalized vector mean,
/// with radius equal to the farthest point's angular distance.
///
/// Proof-backed by `riemann` `Hausdorff/RegionBound`: the max-radius cap about a
/// centre **contains every point** (`boundingCap_contains`, soundness — so the
/// coverer never starts from a cap that omits part of the region) and is
/// **tightest at that centre** (`boundingCap_minimal`). (The normalized-mean
/// centre is the uncertified tightness heuristic; the radius is what makes the
/// bound correct.)
pub(super) fn bounding_cap_of(points: &[S2Point]) -> Cap {
    if points.is_empty() {
        return Cap::empty(S2Point::from_vector(Vector::new(1.0, 0.0, 0.0)));
    }
    let mut sum = Vector::new(0.0, 0.0, 0.0);
    for p in points {
        sum = sum.add(p.vector());
    }
    let mean = sum.normalize();
    let center = if mean.norm2() == 0.0 {
        points[0]
    } else {
        S2Point::from_vector(mean)
    };
    let radius = points
        .iter()
        .map(|p| center.vector().angle(p.vector()).radians())
        .fold(0.0_f64, f64::max);
    Cap::from_center_angle(center, Angle::from_radians(radius))
}

impl Region for Cap {
    fn contains_point(&self, p: S2Point) -> bool {
        Cap::contains_point(*self, p)
    }
    fn bounding_cap(&self) -> Cap {
        *self
    }
    fn bounding_rect(&self) -> Rect {
        self.rect_bound()
    }
    fn intersects_cell(&self, cell: S2CellId) -> bool {
        let verts = cell.vertices();
        verts.iter().any(|&v| self.contains_point(v)) || self.intersects(bounding_cap_of(&verts))
    }
    fn contains_cell(&self, cell: S2CellId) -> bool {
        cell.vertices().iter().all(|&v| self.contains_point(v))
    }
}

impl Region for Loop {
    fn contains_point(&self, p: S2Point) -> bool {
        Loop::contains_point(self, p)
    }
    fn bounding_cap(&self) -> Cap {
        bounding_cap_of(self.vertices())
    }
    fn bounding_rect(&self) -> Rect {
        self.rect_bound()
    }
    fn intersects_cell(&self, cell: S2CellId) -> bool {
        loops_overlap(self, &cell_loop(cell))
    }
    fn contains_cell(&self, cell: S2CellId) -> bool {
        loop_contains_loop(self, &cell_loop(cell))
    }
}

impl Region for Polygon {
    fn contains_point(&self, p: S2Point) -> bool {
        Polygon::contains_point(self, p)
    }
    fn bounding_cap(&self) -> Cap {
        bounding_cap_of(self.shell().vertices())
    }
    fn bounding_rect(&self) -> Rect {
        self.rect_bound()
    }
    fn intersects_cell(&self, cell: S2CellId) -> bool {
        let cl = cell_loop(cell);
        loops_overlap(self.shell(), &cl) && !self.holes().iter().any(|h| loop_contains_loop(h, &cl))
    }
    fn contains_cell(&self, cell: S2CellId) -> bool {
        let cl = cell_loop(cell);
        loop_contains_loop(self.shell(), &cl) && !self.holes().iter().any(|h| loops_overlap(h, &cl))
    }
}

mod coverer;
pub use coverer::RegionCoverer;

#[cfg(test)]
mod tests;
