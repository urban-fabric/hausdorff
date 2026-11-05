//! Native DE-9IM topological relate — the in-house replacement for `geo`'s
//! `Relate`. Computes the 3×3 intersection matrix of two geometries'
//! Interior/Boundary/Exterior point sets, then answers the OGC predicates
//! (`within`/`contains`/`covers`/`coveredby`/`touches`/`crosses`/`equals`/
//! `overlaps`/`disjoint`) from it.
//!
//! The matrix is built by *noding* the two geometries' 1-D edges against each
//! other, then labelling every resulting sub-edge and node with its location in
//! each operand (Interior/Boundary/Exterior) and recording the dimension of each
//! Location×Location cell. Areal interiors are handled with point-in-polygon.

use crate::space::r2::{Coord, Geometry, Rect, StrTree};

mod decompose;
mod geom;
mod matrix;
mod noding;
mod predicate;
mod predicates;
#[cfg(test)]
mod tests;
pub use matrix::IntersectionMatrix;
use matrix::{B, E, I, Loc};
pub use predicate::SpatialPredicate;

use decompose::Decomposed;
use geom::{boxes_strictly_disjoint, dimension, disjoint_matrix, geom_bbox};
use noding::{dedup_coords, node_segments, padded_box, rings_cross};
use predicates::{point_on_segment_interior, sample_interior};

/// Public entry: the DE-9IM matrix of `a` relative to `b`.
#[must_use]
pub fn relate(a: &Geometry, b: &Geometry) -> IntersectionMatrix {
    let _span = crate::profile::span("relate");
    if let (Some(ra), Some(rb)) = (geom_bbox(a), geom_bbox(b))
        && boxes_strictly_disjoint(ra, rb)
        && let Some(im) = disjoint_matrix(a, b)
    {
        return im;
    }

    let (ga, gb) = {
        let _s = crate::profile::span("decompose");
        (Decomposed::new(a), Decomposed::new(b))
    };
    let mut im = IntersectionMatrix::empty();
    im.set_at_least(E, E, 2);

    let mut pts: Vec<Coord> = Vec::new();
    {
        let _s = crate::profile::span("node");
        pts.extend(ga.points.iter().copied());
        pts.extend(gb.points.iter().copied());
        pts.extend(ga.all_vertices());
        pts.extend(gb.all_vertices());
        node_segments(&ga.segments, &gb.segments, &mut pts);
        dedup_coords(&mut pts);
    }

    {
        let _s = crate::profile::span("locate_points");
        for &p in &pts {
            let la = ga.locate(p);
            let lb = gb.locate(p);
            im.set_at_least(la as usize, lb as usize, 0);
        }
    }

    {
        let _s = crate::profile::span("edges");
        let pts_index: StrTree<Coord> =
            StrTree::new(pts.iter().map(|&p| (Rect::new(p, p), p)).collect());
        contribute_edges(&ga, &gb, &pts_index, &mut im, false);
        contribute_edges(&gb, &ga, &pts_index, &mut im, true);
    }

    {
        let _s = crate::profile::span("areas");
        contribute_areas(&ga, &gb, &mut im);
    }

    im
}

fn contribute_edges(
    from: &Decomposed,
    other: &Decomposed,
    pts_index: &StrTree<Coord>,
    im: &mut IntersectionMatrix,
    swapped: bool,
) {
    for seg in &from.segments {
        let (a, b) = (seg.a, seg.b);
        let mut on: Vec<Coord> = vec![a, b];
        pts_index.query_each(padded_box(a, b), |_, &p| {
            if point_on_segment_interior(p, a, b) {
                on.push(p);
            }
        });
        on.sort_by(|p, q| {
            let dp = (p.x - a.x).powi(2) + (p.y - a.y).powi(2);
            let dq = (q.x - a.x).powi(2) + (q.y - a.y).powi(2);
            dp.total_cmp(&dq)
        });
        on.dedup();
        for w in on.windows(2) {
            let mid = Coord::new((w[0].x + w[1].x) / 2.0, (w[0].y + w[1].y) / 2.0);
            let lf = seg.owner;
            let lo = other.locate(mid);
            let (la, lb) = if swapped { (lo, lf) } else { (lf, lo) };
            im.set_at_least(la as usize, lb as usize, 1);
        }
    }
}

/// 2-D contributions: areal interior overlaps (via native overlay), and the
/// interior/boundary ∩ exterior cells.
fn contribute_areas(ga: &Decomposed, gb: &Decomposed, im: &mut IntersectionMatrix) {
    let a_areal = !ga.polygons.is_empty();
    let b_areal = !gb.polygons.is_empty();

    if a_areal && b_areal {
        let cross = rings_cross(ga, gb);
        let (mut a_in, mut a_out) = (false, false);
        for poly in &ga.polygons {
            sample_interior(poly, ga, |p| {
                match gb.locate_areal_only(p) {
                    Loc::Interior => a_in = true,
                    Loc::Exterior => a_out = true,
                    Loc::Boundary => {}
                }
                a_in && a_out
            });
            if a_in && a_out {
                break;
            }
        }
        let (mut b_in, mut b_out) = (false, false);
        for poly in &gb.polygons {
            sample_interior(poly, gb, |p| {
                match ga.locate_areal_only(p) {
                    Loc::Interior => b_in = true,
                    Loc::Exterior => b_out = true,
                    Loc::Boundary => {}
                }
                b_in && b_out
            });
            if b_in && b_out {
                break;
            }
        }
        if cross || a_in || b_in {
            im.set_at_least(I, I, 2);
        }
        if cross || a_out {
            im.set_at_least(I, E, 2);
        }
        if cross || b_out {
            im.set_at_least(E, I, 2);
        }
    } else {
        if a_areal {
            im.set_at_least(I, E, 2);
        }
        if b_areal {
            im.set_at_least(E, I, 2);
        }
    }

    add_exterior_cells(ga, gb, im, false);
    add_exterior_cells(gb, ga, im, true);
}

/// Record Boundary(from) ∩ Exterior(other) and low-D Interior(from) ∩
/// Exterior(other) cells (areal interiors are handled by [`contribute_areas`]).
fn add_exterior_cells(
    from: &Decomposed,
    other: &Decomposed,
    im: &mut IntersectionMatrix,
    swapped: bool,
) {
    let set = |im: &mut IntersectionMatrix, a: usize, dim: i8| {
        let (ra, rb) = if swapped { (E, a) } else { (a, E) };
        im.set_at_least(ra, rb, dim);
    };

    let bdim = match from.dim {
        2 => 1,
        1 => 0,
        _ => -1,
    };
    if bdim >= 0 {
        let bps = boundary_test_points(from);
        if bps.iter().any(|&p| other.locate(p) == Loc::Exterior) {
            set(im, B, bdim);
        }
    }

    if from.dim == 0
        && from
            .points
            .iter()
            .any(|&p| other.locate(p) == Loc::Exterior)
    {
        set(im, I, 0);
    }
    if from.dim == 1
        && from.segments.iter().any(|s| {
            let mid = Coord::new((s.a.x + s.b.x) / 2.0, (s.a.y + s.b.y) / 2.0);
            other.locate(mid) == Loc::Exterior
        })
    {
        set(im, I, 1);
    }
}

fn boundary_test_points(d: &Decomposed) -> Vec<Coord> {
    if d.dim == 2 {
        let mut v = Vec::new();
        for s in &d.segments {
            v.push(s.a);
            v.push(Coord::new((s.a.x + s.b.x) / 2.0, (s.a.y + s.b.y) / 2.0));
        }
        v
    } else {
        d.boundary_points()
    }
}

/// Public convenience mirroring how galois consumes the matrix.
#[must_use]
pub fn geometry_dimension(g: &Geometry) -> i8 {
    dimension(g)
}
