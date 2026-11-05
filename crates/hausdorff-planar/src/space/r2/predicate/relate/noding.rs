//! Segment noding and the areal spatial index: padded bounding boxes and the
//! segment STR-tree ([`padded_box`], [`segment_index`], [`node_segments`]),
//! pairwise segment intersection ([`collect_intersections`], [`proper_cross`],
//! [`rings_cross`]), coordinate dedup, and [`ArealIndex`] — the ring-edge index
//! behind `O(log n)` areal point location.

use super::decompose::Seg;
use super::predicates::point_on_segment;
use super::{Decomposed, Loc};
use crate::space::r2::util::{orient, ring_vertices};
use crate::space::r2::{Coord, Polygon, Rect, StrTree};

/// Do any ring edge of `ga` and ring edge of `gb` properly cross (strict interior
/// crossing, not a shared vertex or collinear touch)?
pub(super) fn rings_cross(ga: &Decomposed, gb: &Decomposed) -> bool {
    let b_edges: Vec<(Rect, (Coord, Coord))> = gb
        .segments
        .iter()
        .filter(|s| s.owner == Loc::Boundary)
        .map(|s| (padded_box(s.a, s.b), (s.a, s.b)))
        .collect();
    if b_edges.is_empty() {
        return false;
    }
    let tree = StrTree::new(b_edges);
    let mut crossed = false;
    for sa in ga.segments.iter().filter(|s| s.owner == Loc::Boundary) {
        if crossed {
            break;
        }
        tree.query_each(padded_box(sa.a, sa.b), |_, &(b1, b2)| {
            if proper_cross(sa.a, sa.b, b1, b2) {
                crossed = true;
            }
        });
    }
    crossed
}

fn proper_cross(a1: Coord, a2: Coord, b1: Coord, b2: Coord) -> bool {
    let (d1, d2, d3, d4) = (
        orient(b1, b2, a1),
        orient(b1, b2, a2),
        orient(a1, a2, b1),
        orient(a1, a2, b2),
    );
    d1 != 0.0
        && d2 != 0.0
        && d3 != 0.0
        && d4 != 0.0
        && (d1 > 0.0) != (d2 > 0.0)
        && (d3 > 0.0) != (d4 > 0.0)
}

/// A segment's bounding box, expanded by a scale-relative epsilon so that the
/// box test never excludes a pair that the exact (eps-tolerant) intersection
/// test would accept — the box filter stays a safe superset.
pub(super) fn padded_box(a: Coord, b: Coord) -> Rect {
    let e = 1e-7 * (1.0 + a.x.abs() + a.y.abs() + b.x.abs() + b.y.abs());
    Rect::new(
        Coord::new(a.x.min(b.x) - e, a.y.min(b.y) - e),
        Coord::new(a.x.max(b.x) + e, a.y.max(b.y) + e),
    )
}

/// Index a segment list by padded bounding box for spatial queries.
fn segment_index(segs: &[Seg]) -> StrTree<usize> {
    StrTree::new(
        segs.iter()
            .enumerate()
            .map(|(i, s)| (padded_box(s.a, s.b), i))
            .collect(),
    )
}

/// Node two segment sets: collect every pairwise intersection point, but only
/// test pairs whose bounding boxes overlap (others cannot intersect, so this is
/// identical to the `O(n·m)` double loop — just without the misses). Indexes the
/// larger set and queries the smaller, so the cost is `O((n+m)·log)`.
pub(super) fn node_segments(sa: &[Seg], sb: &[Seg], out: &mut Vec<Coord>) {
    if sa.is_empty() || sb.is_empty() {
        return;
    }
    if sa.len() <= sb.len() {
        let tree = segment_index(sb);
        for s in sa {
            tree.query_each(padded_box(s.a, s.b), |_, &j| {
                collect_intersections(s.a, s.b, sb[j].a, sb[j].b, out);
            });
        }
    } else {
        let tree = segment_index(sa);
        for t in sb {
            tree.query_each(padded_box(t.a, t.b), |_, &j| {
                collect_intersections(sa[j].a, sa[j].b, t.a, t.b, out);
            });
        }
    }
}

/// Spatial index of a geometry's polygon ring edges for `O(log n)` areal point
/// location. A point-box query answers "on a ring?" (boundary); a rightward-ray
/// query gives each polygon's crossing parity ("inside?"). Reproduces exactly
/// the result of scanning every ring with [`locate_in_polygon`](super::predicates::locate_in_polygon).
pub(super) struct ArealIndex {
    /// Ring edge `(a, b)` tagged with its owning polygon's index.
    tree: StrTree<(Coord, Coord, u32)>,
    n_polys: usize,
}

impl ArealIndex {
    pub(super) fn build(polygons: &[Polygon]) -> Option<ArealIndex> {
        if polygons.is_empty() {
            return None;
        }
        let mut items: Vec<(Rect, (Coord, Coord, u32))> = Vec::new();
        for (pi, poly) in polygons.iter().enumerate() {
            let rings =
                std::iter::once(&poly.exterior().0).chain(poly.interiors().iter().map(|r| &r.0));
            for ring in rings {
                let v = ring_vertices(ring);
                let n = v.len();
                for i in 0..n {
                    let (a, b) = (v[i], v[(i + 1) % n]);
                    if a != b {
                        items.push((padded_box(a, b), (a, b, pi as u32)));
                    }
                }
            }
        }
        (!items.is_empty()).then(|| ArealIndex {
            tree: StrTree::new(items),
            n_polys: polygons.len(),
        })
    }

    /// Boundary if `p` lies on any ring, else Interior if it is inside any
    /// polygon (per-polygon even-odd over that polygon's rings, so holes
    /// subtract), else Exterior. Interior of one polygon wins over Boundary of
    /// another, matching [`locate_areal`](super::predicates::locate_areal).
    pub(super) fn locate(&self, p: Coord) -> Loc {
        let bq = padded_box(p, p);
        let ray = Rect::new(p, Coord::new(f64::INFINITY, p.y));

        if self.n_polys == 1 {
            let mut on = false;
            self.tree
                .query_each(bq, |_, &(a, b, _)| on |= point_on_segment(p, a, b));
            if on {
                return Loc::Boundary;
            }
            let mut crossings = 0u32;
            self.tree.query_each(ray, |_, &(a, b, _)| {
                crossings += u32::from(ray_crosses(p, a, b))
            });
            return if crossings % 2 == 1 {
                Loc::Interior
            } else {
                Loc::Exterior
            };
        }

        let mut on_poly = vec![false; self.n_polys];
        self.tree.query_each(bq, |_, &(a, b, poly)| {
            on_poly[poly as usize] |= point_on_segment(p, a, b);
        });
        let mut cross = vec![0u32; self.n_polys];
        self.tree.query_each(ray, |_, &(a, b, poly)| {
            cross[poly as usize] += u32::from(ray_crosses(p, a, b));
        });
        let mut any_boundary = false;
        for poly in 0..self.n_polys {
            if on_poly[poly] {
                any_boundary = true;
            } else if cross[poly] % 2 == 1 {
                return Loc::Interior;
            }
        }
        if any_boundary {
            Loc::Boundary
        } else {
            Loc::Exterior
        }
    }
}

/// Does a rightward horizontal ray from `p` cross segment `a`–`b`? Half-open in
/// `y` (`[min_y, max_y)`), the crossing-number convention [`point_in_ring`](crate::space::r2::util::point_in_ring) uses
/// (a shared vertex counts once).
fn ray_crosses(p: Coord, a: Coord, b: Coord) -> bool {
    let (lo, hi) = if a.y <= b.y { (a, b) } else { (b, a) };
    if p.y < lo.y || p.y >= hi.y {
        return false;
    }
    let t = (p.y - lo.y) / (hi.y - lo.y);
    lo.x + t * (hi.x - lo.x) > p.x
}

/// Proper + collinear intersection points of two segments.
fn collect_intersections(a1: Coord, a2: Coord, b1: Coord, b2: Coord, out: &mut Vec<Coord>) {
    let d1 = orient(b1, b2, a1);
    let d2 = orient(b1, b2, a2);
    let d3 = orient(a1, a2, b1);
    let d4 = orient(a1, a2, b2);
    if ((d1 > 0.0) != (d2 > 0.0)) && ((d3 > 0.0) != (d4 > 0.0)) && d1 != 0.0 && d2 != 0.0 {
        let denom = (a1.x - a2.x) * (b1.y - b2.y) - (a1.y - a2.y) * (b1.x - b2.x);
        if denom != 0.0 {
            let t = ((a1.x - b1.x) * (b1.y - b2.y) - (a1.y - b1.y) * (b1.x - b2.x)) / denom;
            out.push(Coord::new(
                a1.x + t * (a2.x - a1.x),
                a1.y + t * (a2.y - a1.y),
            ));
        }
    }
    for (p, q1, q2) in [(a1, b1, b2), (a2, b1, b2), (b1, a1, a2), (b2, a1, a2)] {
        if point_on_segment(p, q1, q2) {
            out.push(p);
        }
    }
}

pub(super) fn dedup_coords(pts: &mut Vec<Coord>) {
    pts.retain(|c| c.is_finite());
    pts.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap().then(a.y.total_cmp(&b.y)));
    pts.dedup();
}
