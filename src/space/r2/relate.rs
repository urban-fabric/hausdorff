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

use super::util::ring_vertices;
use super::{Coord, Geometry, Polygon, Rect, StrTree};

/// Location of a point relative to a geometry.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Loc {
    Interior = 0,
    Boundary = 1,
    Exterior = 2,
}

/// A DE-9IM intersection matrix: `m[a][b]` is the dimension (`-1` = none, else
/// `0`/`1`/`2`) of `a`-part(A) ∩ `b`-part(B), with rows/cols ordered
/// Interior, Boundary, Exterior.
#[derive(Clone, Copy, Debug)]
pub struct IntersectionMatrix {
    m: [[i8; 3]; 3],
}

const I: usize = 0;
const B: usize = 1;
const E: usize = 2;

impl IntersectionMatrix {
    fn empty() -> Self {
        Self { m: [[-1; 3]; 3] }
    }

    #[inline]
    fn set_at_least(&mut self, a: usize, b: usize, dim: i8) {
        if dim > self.m[a][b] {
            self.m[a][b] = dim;
        }
    }

    #[inline]
    fn get(&self, a: usize, b: usize) -> i8 {
        self.m[a][b]
    }

    #[inline]
    fn is_t(&self, a: usize, b: usize) -> bool {
        self.m[a][b] >= 0
    }
    #[inline]
    fn is_f(&self, a: usize, b: usize) -> bool {
        self.m[a][b] < 0
    }

    // --- OGC predicates, matching geo's RelateMatrix method names ---

    #[must_use]
    pub fn is_disjoint(&self) -> bool {
        self.is_f(I, I) && self.is_f(I, B) && self.is_f(B, I) && self.is_f(B, B)
    }

    #[must_use]
    pub fn is_intersects(&self) -> bool {
        !self.is_disjoint()
    }

    /// `A` within `B`: `T*F**F***`.
    #[must_use]
    pub fn is_within(&self) -> bool {
        self.is_t(I, I) && self.is_f(I, E) && self.is_f(B, E)
    }

    /// `A` contains `B` — within with the operands swapped.
    #[must_use]
    pub fn is_contains(&self) -> bool {
        self.is_t(I, I) && self.is_f(E, I) && self.is_f(E, B)
    }

    /// `A` covers `B`: `T*****FF*` | `*T****FF*` | `***T**FF*` | `****T*FF*`.
    #[must_use]
    pub fn is_covers(&self) -> bool {
        (self.is_t(I, I) || self.is_t(I, B) || self.is_t(B, I) || self.is_t(B, B))
            && self.is_f(E, I)
            && self.is_f(E, B)
    }

    /// `A` covered by `B` — covers with the operands swapped.
    #[must_use]
    pub fn is_coveredby(&self) -> bool {
        (self.is_t(I, I) || self.is_t(I, B) || self.is_t(B, I) || self.is_t(B, B))
            && self.is_f(I, E)
            && self.is_f(B, E)
    }

    /// `A` equals `B` (same point set): `T*F**FFF*`.
    #[must_use]
    pub fn is_equal_topo(&self) -> bool {
        self.is_t(I, I)
            && self.is_f(I, E)
            && self.is_f(E, I)
            && self.is_f(B, E)
            && self.is_f(E, B)
    }

    /// `A` touches `B`: interiors disjoint, some boundary contact.
    /// `FT*******` | `F**T*****` | `F***T****`.
    #[must_use]
    pub fn is_touches(&self) -> bool {
        self.is_f(I, I) && (self.is_t(B, B) || self.is_t(I, B) || self.is_t(B, I))
    }

    /// `A` crosses `B` — dimension-dependent (OGC).
    // The dimension-pair clauses are enumerated symmetrically to mirror the OGC
    // table; factoring them (as clippy suggests) obscures that structure.
    #[allow(clippy::nonminimal_bool)]
    #[must_use]
    pub fn is_crosses(&self, dim_a: i8, dim_b: i8) -> bool {
        if (dim_a == 0 && dim_b == 1)
            || (dim_a == 0 && dim_b == 2)
            || (dim_a == 1 && dim_b == 2)
        {
            self.is_t(I, I) && self.is_t(I, E)
        } else if (dim_a == 1 && dim_b == 0)
            || (dim_a == 2 && dim_b == 0)
            || (dim_a == 2 && dim_b == 1)
        {
            self.is_t(I, I) && self.is_t(E, I)
        } else if dim_a == 1 && dim_b == 1 {
            self.get(I, I) == 0
        } else {
            false
        }
    }

    /// `A` overlaps `B` — dimension-dependent (OGC). Only defined for equal dims.
    #[must_use]
    pub fn is_overlaps(&self, dim_a: i8, dim_b: i8) -> bool {
        if dim_a != dim_b {
            return false;
        }
        if dim_a == 0 || dim_a == 2 {
            self.is_t(I, I) && self.is_t(I, E) && self.is_t(E, I)
        } else if dim_a == 1 {
            self.get(I, I) == 1 && self.is_t(I, E) && self.is_t(E, I)
        } else {
            false
        }
    }
}

/// Topological dimension of a geometry (max over components): 0/1/2, or -1 empty.
fn dimension(g: &Geometry) -> i8 {
    match g {
        Geometry::Point(_) | Geometry::MultiPoint(_) => 0,
        Geometry::Line(_) | Geometry::LineString(_) | Geometry::MultiLineString(_) => 1,
        Geometry::Polygon(_) | Geometry::MultiPolygon(_) | Geometry::Rect(_) | Geometry::Triangle(_) => 2,
        Geometry::GeometryCollection(gs) => gs.iter().map(dimension).max().unwrap_or(-1),
    }
}

/// Public entry: the DE-9IM matrix of `a` relative to `b`.
#[must_use]
pub fn relate(a: &Geometry, b: &Geometry) -> IntersectionMatrix {
    let _span = crate::profile::span("relate");
    // Fast path: strictly disjoint bounding boxes ⇒ disjoint geometries. The boxes
    // come from raw coordinates (no allocation), and when each operand's boundary
    // dimension is unambiguous we return the matrix *without* decomposing — the
    // common non-overlapping case, the way `geo`'s relate short-circuits.
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
    // Exterior ∩ Exterior is always 2-D in the plane (for non-empty geometries).
    im.set_at_least(E, E, 2);

    // Build the set of test points: all vertices of both, plus all pairwise
    // edge intersection points. Each test point is classified in A and B, and
    // contributes a 0-D cell. Edge interiors contribute 1-D cells; areas 2-D.
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

    // 0-D contributions from test points.
    {
        let _s = crate::profile::span("locate_points");
        for &p in &pts {
            let la = ga.locate(p);
            let lb = gb.locate(p);
            im.set_at_least(la as usize, lb as usize, 0);
        }
    }

    // 1-D contributions: sample the interior of each edge sub-segment between
    // consecutive split points, classifying its midpoint.
    {
        let _s = crate::profile::span("edges");
        contribute_edges(&ga, &gb, &pts, &mut im, false);
        contribute_edges(&gb, &ga, &pts, &mut im, true);
    }

    // 2-D contributions and exterior cells from interiors.
    {
        let _s = crate::profile::span("areas");
        contribute_areas(&ga, &gb, &mut im);
    }

    im
}

/// A 1-D edge tagged with its location in its *own* geometry: a line/curve edge
/// is part of the Interior; a polygon ring edge is part of the Boundary. Tagging
/// avoids re-deriving that location with floating-point geometry.
#[derive(Clone, Copy)]
struct Seg {
    a: Coord,
    b: Coord,
    owner: Loc,
}

/// Bounding box `(min, max)` of a geometry from its raw coordinates (no
/// decomposition), or `None` if it has none. Folds the raw coordinate slices
/// directly — unlike `coords_iter()`, which materializes a `Vec` — so the
/// disjoint-box fast path in [`relate`] allocates nothing.
fn geom_bbox(g: &Geometry) -> Option<(Coord, Coord)> {
    // Accumulate into `Coord`s rather than four scalars: a `Coord` is two packed
    // f64s, so `Coord::new(lo.x.min(c.x), lo.y.min(c.y))` lowers to a single
    // packed-min instruction (both lanes at once). Emptiness falls out of
    // `lo.x > hi.x` (the untouched ±∞ initial values), avoiding an `any` flag
    // whose per-step store would serialize the reduction.
    let mut lo = Coord::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Coord::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    fold_coords(g, &mut |c| {
        lo = Coord::new(lo.x.min(c.x), lo.y.min(c.y));
        hi = Coord::new(hi.x.max(c.x), hi.y.max(c.y));
    });
    (lo.x <= hi.x).then_some((lo, hi))
}

/// Visit every coordinate of `g` by value, without allocating an intermediate
/// `Vec` (the non-allocating counterpart of `coords_iter`). Generic over the
/// visitor so it inlines — a trait object here would make every coordinate a
/// virtual call and erase the win.
fn fold_coords<F: FnMut(Coord)>(g: &Geometry, f: &mut F) {
    match g {
        Geometry::Point(p) => f(p.0),
        Geometry::Line(l) => {
            f(l.start);
            f(l.end);
        }
        Geometry::LineString(ls) => {
            for &c in &ls.0 {
                f(c);
            }
        }
        Geometry::Polygon(p) => fold_polygon(p, f),
        Geometry::MultiPoint(mp) => {
            for p in mp.iter() {
                f(p.0);
            }
        }
        Geometry::MultiLineString(mls) => {
            for ls in mls.iter() {
                for &c in &ls.0 {
                    f(c);
                }
            }
        }
        Geometry::MultiPolygon(mp) => {
            for p in mp.iter() {
                fold_polygon(p, f);
            }
        }
        Geometry::Rect(r) => {
            let (lo, hi) = (r.min(), r.max());
            f(lo);
            f(Coord::new(hi.x, lo.y));
            f(hi);
            f(Coord::new(lo.x, hi.y));
        }
        Geometry::Triangle(t) => {
            f(t.0);
            f(t.1);
            f(t.2);
        }
        Geometry::GeometryCollection(gs) => {
            for sub in gs.iter() {
                fold_coords(sub, f);
            }
        }
    }
}

fn fold_polygon<F: FnMut(Coord)>(p: &Polygon, f: &mut F) {
    for &c in &p.exterior().0 {
        f(c);
    }
    for h in p.interiors() {
        for &c in &h.0 {
            f(c);
        }
    }
}

/// Do two boxes have a strict gap (so the geometries cannot even touch)?
fn boxes_strictly_disjoint(a: (Coord, Coord), b: (Coord, Coord)) -> bool {
    a.1.x < b.0.x || b.1.x < a.0.x || a.1.y < b.0.y || b.1.y < a.0.y
}

/// The DE-9IM matrix of two disjoint geometries — each operand's interior and
/// boundary lie wholly in the other's exterior. `None` when an operand's boundary
/// dimension is ambiguous (mod-2 line collections), signalling that the exact
/// decomposition path should run instead.
fn disjoint_matrix(a: &Geometry, b: &Geometry) -> Option<IntersectionMatrix> {
    let (bda, bdb) = (boundary_dim(a)?, boundary_dim(b)?);
    let mut im = IntersectionMatrix::empty();
    im.set_at_least(E, E, 2);
    im.set_at_least(I, E, dimension(a));
    if bda >= 0 {
        im.set_at_least(B, E, bda);
    }
    im.set_at_least(E, I, dimension(b));
    if bdb >= 0 {
        im.set_at_least(E, B, bdb);
    }
    Some(im)
}

/// Dimension of `g`'s boundary: `1` areal rings, `0` line endpoints, `-1` none
/// (point or closed line). `None` for mod-2 line collections, where it can't be
/// known without decomposing.
fn boundary_dim(g: &Geometry) -> Option<i8> {
    match g {
        Geometry::Polygon(_)
        | Geometry::MultiPolygon(_)
        | Geometry::Rect(_)
        | Geometry::Triangle(_) => Some(1),
        Geometry::Point(_) | Geometry::MultiPoint(_) => Some(-1),
        Geometry::Line(l) => Some(if l.start != l.end { 0 } else { -1 }),
        Geometry::LineString(ls) => Some(line_string_boundary_dim(&ls.0)),
        Geometry::MultiLineString(_) | Geometry::GeometryCollection(_) => None,
    }
}

/// An open linestring (distinct endpoints) has a 0-D boundary; a closed one none.
fn line_string_boundary_dim(coords: &[Coord]) -> i8 {
    if coords.len() >= 2 && coords[0] != coords[coords.len() - 1] {
        0
    } else {
        -1
    }
}

/// Decomposed view of a geometry: its 0-D points, boundary points (mod-2 for
/// lines), 1-D segments, and 2-D polygons.
struct Decomposed {
    points: Vec<Coord>,
    boundary: Vec<Coord>,
    segments: Vec<Seg>,
    polygons: Vec<Polygon>,
    dim: i8,
}

impl Decomposed {
    fn new(g: &Geometry) -> Self {
        let mut d = Decomposed {
            points: Vec::new(),
            boundary: Vec::new(),
            segments: Vec::new(),
            polygons: Vec::new(),
            dim: dimension(g),
        };
        d.add(g);
        d
    }

    fn add(&mut self, g: &Geometry) {
        match g {
            Geometry::Point(p) => self.points.push(p.0),
            Geometry::MultiPoint(mp) => self.points.extend(mp.0.iter().map(|p| p.0)),
            Geometry::Line(l) => {
                if l.start != l.end {
                    self.segments.push(Seg { a: l.start, b: l.end, owner: Loc::Interior });
                    self.add_line_boundary(&[l.start, l.end]);
                } else {
                    // degenerate zero-length line behaves as a point
                    self.points.push(l.start);
                }
            }
            Geometry::LineString(ls) => self.add_linestring(&ls.0),
            Geometry::MultiLineString(mls) => {
                for ls in mls.iter() {
                    self.add_linestring(&ls.0);
                }
            }
            Geometry::Polygon(p) => self.add_polygon(p),
            Geometry::MultiPolygon(mp) => {
                for p in mp.iter() {
                    self.add_polygon(p);
                }
            }
            Geometry::Rect(r) => self.add_polygon(&r.to_polygon()),
            Geometry::Triangle(t) => self.add_polygon(&t.to_polygon()),
            Geometry::GeometryCollection(gs) => {
                for sub in gs {
                    self.add(sub);
                }
            }
        }
    }

    fn add_linestring(&mut self, coords: &[Coord]) {
        for w in coords.windows(2) {
            if w[0] != w[1] {
                self.segments.push(Seg { a: w[0], b: w[1], owner: Loc::Interior });
            }
        }
        if coords.len() >= 2 {
            self.add_line_boundary(&[coords[0], coords[coords.len() - 1]]);
        }
    }

    /// Boundary endpoints follow the mod-2 rule: a point that is an endpoint of
    /// an odd number of curve components is on the boundary.
    fn add_line_boundary(&mut self, endpoints: &[Coord]) {
        for &e in endpoints {
            if e.is_finite() {
                self.boundary.push(e);
            }
        }
    }

    fn add_polygon(&mut self, p: &Polygon) {
        ring_segments(&p.exterior().0, &mut self.segments);
        for h in p.interiors() {
            ring_segments(&h.0, &mut self.segments);
        }
        self.polygons.push(p.clone());
    }

    fn all_vertices(&self) -> Vec<Coord> {
        let mut v: Vec<Coord> = self.points.clone();
        for s in &self.segments {
            v.push(s.a);
            v.push(s.b);
        }
        v
    }

    /// Effective boundary points after the mod-2 rule.
    fn boundary_points(&self) -> Vec<Coord> {
        let mut counts: Vec<(Coord, usize)> = Vec::new();
        for &b in &self.boundary {
            if let Some(c) = counts.iter_mut().find(|(p, _)| *p == b) {
                c.1 += 1;
            } else {
                counts.push((b, 1));
            }
        }
        counts
            .into_iter()
            .filter(|(_, n)| n % 2 == 1)
            .map(|(p, _)| p)
            .collect()
    }

    /// Location of `p` relative to this geometry.
    fn locate(&self, p: Coord) -> Loc {
        // areal: inside / on-ring / outside dominates
        if !self.polygons.is_empty() {
            let mut on_boundary = false;
            let mut inside = false;
            for poly in &self.polygons {
                match locate_in_polygon(p, poly) {
                    Loc::Boundary => on_boundary = true,
                    Loc::Interior => inside = true,
                    Loc::Exterior => {}
                }
            }
            if inside {
                return Loc::Interior;
            }
            if on_boundary {
                return Loc::Boundary;
            }
        }
        // linear: boundary endpoints, then segment interiors / vertices
        if !self.segments.is_empty() {
            let bpts = self.boundary_points();
            if bpts.contains(&p) {
                return Loc::Boundary;
            }
            for s in &self.segments {
                if p == s.a || p == s.b || point_on_segment_interior(p, s.a, s.b) {
                    return s.owner;
                }
            }
        }
        // puntal
        if self.points.contains(&p) {
            return Loc::Interior;
        }
        Loc::Exterior
    }
}

fn ring_segments(coords: &[Coord], out: &mut Vec<Seg>) {
    let v = ring_vertices(coords);
    let n = v.len();
    if n < 2 {
        return;
    }
    for i in 0..n {
        let (a, b) = (v[i], v[(i + 1) % n]);
        if a != b {
            out.push(Seg { a, b, owner: Loc::Boundary });
        }
    }
}

/// Classify a point against a set of polygons (areal location only).
fn locate_areal(p: Coord, polys: &[Polygon]) -> Loc {
    let mut inside = false;
    let mut on = false;
    for poly in polys {
        match locate_in_polygon(p, poly) {
            Loc::Interior => inside = true,
            Loc::Boundary => on = true,
            Loc::Exterior => {}
        }
    }
    if inside {
        Loc::Interior
    } else if on {
        Loc::Boundary
    } else {
        Loc::Exterior
    }
}

/// Do any ring edge of `ga` and ring edge of `gb` properly cross (strict interior
/// crossing, not a shared vertex or collinear touch)?
fn rings_cross(ga: &Decomposed, gb: &Decomposed) -> bool {
    // A proper crossing needs overlapping bounding boxes, so index gb's ring
    // edges and probe with ga's — only box-overlapping pairs are tested.
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

/// Several points strictly inside the polygon, spread around it (bottom-vertex
/// bisector, centroid, and each edge midpoint nudged inward).
fn interior_samples(poly: &Polygon) -> Vec<Coord> {
    let v = ring_vertices(&poly.exterior().0);
    let n = v.len();
    if n < 3 {
        return Vec::new();
    }
    let mut out = Vec::new();
    if let Some(p) = representative_interior_point(poly) {
        out.push(p);
    }
    let (cx, cy) = (
        v.iter().map(|c| c.x).sum::<f64>() / n as f64,
        v.iter().map(|c| c.y).sum::<f64>() / n as f64,
    );
    let centroid = Coord::new(cx, cy);
    if locate_in_polygon(centroid, poly) == Loc::Interior {
        out.push(centroid);
    }
    for i in 0..n {
        let (a, b) = (v[i], v[(i + 1) % n]);
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let len = dx.hypot(dy);
        if len < 1e-12 {
            continue;
        }
        let mid = Coord::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);
        let (nx, ny) = (-dy / len, dx / len);
        let step = len * 1e-3;
        for s in [step, -step] {
            let p = Coord::new(mid.x + nx * s, mid.y + ny * s);
            if locate_in_polygon(p, poly) == Loc::Interior {
                out.push(p);
                break;
            }
        }
    }
    out
}

/// A point guaranteed strictly inside the polygon (bottom-vertex bisector).
fn representative_interior_point(poly: &Polygon) -> Option<Coord> {
    let v = ring_vertices(&poly.exterior().0);
    let n = v.len();
    if n < 3 {
        return None;
    }
    let mut k = 0;
    for i in 1..n {
        if v[i].y < v[k].y || (v[i].y == v[k].y && v[i].x < v[k].x) {
            k = i;
        }
    }
    let cur = v[k];
    let prev = v[(k + n - 1) % n];
    let next = v[(k + 1) % n];
    let u = |a: Coord| {
        let l = (a.x * a.x + a.y * a.y).sqrt();
        if l < 1e-12 { (0.0, 0.0) } else { (a.x / l, a.y / l) }
    };
    let (u1x, u1y) = u(Coord::new(prev.x - cur.x, prev.y - cur.y));
    let (u2x, u2y) = u(Coord::new(next.x - cur.x, next.y - cur.y));
    let (mut bx, mut by) = (u1x + u2x, u1y + u2y);
    let bl = (bx * bx + by * by).sqrt();
    if bl < 1e-12 {
        bx = -u1y;
        by = u1x;
    } else {
        bx /= bl;
        by /= bl;
    }
    let elen = (prev.x - cur.x).hypot(prev.y - cur.y).min((next.x - cur.x).hypot(next.y - cur.y));
    let step = elen * 1e-3;
    let cand = Coord::new(cur.x + bx * step, cur.y + by * step);
    if locate_in_polygon(cand, poly) == Loc::Interior {
        Some(cand)
    } else {
        let alt = Coord::new(cur.x - bx * step, cur.y - by * step);
        (locate_in_polygon(alt, poly) == Loc::Interior).then_some(alt)
    }
}

/// Classify a point against a single polygon (ring boundary vs interior).
fn locate_in_polygon(p: Coord, poly: &Polygon) -> Loc {
    if point_on_ring(p, &poly.exterior().0) {
        return Loc::Boundary;
    }
    for h in poly.interiors() {
        if point_on_ring(p, &h.0) {
            return Loc::Boundary;
        }
    }
    if !point_in_ring(p, &poly.exterior().0) {
        return Loc::Exterior;
    }
    for h in poly.interiors() {
        if point_in_ring(p, &h.0) {
            return Loc::Exterior; // inside a hole
        }
    }
    Loc::Interior
}

fn orient(a: Coord, b: Coord, c: Coord) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

fn on_bbox(a: Coord, b: Coord, p: Coord) -> bool {
    let eps = 1e-9 * (1.0 + a.x.abs() + b.x.abs() + a.y.abs() + b.y.abs());
    p.x <= a.x.max(b.x) + eps
        && p.x >= a.x.min(b.x) - eps
        && p.y <= a.y.max(b.y) + eps
        && p.y >= a.y.min(b.y) - eps
}

/// Is `p` on segment `a–b`? Collinearity uses a scale-relative tolerance so
/// computed intersection points and midpoints land on the edges they came from.
fn point_on_segment(p: Coord, a: Coord, b: Coord) -> bool {
    let o = orient(a, b, p);
    let scale = (b.x - a.x).hypot(b.y - a.y) * (p.x - a.x).hypot(p.y - a.y);
    o.abs() <= 1e-9 * scale.max(1.0) && on_bbox(a, b, p)
}

fn point_on_segment_interior(p: Coord, a: Coord, b: Coord) -> bool {
    p != a && p != b && point_on_segment(p, a, b)
}

fn point_on_ring(p: Coord, ring: &[Coord]) -> bool {
    let v = ring_vertices(ring);
    let n = v.len();
    for i in 0..n {
        if point_on_segment(p, v[i], v[(i + 1) % n]) {
            return true;
        }
    }
    false
}

fn point_in_ring(p: Coord, ring: &[Coord]) -> bool {
    let v = ring_vertices(ring);
    let n = v.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (pi, pj) = (v[i], v[j]);
        if (pi.y > p.y) != (pj.y > p.y)
            && p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// A segment's bounding box, expanded by a scale-relative epsilon so that the
/// box test never excludes a pair that the exact (eps-tolerant) intersection
/// test would accept — the box filter stays a safe superset.
fn padded_box(a: Coord, b: Coord) -> Rect {
    let e = 1e-7 * (1.0 + a.x.abs() + a.y.abs() + b.x.abs() + b.y.abs());
    Rect::new(
        Coord::new(a.x.min(b.x) - e, a.y.min(b.y) - e),
        Coord::new(a.x.max(b.x) + e, a.y.max(b.y) + e),
    )
}

/// Index a segment list by padded bounding box for spatial queries.
fn segment_index(segs: &[Seg]) -> StrTree<usize> {
    StrTree::new(segs.iter().enumerate().map(|(i, s)| (padded_box(s.a, s.b), i)).collect())
}

/// Node two segment sets: collect every pairwise intersection point, but only
/// test pairs whose bounding boxes overlap (others cannot intersect, so this is
/// identical to the `O(n·m)` double loop — just without the misses). Indexes the
/// larger set and queries the smaller, so the cost is `O((n+m)·log)`.
fn node_segments(sa: &[Seg], sb: &[Seg], out: &mut Vec<Coord>) {
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

/// Proper + collinear intersection points of two segments.
fn collect_intersections(a1: Coord, a2: Coord, b1: Coord, b2: Coord, out: &mut Vec<Coord>) {
    let d1 = orient(b1, b2, a1);
    let d2 = orient(b1, b2, a2);
    let d3 = orient(a1, a2, b1);
    let d4 = orient(a1, a2, b2);
    if ((d1 > 0.0) != (d2 > 0.0)) && ((d3 > 0.0) != (d4 > 0.0)) && d1 != 0.0 && d2 != 0.0 {
        // proper crossing — solve for the point
        let denom = (a1.x - a2.x) * (b1.y - b2.y) - (a1.y - a2.y) * (b1.x - b2.x);
        if denom != 0.0 {
            let t = ((a1.x - b1.x) * (b1.y - b2.y) - (a1.y - b1.y) * (b1.x - b2.x)) / denom;
            out.push(Coord::new(a1.x + t * (a2.x - a1.x), a1.y + t * (a2.y - a1.y)));
        }
    }
    // collinear / touching endpoints
    for (p, q1, q2) in [(a1, b1, b2), (a2, b1, b2), (b1, a1, a2), (b2, a1, a2)] {
        if point_on_segment(p, q1, q2) {
            out.push(p);
        }
    }
}

fn dedup_coords(pts: &mut Vec<Coord>) {
    pts.retain(|c| c.is_finite());
    pts.sort_by(|a, b| {
        a.x.partial_cmp(&b.x)
            .unwrap()
            .then(a.y.partial_cmp(&b.y).unwrap())
    });
    pts.dedup();
}

/// 1-D contributions: for each edge of `from`, split it at all test points lying
/// on it, then classify each sub-segment's midpoint in both geometries.
fn contribute_edges(
    from: &Decomposed,
    other: &Decomposed,
    pts: &[Coord],
    im: &mut IntersectionMatrix,
    swapped: bool,
) {
    for seg in &from.segments {
        let (a, b) = (seg.a, seg.b);
        // split points on this edge, ordered along it
        let mut on: Vec<Coord> = vec![a, b];
        for &p in pts {
            if point_on_segment_interior(p, a, b) {
                on.push(p);
            }
        }
        on.sort_by(|p, q| {
            let dp = (p.x - a.x).powi(2) + (p.y - a.y).powi(2);
            let dq = (q.x - a.x).powi(2) + (q.y - a.y).powi(2);
            dp.partial_cmp(&dq).unwrap()
        });
        on.dedup();
        for w in on.windows(2) {
            let mid = Coord::new((w[0].x + w[1].x) / 2.0, (w[0].y + w[1].y) / 2.0);
            // the sub-edge interior is in `from`'s own Interior/Boundary by
            // construction (tagged) — don't re-derive it with float geometry.
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
        // Probe each areal interior at many points (vertices, edge-midpoints,
        // centroids), and detect boundary crossings. A proper ring crossing means
        // each interior straddles the other's boundary.
        let a_samples: Vec<Coord> = ga.polygons.iter().flat_map(interior_samples).collect();
        let b_samples: Vec<Coord> = gb.polygons.iter().flat_map(interior_samples).collect();
        let cross = rings_cross(ga, gb);
        let a_in = a_samples.iter().any(|&p| locate_areal(p, &gb.polygons) == Loc::Interior);
        let a_out = a_samples.iter().any(|&p| locate_areal(p, &gb.polygons) == Loc::Exterior);
        let b_in = b_samples.iter().any(|&p| locate_areal(p, &ga.polygons) == Loc::Interior);
        let b_out = b_samples.iter().any(|&p| locate_areal(p, &ga.polygons) == Loc::Exterior);
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
        // an areal operand's 2-D interior is always (partly) outside a lower-D one
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

    // boundary of `from` outside other (rings for areal, endpoints for linear)
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

    // 0-D / 1-D interior of `from` outside other (points & lines)
    if from.dim == 0 && from.points.iter().any(|&p| other.locate(p) == Loc::Exterior) {
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
        // ring edges represent the boundary — sample vertices and edge midpoints
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
