//! Coordinate- and shape-level analysis used before and during the relate: a
//! geometry's topological [`dimension`], its raw bounding box ([`geom_bbox`] via
//! the non-allocating [`fold_coords`]), the strict box-disjointness test, and
//! the fast [`disjoint_matrix`] / [`boundary_dim`] shortcuts.

use super::super::{Coord, Geometry, Polygon};
use super::{B, E, I, IntersectionMatrix};

/// Topological dimension of a geometry (max over components): 0/1/2, or -1 empty.
pub(super) fn dimension(g: &Geometry) -> i8 {
    match g {
        Geometry::Point(_) | Geometry::MultiPoint(_) => 0,
        Geometry::Line(_) | Geometry::LineString(_) | Geometry::MultiLineString(_) => 1,
        Geometry::Polygon(_)
        | Geometry::MultiPolygon(_)
        | Geometry::Rect(_)
        | Geometry::Triangle(_) => 2,
        Geometry::GeometryCollection(gs) => gs.iter().map(dimension).max().unwrap_or(-1),
    }
}

/// Bounding box `(min, max)` of a geometry from its raw coordinates (no
/// decomposition), or `None` if it has none. Folds the raw coordinate slices
/// directly — unlike `coords_iter()`, which materializes a `Vec` — so the
/// disjoint-box fast path in [`relate`](super::relate) allocates nothing.
pub(super) fn geom_bbox(g: &Geometry) -> Option<(Coord, Coord)> {
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
pub(super) fn boxes_strictly_disjoint(a: (Coord, Coord), b: (Coord, Coord)) -> bool {
    a.1.x < b.0.x || b.1.x < a.0.x || a.1.y < b.0.y || b.1.y < a.0.y
}

/// The DE-9IM matrix of two disjoint geometries — each operand's interior and
/// boundary lie wholly in the other's exterior. `None` when an operand's boundary
/// dimension is ambiguous (mod-2 line collections), signalling that the exact
/// decomposition path should run instead.
pub(super) fn disjoint_matrix(a: &Geometry, b: &Geometry) -> Option<IntersectionMatrix> {
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
