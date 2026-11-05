//! Delaunay and Voronoi, against the definition and against Euler's formula.
//!
//! Three independent oracles, none of which needs a reference implementation:
//!
//! 1. **The Delaunay property itself** — no point strictly inside any triangle's
//!    circumcircle. That is the definition, so checking it is not a proxy for
//!    correctness, it *is* correctness. `Triangulation::is_delaunay` exposes it.
//! 2. **Euler's formula.** For a triangulation of `n` points with `h` on the
//!    convex hull, the triangle count is exactly `2n − 2 − h`. A dropped or
//!    duplicated triangle changes that integer, so it catches structural errors
//!    the circumcircle test would not — an empty triangulation trivially satisfies
//!    the Delaunay property.
//! 3. **Area conservation.** The triangles must tile the convex hull exactly, so
//!    their areas sum to the hull's. This catches overlaps and gaps, which Euler's
//!    count alone can miss if two errors cancel.
//!
//! Together they pin the result: the right *number* of triangles, covering the
//! right *area*, satisfying the right *property*. For Voronoi the oracle is the
//! defining property instead — every point in a cell is nearer that cell's site
//! than any other.
//!
//! Split by sub-topic because the single file had grown past the 400-line house
//! rule. The shared fixtures and the soundness battery live here; the three
//! children are the Delaunay properties, the adversarial inputs that the
//! cocircular defect was found with, and the Voronoi dual.

use crate::space::r2::{
    Area, Coord, CoordsIter, Geometry, LineString, MultiPoint, Point, Polygon, Triangulate,
    delaunay, hull_area, nearest_site, voronoi,
};

fn pts(v: &[(f64, f64)]) -> Vec<Coord> {
    v.iter().map(|&(x, y)| Coord::new(x, y)).collect()
}

/// A deterministic scatter, so a failure is reproducible.
fn scatter(n: usize, scale: f64, seed: u64) -> Vec<Coord> {
    let mut s = seed | 1;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let x = ((s >> 11) as f64 / (1u64 << 53) as f64) * scale;
        s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let y = ((s >> 11) as f64 / (1u64 << 53) as f64) * scale;
        out.push(Coord::new(x, y));
    }
    out
}

/// The number of input points on the convex hull, needed for Euler's count.
///
/// Counted off the hull ring directly rather than by rebuilding a `Polygon` from
/// it: `Polygon::new` re-closes a ring, so a round trip through it reports one
/// vertex more than the hull has. The first version of this helper did that and
/// claimed a triangle had four hull points.
fn hull_count(points: &[Coord]) -> usize {
    let ring = crate::space::r2::ConvexHull::convex_hull(&MultiPoint::new(
        points.iter().map(|&c| Point(c)).collect(),
    ))
    .exterior()
    .0
    .clone();
    // Euler's formula counts points **on the hull boundary**, not hull *vertices*.
    // The two differ exactly where a hull edge has input points along it, which is
    // every lattice: a 7x7 grid has 4 hull vertices and 24 boundary points, and
    // using 4 predicts 92 triangles where 72 is right. Counting vertices was the
    // first version of this helper and it accused a correct triangulation of
    // losing twenty triangles.
    // Tolerances scale with the cloud's **extent**, not with its distance from the
    // origin. Those are different numbers and confusing them breaks at both ends:
    //
    //  - a fixed absolute epsilon is far too loose at an extent of 1e-9, marking
    //    interior points as lying on the boundary; and
    //  - an epsilon taken from the coordinate *magnitude* is far too loose for a
    //    cloud one unit wide sitting at 1e9, where it becomes 1.0 and swallows the
    //    whole point set.
    //
    // Both were in earlier versions of this helper and both accused a correct
    // triangulation of being wrong. The extent is the only scale a geometric
    // tolerance can sensibly be relative to.
    let (mut lo_x, mut hi_x, mut lo_y, mut hi_y) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for p in points {
        lo_x = lo_x.min(p.x);
        hi_x = hi_x.max(p.x);
        lo_y = lo_y.min(p.y);
        hi_y = hi_y.max(p.y);
    }
    let span = (hi_x - lo_x).max(hi_y - lo_y).max(f64::MIN_POSITIVE);
    let eps = span * 1e-9;
    let mut on_boundary = 0;
    for &p in points {
        let mut touches = false;
        for w in ring.windows(2) {
            let (a, b) = (w[0], w[1]);
            let cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
            let edge = (b.x - a.x).abs().max((b.y - a.y).abs()).max(span);
            if cross.abs() <= eps * edge
                && p.x >= a.x.min(b.x) - eps
                && p.x <= a.x.max(b.x) + eps
                && p.y >= a.y.min(b.y) - eps
                && p.y <= a.y.max(b.y) + eps
            {
                touches = true;
                break;
            }
        }
        if touches {
            on_boundary += 1;
        }
    }
    on_boundary
}

/// Sum of triangle areas and the hull area, both computed on coordinates
/// translated to the centroid.
///
/// Translating matters for the same reason it matters inside the triangulator: at
/// a magnitude of 1e9 a triangle's area is a difference of products near 1e18, and
/// comparing two such sums measures the *comparison's* rounding rather than the
/// triangulation's correctness.
fn areas_translated(t: &crate::space::r2::Triangulation) -> (f64, f64) {
    let n = t.points.len() as f64;
    let shift = Coord::new(
        t.points.iter().map(|p| p.x).sum::<f64>() / n,
        t.points.iter().map(|p| p.y).sum::<f64>() / n,
    );
    let at = |c: Coord| Coord::new(c.x - shift.x, c.y - shift.y);
    let sum: f64 = t
        .triangles
        .iter()
        .map(|&[a, b, c]| {
            let (pa, pb, pc) = (at(t.points[a]), at(t.points[b]), at(t.points[c]));
            ((pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x)).abs() / 2.0
        })
        .sum();
    let moved: Vec<Coord> = t.points.iter().map(|&c| at(c)).collect();
    (sum, hull_area(&moved))
}

/// Assert the three structural invariants at once: Euler's triangle count, exact
/// tiling of the convex hull, and the empty-circumcircle property.
#[track_caller]
fn assert_sound(name: &str, points: &[Coord]) -> crate::space::r2::Triangulation {
    let t = delaunay(points);
    let n = t.points.len();
    if n < 3 {
        assert!(
            t.triangles.is_empty(),
            "{name}: fewer than 3 points cannot triangulate"
        );
        return t;
    }
    let h = hull_count(&t.points);
    assert_eq!(
        t.triangles.len(),
        2 * n - 2 - h,
        "{name}: Euler's count for {n} points with {h} on the hull boundary"
    );
    let (sum, hull) = areas_translated(&t);
    assert!(
        (sum - hull).abs() <= 1e-9 * hull.max(f64::MIN_POSITIVE),
        "{name}: triangles total {sum} but the hull is {hull} -- an overlap or a gap"
    );
    assert!(
        t.is_delaunay(),
        "{name}: the empty-circumcircle property must hold"
    );
    t
}

mod adversarial;
mod delaunay;
mod voronoi;
