//! Self-contained set-algebra conformance for the native Martinez–Rueda BooleanOps.
//!
//! Formerly a differential test against `geo`. The reference is now provided by
//! dependency-free independent oracles:
//!   * a closed-form axis-aligned rectangle-overlap formula (rectangle inputs),
//!   * a Sutherland–Hodgman convex clipper (triangle inputs),
//!   * golden constants captured from the validated implementation (the handful of
//!     concave/holed fixtures whose intersection has no simple closed form),
//! plus the full set-algebra invariants (inclusion–exclusion, difference,
//! symmetric difference, commutativity, idempotence, monotone bounds), which tie
//! the four operations to one another and to `unsigned_area`. No geometry crate.

use hausdorff::space::r2::bool_ops::BooleanOps;
use hausdorff::{Area, LineString, MultiPolygon, Polygon};

fn hpoly(rings: &[&[(f64, f64)]]) -> Polygon {
    let ext = LineString::from(rings[0].to_vec());
    let holes = rings[1..]
        .iter()
        .map(|r| LineString::from(r.to_vec()))
        .collect();
    Polygon::new(ext, holes)
}

fn hmp(polys: &[&[&[(f64, f64)]]]) -> MultiPolygon {
    MultiPolygon::new(polys.iter().map(|p| hpoly(p)).collect())
}

fn ha(mp: &MultiPolygon) -> f64 {
    mp.unsigned_area()
}

const EPS: f64 = 1e-6;

// ---- independent, dependency-free reference oracles ----

fn shoelace(poly: &[(f64, f64)]) -> f64 {
    if poly.len() < 3 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..poly.len() {
        let (x0, y0) = poly[i];
        let (x1, y1) = poly[(i + 1) % poly.len()];
        s += x0 * y1 - x1 * y0;
    }
    s.abs() * 0.5
}

/// Closed-form overlap area of two axis-aligned rectangles `(x0, y0, x1, y1)`
/// with `x0 < x1`, `y0 < y1`.
fn rect_overlap_area(a: (f64, f64, f64, f64), b: (f64, f64, f64, f64)) -> f64 {
    let ix = (a.2.min(b.2) - a.0.max(b.0)).max(0.0);
    let iy = (a.3.min(b.3) - a.1.max(b.1)).max(0.0);
    ix * iy
}

/// Area of the intersection of two convex polygons via Sutherland–Hodgman
/// (an algorithm wholly independent of the crate's Martinez–Rueda clipper).
/// `clip` is reoriented counter-clockwise; `subject` may be any convex polygon.
fn convex_intersection_area(subject: &[(f64, f64)], clip: &[(f64, f64)]) -> f64 {
    fn signed_area(p: &[(f64, f64)]) -> f64 {
        let mut s = 0.0;
        for i in 0..p.len() {
            let (x0, y0) = p[i];
            let (x1, y1) = p[(i + 1) % p.len()];
            s += x0 * y1 - x1 * y0;
        }
        s * 0.5
    }
    let mut clip = clip.to_vec();
    if signed_area(&clip) < 0.0 {
        clip.reverse();
    }
    let mut out: Vec<(f64, f64)> = subject.to_vec();
    let n = clip.len();
    for i in 0..n {
        if out.is_empty() {
            break;
        }
        let a = clip[i];
        let b = clip[(i + 1) % n];
        // For CCW `clip`, "inside" is left of (or on) the directed edge a->b.
        let inside =
            |p: (f64, f64)| (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0) >= -1e-12;
        let isect = |p: (f64, f64), q: (f64, f64)| {
            let r = (q.0 - p.0, q.1 - p.1);
            let d = (b.0 - a.0, b.1 - a.1);
            let denom = r.0 * d.1 - r.1 * d.0;
            let t = ((a.0 - p.0) * d.1 - (a.1 - p.1) * d.0) / denom;
            (p.0 + t * r.0, p.1 + t * r.1)
        };
        let input = std::mem::take(&mut out);
        let m = input.len();
        for j in 0..m {
            let cur = input[j];
            let prev = input[(j + m - 1) % m];
            let (cur_in, prev_in) = (inside(cur), inside(prev));
            if cur_in {
                if !prev_in {
                    out.push(isect(prev, cur));
                }
                out.push(cur);
            } else if prev_in {
                out.push(isect(prev, cur));
            }
        }
    }
    shoelace(&out)
}

/// Assert the four boolean operations agree with an independently-computed
/// intersection area and satisfy the set-algebra invariants among themselves.
fn check_ops(a: &MultiPolygon, b: &MultiPolygon, expected_inter: f64) {
    let inter = ha(&a.intersection(b));
    let uni = ha(&a.union(b));
    let diff = ha(&a.difference(b));
    let xor = ha(&a.xor(b));
    let (aa, ba) = (ha(a), ha(b));

    // Absolute anchor from an independent oracle.
    assert!(
        (inter - expected_inter).abs() < EPS,
        "intersection {inter} vs expected {expected_inter}"
    );
    // Inclusion–exclusion and the derived identities.
    assert!((uni - (aa + ba - inter)).abs() < EPS, "A∪B = A+B-A∩B");
    assert!((diff - (aa - inter)).abs() < EPS, "A\\B = A-A∩B");
    assert!((xor - (uni - inter)).abs() < EPS, "A⊕B = A∪B-A∩B");
    // Monotone bounds.
    assert!(
        inter >= -EPS && inter <= aa.min(ba) + EPS,
        "0 ≤ A∩B ≤ min(A,B): {inter}"
    );
    assert!(
        uni >= aa.max(ba) - EPS && uni <= aa + ba + EPS,
        "max(A,B) ≤ A∪B ≤ A+B: {uni}"
    );
    // Commutativity (area level).
    assert!((ha(&b.intersection(a)) - inter).abs() < EPS, "A∩B = B∩A");
    assert!((ha(&b.union(a)) - uni).abs() < EPS, "A∪B = B∪A");
    // Idempotence.
    assert!((ha(&a.intersection(a)) - aa).abs() < EPS, "A∩A = A");
    assert!((ha(&a.union(a)) - aa).abs() < EPS, "A∪A = A");
}

fn sq(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
}

/// Build multipolygons from axis-aligned rectangles and verify the boolean ops.
/// Rectangles within each operand MUST be mutually disjoint, so the intersection
/// area is exactly the sum of pairwise overlaps (the independent oracle).
fn check_rects(a_rects: &[(f64, f64, f64, f64)], b_rects: &[(f64, f64, f64, f64)]) {
    let to_mp = |rs: &[(f64, f64, f64, f64)]| {
        MultiPolygon::new(
            rs.iter()
                .map(|&(x0, y0, x1, y1)| hpoly(&[&sq(x0, y0, x1, y1)]))
                .collect(),
        )
    };
    let expected_inter: f64 = a_rects
        .iter()
        .flat_map(|&a| b_rects.iter().map(move |&b| rect_overlap_area(a, b)))
        .sum();
    check_ops(&to_mp(a_rects), &to_mp(b_rects), expected_inter);
}

#[path = "bool_ops/basic.rs"]
mod basic;

#[path = "bool_ops/isochrone_bands.rs"]
mod isochrone_bands;

#[path = "bool_ops/pinched.rs"]
mod pinched;

#[path = "bool_ops/random_rects.rs"]
mod random_rects;

#[path = "bool_ops/random_triangles.rs"]
mod random_triangles;

#[path = "bool_ops/unary_union.rs"]
mod unary_union;
