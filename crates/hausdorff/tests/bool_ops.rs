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

#[test]
fn disjoint() {
    check_rects(&[(0.0, 0.0, 1.0, 1.0)], &[(5.0, 5.0, 6.0, 6.0)]);
}

#[test]
fn nested() {
    check_rects(&[(0.0, 0.0, 10.0, 10.0)], &[(2.0, 2.0, 4.0, 4.0)]);
}

#[test]
fn crossing_squares() {
    check_rects(&[(0.0, 0.0, 4.0, 4.0)], &[(2.0, 2.0, 6.0, 6.0)]);
}

#[test]
fn identical() {
    check_rects(&[(0.0, 0.0, 3.0, 3.0)], &[(0.0, 0.0, 3.0, 3.0)]);
}

#[test]
fn shared_edge() {
    check_rects(&[(0.0, 0.0, 2.0, 2.0)], &[(2.0, 0.0, 4.0, 2.0)]);
}

#[test]
fn partial_edge_overlap() {
    check_rects(&[(0.0, 0.0, 4.0, 2.0)], &[(2.0, 0.0, 6.0, 2.0)]);
}

#[test]
fn t_vertex_touching() {
    check_rects(&[(0.0, 0.0, 2.0, 2.0)], &[(2.0, 0.0, 4.0, 1.0)]);
}

#[test]
fn t_vertex_grid_tiles() {
    // left rectangle vs two stacked right rectangles (all touching along x = 4).
    check_rects(
        &[(0.0, 0.0, 4.0, 2.0)],
        &[(4.0, 0.0, 6.0, 1.0), (4.0, 1.0, 6.0, 2.0)],
    );
}

#[test]
fn t_vertex_overlap() {
    check_rects(&[(0.0, 0.0, 4.0, 4.0)], &[(2.0, 1.0, 6.0, 3.0)]);
}

#[test]
fn multipolygon_inputs() {
    check_rects(
        &[(0.0, 0.0, 2.0, 2.0), (4.0, 4.0, 6.0, 6.0)],
        &[(1.0, 1.0, 5.0, 5.0)],
    );
}

// ---- concave / holed fixtures: golden intersection areas captured from the
// validated implementation (HOLE and STAIR additionally cross-checked by hand). ----

#[test]
fn t_vertex_staircase_union() {
    let big = sq(0.0, 0.0, 6.0, 6.0);
    let comb = vec![
        (3.0, 0.0),
        (6.0, 0.0),
        (6.0, 2.0),
        (4.0, 2.0),
        (4.0, 4.0),
        (6.0, 4.0),
        (6.0, 6.0),
        (3.0, 6.0),
    ];
    check_ops(&hmp(&[&[&big]]), &hmp(&[&[&comb]]), 14.0);
}

#[test]
fn concave_l_shapes() {
    let l1 = [
        (0.0, 0.0),
        (6.0, 0.0),
        (6.0, 2.0),
        (2.0, 2.0),
        (2.0, 6.0),
        (0.0, 6.0),
    ];
    let l2 = [
        (1.0, 1.0),
        (5.0, 1.0),
        (5.0, 5.0),
        (3.0, 5.0),
        (3.0, 3.0),
        (1.0, 3.0),
    ];
    check_ops(&hmp(&[&[&l1]]), &hmp(&[&[&l2]]), 5.0);
}

#[test]
fn star_into_square() {
    let star = [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (3.0, 3.0), (0.0, 6.0)];
    check_ops(&hmp(&[&[&star]]), &hmp(&[&[&sq(1.0, 1.0, 5.0, 5.0)]]), 12.0);
}

#[test]
fn subject_with_hole() {
    let outer = sq(0.0, 0.0, 10.0, 10.0);
    let hole = sq(3.0, 3.0, 7.0, 7.0);
    check_ops(
        &hmp(&[&[&outer, &hole]]),
        &hmp(&[&[&sq(5.0, 5.0, 12.0, 12.0)]]),
        21.0,
    );
}

/// Small deterministic LCG so the stress battery is reproducible without rand.
struct Lcg(u64);
impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
    fn coord(&mut self, n: i64) -> f64 {
        (self.next_u32() as i64 % n) as f64
    }
}

/// A random axis-aligned rectangle on an integer grid.
fn rand_rect(rng: &mut Lcg, grid: i64) -> (f64, f64, f64, f64) {
    let x0 = rng.coord(grid);
    let y0 = rng.coord(grid);
    let w = 1.0 + rng.coord(grid / 2);
    let h = 1.0 + rng.coord(grid / 2);
    (x0, y0, x0 + w, y0 + h)
}

#[test]
fn stress_random_rectangles() {
    let mut rng = Lcg(0x1234_5678_9abc_def0);
    for _ in 0..400 {
        let a = rand_rect(&mut rng, 12);
        let b = rand_rect(&mut rng, 12);
        check_rects(&[a], &[b]);
    }
}

#[test]
fn stress_random_multirects() {
    let mut rng = Lcg(0x0fed_cba9_8765_4321);
    for _ in 0..200 {
        let a1 = rand_rect(&mut rng, 10);
        let r2 = rand_rect(&mut rng, 10);
        let a2 = (r2.0 + 20.0, r2.1, r2.2 + 20.0, r2.3); // shifted clear of a1 (x < 10)
        let (bx, by) = (rng.coord(6), rng.coord(8));
        let b1 = (bx, by, bx + 18.0 + rng.coord(14), by + 2.0 + rng.coord(8));
        check_rects(&[a1, a2], &[b1]);
    }
}

fn rand_tri(rng: &mut Lcg, grid: i64) -> Vec<(f64, f64)> {
    vec![
        (rng.coord(grid), rng.coord(grid)),
        (rng.coord(grid), rng.coord(grid)),
        (rng.coord(grid), rng.coord(grid)),
    ]
}

fn tri_area2(t: &[(f64, f64)]) -> f64 {
    let (a, b, c) = (t[0], t[1], t[2]);
    ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs()
}

#[test]
fn stress_random_triangles() {
    let mut rng = Lcg(0xdead_beef_cafe_0001);
    let mut tested = 0;
    while tested < 300 {
        let a = rand_tri(&mut rng, 16);
        let b = rand_tri(&mut rng, 16);
        if tri_area2(&a) < 2.0 || tri_area2(&b) < 2.0 {
            continue;
        }
        // Sutherland–Hodgman gives the independent intersection-area oracle.
        let expected = convex_intersection_area(&a, &b);
        check_ops(&hmp(&[&[&a]]), &hmp(&[&[&b]]), expected);
        tested += 1;
    }
}

mod unary_union_conformance {
    use super::*;
    use hausdorff::space::r2::bool_ops::unary_union;

    /// Independent oracle for the cascaded `unary_union`: a linear left-fold of
    /// pairwise unions (a different code path through the same clipper).
    fn linear_fold(polys: &[Polygon]) -> MultiPolygon {
        polys.iter().fold(MultiPolygon::default(), |acc, p| {
            acc.union(&MultiPolygon::new(vec![p.clone()]))
        })
    }

    fn check_unary(polys: &[Vec<(f64, f64)>]) -> MultiPolygon {
        let ours: Vec<Polygon> = polys.iter().map(|p| hpoly(&[p])).collect();
        let our_union = unary_union(&ours);
        let folded = linear_fold(&ours);
        assert!(
            (ha(&our_union) - ha(&folded)).abs() < EPS,
            "unary_union area {} vs linear fold {}",
            ha(&our_union),
            ha(&folded)
        );
        our_union
    }

    #[test]
    fn overlapping_grid_of_squares() {
        let mut polys = Vec::new();
        for i in 0..4 {
            for j in 0..4 {
                let (x, y) = (i as f64, j as f64);
                polys.push(sq(x, y, x + 1.5, y + 1.5));
            }
        }
        let out = check_unary(&polys);
        assert_eq!(out.len(), 1);
        assert!((ha(&out) - 4.5 * 4.5).abs() < EPS);
    }

    #[test]
    fn touching_grid_of_squares() {
        let mut polys = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                let (x, y) = (i as f64, j as f64);
                polys.push(sq(x, y, x + 1.0, y + 1.0));
            }
        }
        let out = check_unary(&polys);
        assert_eq!(out.len(), 1);
        assert!((ha(&out) - 9.0).abs() < EPS);
    }

    #[test]
    fn disjoint_squares() {
        let polys: Vec<_> = (0..6)
            .map(|i| {
                let x = 3.0 * i as f64;
                sq(x, 0.0, x + 1.0, 1.0)
            })
            .collect();
        let out = check_unary(&polys);
        assert_eq!(out.len(), 6);
        assert!((ha(&out) - 6.0).abs() < EPS);
    }

    #[test]
    fn nested_rings() {
        let donut = |lo: f64, hi: f64| {
            let outer = sq(lo, lo, hi, hi);
            let inner: Vec<(f64, f64)> = sq(lo + 1.0, lo + 1.0, hi - 1.0, hi - 1.0)
                .into_iter()
                .rev()
                .collect();
            hpoly(&[&outer, &inner])
        };
        let polys = vec![
            donut(0.0, 12.0),
            donut(2.0, 10.0),
            hpoly(&[&sq(5.0, 5.0, 7.0, 7.0)]),
        ];
        let our_union = unary_union(&polys);
        assert!(
            (ha(&our_union) - ha(&linear_fold(&polys))).abs() < EPS,
            "nested rings area {} vs linear fold {}",
            ha(&our_union),
            ha(&linear_fold(&polys))
        );
    }

    #[test]
    fn multipolygon_operands() {
        let polys = vec![
            hpoly(&[&sq(0.0, 0.0, 2.0, 2.0)]),
            hpoly(&[&sq(1.0, 1.0, 3.0, 3.0)]),
            hpoly(&[&sq(10.0, 0.0, 11.0, 1.0)]),
        ];
        let our_union = unary_union(&polys);
        assert!((ha(&our_union) - ha(&linear_fold(&polys))).abs() < EPS);
    }

    #[test]
    fn area_is_conserved_on_known_fixture() {
        let a = sq(0.0, 0.0, 4.0, 4.0);
        let b = sq(2.0, 2.0, 6.0, 6.0);
        let c = sq(10.0, 10.0, 12.0, 12.0);
        let polys = vec![hpoly(&[&a]), hpoly(&[&b]), hpoly(&[&c])];
        let out = unary_union(&polys);
        let expected = 16.0 + 16.0 - 4.0 + 4.0;
        assert!(
            (ha(&out) - expected).abs() < EPS,
            "area {} vs expected {expected}",
            ha(&out)
        );
    }

    #[test]
    fn matches_linear_pairwise_fold() {
        let mut rng = Lcg(0x5151_5151_5151_5151);
        for _ in 0..20 {
            let polys: Vec<Polygon> = (0..8)
                .map(|_| hpoly(&[&sq2(&rand_rect(&mut rng, 10))]))
                .collect();
            let cascaded = unary_union(&polys);
            assert!(
                (ha(&cascaded) - ha(&linear_fold(&polys))).abs() < EPS,
                "cascaded {} vs linear {}",
                ha(&cascaded),
                ha(&linear_fold(&polys))
            );
        }
    }

    fn sq2(r: &(f64, f64, f64, f64)) -> Vec<(f64, f64)> {
        sq(r.0, r.1, r.2, r.3)
    }
}
