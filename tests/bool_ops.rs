//! Conformance for the native Martinez–Rueda BooleanOps against geo's BooleanOps.
//! Result *areas* must match geo for intersection/union/difference/xor across a
//! battery of configurations (disjoint, nested, crossing, concave, touching,
//! identical, holes, multipolygons), and the set-algebra invariants must hold.

use geo::BooleanOps as _;
use hausdorff::space::r2::bool_ops::BooleanOps;
use hausdorff::{Area, LineString, MultiPolygon, Polygon};

fn hpoly(rings: &[&[(f64, f64)]]) -> Polygon {
    let ext = LineString::from(rings[0].to_vec());
    let holes = rings[1..].iter().map(|r| LineString::from(r.to_vec())).collect();
    Polygon::new(ext, holes)
}
fn gpoly(rings: &[&[(f64, f64)]]) -> geo::Polygon<f64> {
    let ring = |r: &[(f64, f64)]| geo::LineString(r.iter().map(|&(x, y)| geo::Coord { x, y }).collect());
    geo::Polygon::new(ring(rings[0]), rings[1..].iter().map(|r| ring(r)).collect())
}

fn hmp(polys: &[&[&[(f64, f64)]]]) -> MultiPolygon {
    MultiPolygon::new(polys.iter().map(|p| hpoly(p)).collect())
}
fn gmp(polys: &[&[&[(f64, f64)]]]) -> geo::MultiPolygon<f64> {
    geo::MultiPolygon(polys.iter().map(|p| gpoly(p)).collect())
}

fn ha(mp: &MultiPolygon) -> f64 {
    mp.unsigned_area()
}
fn ga(mp: &geo::MultiPolygon<f64>) -> f64 {
    use geo::Area;
    mp.unsigned_area()
}

const EPS: f64 = 1e-6;

fn check(a: &[&[&[(f64, f64)]]], b: &[&[&[(f64, f64)]]]) {
    let (ha_, hb_) = (hmp(a), hmp(b));
    let (ga_, gb_) = (gmp(a), gmp(b));

    let inter = ha(&ha_.intersection(&hb_));
    let uni = ha(&ha_.union(&hb_));
    let diff = ha(&ha_.difference(&hb_));
    let xor = ha(&ha_.xor(&hb_));

    let (gi, gu, gd, gx) = (
        ga(&ga_.intersection(&gb_)),
        ga(&ga_.union(&gb_)),
        ga(&ga_.difference(&gb_)),
        ga(&ga_.xor(&gb_)),
    );

    assert!((inter - gi).abs() < EPS, "intersection {inter} vs geo {gi}");
    assert!((uni - gu).abs() < EPS, "union {uni} vs geo {gu}");
    assert!((diff - gd).abs() < EPS, "difference {diff} vs geo {gd}");
    assert!((xor - gx).abs() < EPS, "xor {xor} vs geo {gx}");

    // set-algebra invariants (independent of geo)
    let area_a = ha(&ha_);
    assert!((uni - (area_a + ha(&hb_) - inter)).abs() < EPS, "A∪B = A+B-A∩B");
    assert!((diff - (area_a - inter)).abs() < EPS, "A\\B = A - A∩B");
    assert!((xor - (uni - inter)).abs() < EPS, "A⊕B = A∪B - A∩B");
}

fn sq(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
}

#[test]
fn disjoint() {
    check(&[&[&sq(0.0, 0.0, 1.0, 1.0)]], &[&[&sq(5.0, 5.0, 6.0, 6.0)]]);
}

#[test]
fn nested() {
    check(&[&[&sq(0.0, 0.0, 10.0, 10.0)]], &[&[&sq(2.0, 2.0, 4.0, 4.0)]]);
}

#[test]
fn crossing_squares() {
    check(&[&[&sq(0.0, 0.0, 4.0, 4.0)]], &[&[&sq(2.0, 2.0, 6.0, 6.0)]]);
}

#[test]
fn identical() {
    check(&[&[&sq(0.0, 0.0, 3.0, 3.0)]], &[&[&sq(0.0, 0.0, 3.0, 3.0)]]);
}

#[test]
fn shared_edge() {
    check(&[&[&sq(0.0, 0.0, 2.0, 2.0)]], &[&[&sq(2.0, 0.0, 4.0, 2.0)]]);
}

#[test]
fn partial_edge_overlap() {
    check(&[&[&sq(0.0, 0.0, 4.0, 2.0)]], &[&[&sq(2.0, 0.0, 6.0, 2.0)]]);
}

// ---- T-vertex configurations (the documented overlay limitation) -----------
// A vertex of one polygon lying on the *interior* of the other's edge. Common on
// integer grids / tiled data. These probe whether the sweep mislabels in/out.

#[test]
fn t_vertex_touching_matches_geo() {
    // B's corner (2,1) lies on the interior of A's right edge (x=2, y∈[0,2]).
    // Interiors are disjoint; they touch along x=2, y∈[0,1]. Union = 4 + 2 = 6.
    let a = sq(0.0, 0.0, 2.0, 2.0);
    let b = vec![(2.0, 0.0), (4.0, 0.0), (4.0, 1.0), (2.0, 1.0)];
    check(&[&[&a]], &[&[&b]]);
}

#[test]
fn t_vertex_grid_tiles_matches_geo() {
    // A big left tile abutting two stacked right tiles: the right tiles' shared
    // corner (4,1) is a T-vertex on the big tile's right edge (x=4, y∈[0,2]).
    let left = sq(0.0, 0.0, 4.0, 2.0);
    let right_lo = sq(4.0, 0.0, 6.0, 1.0);
    let right_hi = sq(4.0, 1.0, 6.0, 2.0);
    check(&[&[&left]], &[&[&right_lo], &[&right_hi]]);
}

#[test]
fn t_vertex_overlap_matches_geo() {
    // Overlapping T-vertex: B's left edge runs x=2, y∈[1,3]; its corners (2,1)
    // and (2,3) lie on the interior of A's edges, and the interiors overlap.
    let a = sq(0.0, 0.0, 4.0, 4.0);
    let b = sq(2.0, 1.0, 6.0, 3.0);
    check(&[&[&a]], &[&[&b]]);
}

#[test]
fn t_vertex_staircase_union_matches_geo() {
    // A comb/staircase whose steps put several vertices (6,2) and (6,4) on the
    // interior of the big square's right edge (x=6, y∈[0,6]) — many T-vertices.
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
    check(&[&[&big]], &[&[&comb]]);
}

#[test]
fn concave_l_shapes() {
    let l1 = [(0.0, 0.0), (6.0, 0.0), (6.0, 2.0), (2.0, 2.0), (2.0, 6.0), (0.0, 6.0)];
    let l2 = [(1.0, 1.0), (5.0, 1.0), (5.0, 5.0), (3.0, 5.0), (3.0, 3.0), (1.0, 3.0)];
    check(&[&[&l1]], &[&[&l2]]);
}

#[test]
fn star_into_square() {
    let star = [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (3.0, 3.0), (0.0, 6.0)];
    check(&[&[&star]], &[&[&sq(1.0, 1.0, 5.0, 5.0)]]);
}

#[test]
fn subject_with_hole() {
    let outer = sq(0.0, 0.0, 10.0, 10.0);
    let hole = sq(3.0, 3.0, 7.0, 7.0);
    check(&[&[&outer, &hole]], &[&[&sq(5.0, 5.0, 12.0, 12.0)]]);
}

#[test]
fn multipolygon_inputs() {
    let a = sq(0.0, 0.0, 2.0, 2.0);
    let b = sq(4.0, 4.0, 6.0, 6.0);
    let c = sq(1.0, 1.0, 5.0, 5.0);
    check(&[&[&a], &[&b]], &[&[&c]]);
}

// ---- randomized stress: native vs geo over many pseudo-random polygon pairs ----

/// Small deterministic LCG so the stress battery is reproducible without rand.
struct Lcg(u64);
impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
    fn coord(&mut self, n: i64) -> f64 {
        (self.next_u32() as i64 % n) as f64
    }
}

/// A random axis-aligned rectangle on an integer grid (keeps geometry valid and
/// exercises collinear/touching cases when rectangles align).
fn rand_rect(rng: &mut Lcg, grid: i64) -> Vec<(f64, f64)> {
    let x0 = rng.coord(grid);
    let y0 = rng.coord(grid);
    let w = 1.0 + rng.coord(grid / 2);
    let h = 1.0 + rng.coord(grid / 2);
    sq(x0, y0, x0 + w, y0 + h)
}

#[test]
fn stress_random_rectangles_match_geo() {
    let mut rng = Lcg(0x1234_5678_9abc_def0);
    for _ in 0..400 {
        let a = rand_rect(&mut rng, 12);
        let b = rand_rect(&mut rng, 12);
        check(&[&[&a]], &[&[&b]]);
    }
}

fn shift_x(r: &[(f64, f64)], dx: f64) -> Vec<(f64, f64)> {
    r.iter().map(|&(x, y)| (x + dx, y)).collect()
}

#[test]
fn stress_random_multirects_match_geo() {
    // Keep the subject's two components in disjoint x-bands so the MultiPolygon
    // stays valid (no self-overlap); the clip rect is wide enough to cross both.
    let mut rng = Lcg(0x0fed_cba9_8765_4321);
    for _ in 0..200 {
        let a1 = rand_rect(&mut rng, 10);
        let a2 = shift_x(&rand_rect(&mut rng, 10), 20.0);
        let (bx, by) = (rng.coord(6), rng.coord(8));
        let b1 = sq(bx, by, bx + 18.0 + rng.coord(14), by + 2.0 + rng.coord(8));
        check(&[&[&a1], &[&a2]], &[&[&b1]]);
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
fn stress_random_triangles_match_geo() {
    // general-position crossings (non-axis-aligned); skip near-degenerate tris.
    let mut rng = Lcg(0xdead_beef_cafe_0001);
    let mut tested = 0;
    while tested < 300 {
        let a = rand_tri(&mut rng, 16);
        let b = rand_tri(&mut rng, 16);
        if tri_area2(&a) < 2.0 || tri_area2(&b) < 2.0 {
            continue;
        }
        check(&[&[&a]], &[&[&b]]);
        tested += 1;
    }
}
