//! Tests for the STR-tree spatial index: range queries and nearest-neighbour
//! must agree with brute force over random data.

use hausdorff::space::r2::{Coord, Rect, StrTree};

struct Rng(u64);
impl Rng {
    fn f(&mut self, half: f64) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let r = ((self.0 >> 33) & 0xFFFF_FFFF) as f64 / 4_294_967_295.0;
        (r * 2.0 - 1.0) * half
    }
}

fn rect(x: f64, y: f64, w: f64, h: f64) -> Rect {
    Rect::new(Coord::new(x, y), Coord::new(x + w, y + h))
}

fn overlaps(a: Rect, b: Rect) -> bool {
    a.min().x <= b.max().x
        && b.min().x <= a.max().x
        && a.min().y <= b.max().y
        && b.min().y <= a.max().y
}

fn dist_sq(p: Coord, q: Coord) -> f64 {
    (p.x - q.x).powi(2) + (p.y - q.y).powi(2)
}

#[test]
fn empty_and_singleton() {
    let empty: StrTree<u32> = StrTree::new(vec![]);
    assert!(empty.is_empty());
    assert_eq!(empty.len(), 0);
    assert!(empty.query(rect(0.0, 0.0, 1.0, 1.0)).is_empty());
    assert!(empty.nearest(Coord::new(0.0, 0.0)).is_none());

    let one = StrTree::new(vec![(rect(0.0, 0.0, 2.0, 2.0), 7u32)]);
    assert_eq!(one.len(), 1);
    assert_eq!(one.query(rect(1.0, 1.0, 5.0, 5.0)), vec![&7]);
    assert!(one.query(rect(9.0, 9.0, 1.0, 1.0)).is_empty());
    assert_eq!(one.nearest(Coord::new(10.0, 10.0)), Some(&7));
}

#[test]
fn query_matches_brute_force() {
    let mut rng = Rng(0xA11C_E202_4EC6);
    let items: Vec<(Rect, usize)> = (0..1500)
        .map(|i| {
            (
                rect(
                    rng.f(500.0),
                    rng.f(500.0),
                    1.0 + rng.f(10.0).abs(),
                    1.0 + rng.f(10.0).abs(),
                ),
                i,
            )
        })
        .collect();
    let tree = StrTree::new(items.clone());
    assert_eq!(tree.len(), items.len());

    for _ in 0..300 {
        let q = rect(
            rng.f(500.0),
            rng.f(500.0),
            rng.f(60.0).abs(),
            rng.f(60.0).abs(),
        );
        let mut got: Vec<usize> = tree.query(q).into_iter().copied().collect();
        got.sort_unstable();
        let mut want: Vec<usize> = items
            .iter()
            .filter(|(r, _)| overlaps(*r, q))
            .map(|(_, i)| *i)
            .collect();
        want.sort_unstable();
        assert_eq!(got, want, "query mismatch for {q:?}");
    }
}

#[test]
fn nearest_matches_brute_force_on_points() {
    let mut rng = Rng(0xBEEF_F00D);
    // degenerate (point) boxes → nearest is exact 1-NN
    let pts: Vec<(Rect, usize)> = (0..800)
        .map(|i| {
            let p = Coord::new(rng.f(1000.0), rng.f(1000.0));
            (Rect::new(p, p), i)
        })
        .collect();
    let tree = StrTree::new(pts.clone());

    for _ in 0..300 {
        let q = Coord::new(rng.f(1000.0), rng.f(1000.0));
        let got = *tree.nearest(q).expect("non-empty");
        let got_d = dist_sq(q, pts[got].0.min());
        let brute_d = pts
            .iter()
            .map(|(r, _)| dist_sq(q, r.min()))
            .fold(f64::INFINITY, f64::min);
        assert!(
            (got_d - brute_d).abs() < 1e-9,
            "nearest {got_d} vs brute {brute_d}"
        );
    }
}

#[test]
fn query_each_visits_the_same_items() {
    let items: Vec<(Rect, usize)> = (0..200)
        .map(|i| (rect(i as f64, 0.0, 1.0, 1.0), i))
        .collect();
    let tree = StrTree::new(items);
    let q = rect(10.5, -1.0, 5.0, 3.0); // spans x in [10.5, 15.5]
    let mut via_each: Vec<usize> = Vec::new();
    tree.query_each(q, |_, &i| via_each.push(i));
    via_each.sort_unstable();
    let mut via_query: Vec<usize> = tree.query(q).into_iter().copied().collect();
    via_query.sort_unstable();
    assert_eq!(via_each, via_query);
    assert_eq!(via_query, vec![10, 11, 12, 13, 14, 15]);
}
