//! Spatial join (STR-tree build-and-probe) must produce exactly the brute-force
//! O(n·m) all-pairs result: the broad phase must match an all-pairs box test,
//! and the predicate-filtered join must match an all-pairs predicate scan.

use hausdorff::space::r2::{spatial_join, spatial_join_candidates};
use hausdorff::{BoundingRect, Geometry, Intersects, LineString, Polygon, Rect};

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn coord(&mut self, n: u64) -> f64 {
        (self.next() % n) as f64 - n as f64 / 2.0
    }
}

fn sq(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    Geometry::Polygon(Polygon::new(
        LineString::from(vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]),
        vec![],
    ))
}

fn rects_overlap(a: Rect, b: Rect) -> bool {
    a.min().x <= b.max().x
        && b.min().x <= a.max().x
        && a.min().y <= b.max().y
        && b.min().y <= a.max().y
}

fn column(rng: &mut Lcg, n: usize) -> Vec<Geometry> {
    (0..n)
        .map(|_| {
            let (x, y) = (rng.coord(60), rng.coord(60));
            let (w, h) = (1.0 + rng.coord(10).abs(), 1.0 + rng.coord(10).abs());
            sq(x, y, x + w, y + h)
        })
        .collect()
}

#[test]
fn candidates_match_brute_bbox_overlap() {
    let mut rng = Lcg(0xABCD_1234_5678_9EF0);
    let left = column(&mut rng, 120);
    let right = column(&mut rng, 150);

    let mut got = spatial_join_candidates(&left, &right);
    got.sort_unstable();

    let mut want = Vec::new();
    for (i, l) in left.iter().enumerate() {
        for (j, r) in right.iter().enumerate() {
            if let (Some(lb), Some(rb)) = (l.bounding_rect(), r.bounding_rect())
                && rects_overlap(lb, rb)
            {
                want.push((i, j));
            }
        }
    }
    want.sort_unstable();
    assert_eq!(got, want);
}

#[test]
fn join_intersects_matches_brute() {
    let mut rng = Lcg(0x0F1E_2D3C_4B5A_6978);
    let left = column(&mut rng, 120);
    let right = column(&mut rng, 150);

    let mut got = spatial_join(&left, &right, |a, b| a.intersects(b));
    got.sort_unstable();

    let mut want = Vec::new();
    for (i, l) in left.iter().enumerate() {
        for (j, r) in right.iter().enumerate() {
            if l.intersects(r) {
                want.push((i, j));
            }
        }
    }
    want.sort_unstable();
    assert_eq!(got, want);
}
