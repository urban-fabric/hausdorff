//! `PreparedGeometry`: indexed repeated predicates must agree with the one-shot
//! (non-prepared) kernels they accelerate.

use hausdorff::space::r2::{PointBatch, PreparedGeometry, euclidean_distance, euclidean_within};
use hausdorff::{Contains, Coord, Geometry, Intersects, LineString, Polygon};

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn coord(&mut self, n: u64) -> f64 {
        (self.next() % n) as f64 - n as f64 / 2.0
    }
}

fn poly(coords: &[(f64, f64)]) -> Geometry {
    Geometry::Polygon(Polygon::new(LineString::from(coords.to_vec()), vec![]))
}
fn sq(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    poly(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
}

#[test]
fn prepared_distance_within_intersects_match_oneshot() {
    let base = sq(0.0, 0.0, 10.0, 10.0);
    let prepared = PreparedGeometry::new(&base);
    let mut rng = Lcg(0x9E37_79B9_7F4A_7C15);
    for _ in 0..500 {
        let (x, y) = (rng.coord(40), rng.coord(40));
        let (w, h) = (1.0 + rng.coord(8).abs(), 1.0 + rng.coord(8).abs());
        let other = sq(x, y, x + w, y + h);

        let pd = prepared.distance(&other);
        let od = euclidean_distance(&base, &other);
        assert!((pd - od).abs() <= 1e-9 * od.max(1.0), "distance prepared={pd} one-shot={od}");

        assert_eq!(prepared.intersects(&other), base.intersects(&other), "intersects");

        for d in [0.0, 2.0, 5.0, 50.0] {
            assert_eq!(
                prepared.within(&other, d),
                euclidean_within(&base, &other, d),
                "within {d}"
            );
        }
    }
}

#[test]
fn prepared_batch_predicates_match_per_call() {
    // The batch (spatial-join broadcast) methods must return, element-for-element,
    // exactly what the per-call methods do.
    let base = sq(0.0, 0.0, 10.0, 10.0);
    let prepared = PreparedGeometry::new(&base);
    let mut rng = Lcg(0xD1B5_4A32_D192_ED03);

    let others: Vec<Geometry> = (0..400)
        .map(|_| {
            let (x, y) = (rng.coord(40), rng.coord(40));
            let (w, h) = (1.0 + rng.coord(8).abs(), 1.0 + rng.coord(8).abs());
            sq(x, y, x + w, y + h)
        })
        .collect();

    let d = 3.0;
    let within = prepared.within_each(&others, d);
    let inter = prepared.intersects_each(&others);
    let dist = prepared.distance_each(&others);
    assert_eq!(within.len(), others.len());
    for (i, o) in others.iter().enumerate() {
        assert_eq!(within[i], prepared.within(o, d), "within_each[{i}]");
        assert_eq!(inter[i], prepared.intersects(o), "intersects_each[{i}]");
        let pd = prepared.distance(o);
        assert!((dist[i] - pd).abs() <= 1e-12 * pd.max(1.0), "distance_each[{i}]");
    }

    let coords: Vec<Coord> = (0..1000).map(|_| Coord::new(rng.coord(24), rng.coord(24))).collect();
    let batch = PointBatch::from_coords(&coords);
    let mask = prepared.contains_points(&batch);
    assert_eq!(mask.len(), coords.len());
    for (i, &c) in coords.iter().enumerate() {
        assert_eq!(mask[i], prepared.contains_point(c), "contains_points[{i}]");
    }
}

#[test]
fn prepared_contains_point_matches_polygon_with_hole() {
    let shell = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)];
    let hole = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)];
    let donut = Geometry::Polygon(Polygon::new(
        LineString::from(shell.to_vec()),
        vec![LineString::from(hole.to_vec())],
    ));
    let donut_poly = match &donut {
        Geometry::Polygon(p) => p.clone(),
        _ => unreachable!(),
    };
    let prepared = PreparedGeometry::new(&donut);

    let mut rng = Lcg(0x1234_5678_9ABC_DEF0);
    for _ in 0..3000 {
        // spans the shell interior, the hole, and the exterior
        let c = Coord::new(rng.coord(24), rng.coord(24));
        assert_eq!(
            prepared.contains_point(c),
            donut_poly.contains(&c),
            "contains_point at {c:?}"
        );
    }
}
