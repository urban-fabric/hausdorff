use super::*;
use hausdorff_base::earth::GeoPoint;
use hausdorff_planar::algebra::Orientation;
use std::f64::consts::PI;

fn sp(lat: f64, lon: f64) -> S2Point {
    S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
}

fn square(half: f64) -> Vec<S2Point> {
    vec![
        sp(-half, -half),
        sp(-half, half),
        sp(half, half),
        sp(half, -half),
    ]
}

#[test]
fn validation_and_normalization() {
    // CCW square accepted as-is.
    let c = ConvexLoop::from_vertices(square(5.0)).expect("convex");
    assert_eq!(c.num_vertices(), 4);
    assert!(c.area() < 2.0 * PI);

    // CW input (reversed) is normalized to CCW — same convex region.
    let mut rev = square(5.0);
    rev.reverse();
    let c2 = ConvexLoop::from_vertices(rev).expect("convex after normalize");
    assert!((c.area() - c2.area()).abs() < 1e-9);
    assert!(c2.contains_point(sp(0.0, 0.0)));

    // Concave L-shape rejected.
    assert!(
        ConvexLoop::from_vertices(vec![
            sp(0.0, 0.0),
            sp(6.0, 0.0),
            sp(6.0, 6.0),
            sp(3.0, 6.0),
            sp(3.0, 3.0),
            sp(0.0, 3.0),
        ])
        .is_none()
    );

    // Degenerate (< 3 vertices) rejected.
    assert!(ConvexLoop::from_vertices(vec![sp(0.0, 0.0), sp(0.0, 5.0)]).is_none());
}

#[test]
fn binary_search_agrees_with_general_loop() {
    // A many-vertex convex polygon (regular-ish ring near the equator) so the
    // O(log n) search exercises more than the trivial base case.
    let n = 16;
    let mut verts = Vec::new();
    for k in 0..n {
        let theta = 2.0 * PI * f64::from(k) / f64::from(n);
        // small circle of angular radius ~10° about (0,0)
        verts.push(sp(10.0 * theta.sin(), 10.0 * theta.cos()));
    }
    let cvx = ConvexLoop::from_vertices(verts.clone()).expect("convex ring");
    let general = cvx.as_loop().clone();

    let tol = 1e-3;
    let (mut inside, mut outside) = (0, 0);
    for i in -200..=200 {
        for j in -200..=200 {
            let p = sp(f64::from(i) * 0.1, f64::from(j) * 0.1);
            // Margin from the proven O(n) hemisphere predicate.
            let v = cvx.vertices();
            let m = v.len();
            let margin = (0..m)
                .map(|e| S2Point::orient(v[e], v[(e + 1) % m], p))
                .fold(f64::INFINITY, f64::min);
            if margin > tol {
                assert!(cvx.contains_point(p), "log-n missed an interior point");
                assert!(general.contains_point(p), "general loop disagrees inside");
                inside += 1;
            } else if margin < -tol {
                assert!(!cvx.contains_point(p), "log-n kept an exterior point");
                assert!(!general.contains_point(p), "general loop disagrees outside");
                outside += 1;
            }
        }
    }
    assert!(inside > 100 && outside > 100, "in={inside} out={outside}");
}

#[test]
fn convex_containment_of_loops() {
    let outer = ConvexLoop::from_vertices(square(10.0)).unwrap();
    let inner = ConvexLoop::from_vertices(square(3.0)).unwrap();
    assert!(outer.contains(&inner), "small box inside big box");
    assert!(!inner.contains(&outer), "big box not inside small box");

    // A box straddling the boundary is not contained.
    let straddle = ConvexLoop::from_vertices(vec![
        sp(8.0, 8.0),
        sp(8.0, 14.0),
        sp(12.0, 14.0),
        sp(12.0, 8.0),
    ])
    .unwrap();
    assert!(!outer.contains(&straddle));
}

#[test]
fn sat_disjointness() {
    let a = ConvexLoop::from_vertices(square(3.0)).unwrap();

    // Overlapping: a box straddling a's NE corner ([0,6]² ∩ [-3,3]² = [0,3]²).
    let over =
        ConvexLoop::from_vertices(vec![sp(0.0, 0.0), sp(0.0, 6.0), sp(6.0, 6.0), sp(6.0, 0.0)])
            .unwrap();
    assert!(a.intersects(&over));
    assert!(!a.disjoint(&over));

    // Disjoint: a box well to the east (lon 40–46).
    let far = ConvexLoop::from_vertices(vec![
        sp(-3.0, 40.0),
        sp(-3.0, 46.0),
        sp(3.0, 46.0),
        sp(3.0, 40.0),
    ])
    .unwrap();
    assert!(a.disjoint(&far));
    assert!(!a.intersects(&far));
    assert_eq!(a.disjoint(&far), far.disjoint(&a)); // symmetric

    // Containment ⟹ intersecting (not disjoint).
    let inner = ConvexLoop::from_vertices(square(1.0)).unwrap();
    assert!(a.contains(&inner));
    assert!(a.intersects(&inner));
}

/// True iff `v` is in `cp` or within tolerance of its boundary (clip outputs
/// land on parent boundaries, where strict `contains_point` is false).
fn in_or_on(cp: &ConvexLoop, v: S2Point) -> bool {
    if cp.contains_point(v) {
        return true;
    }
    let r = cp.vertices();
    let m = r.len();
    (0..m).any(|i| S2Point::orient(r[i], r[(i + 1) % m], v).abs() < 1e-9)
}

#[test]
fn intersection_overlap_containment_disjoint() {
    let a = ConvexLoop::from_vertices(square(3.0)).unwrap();

    // Overlap: a box straddling a's NE corner ([0,6]² ∩ [-3,3]² = [0,3]²).
    let b = ConvexLoop::from_vertices(vec![sp(0.0, 0.0), sp(0.0, 6.0), sp(6.0, 6.0), sp(6.0, 0.0)])
        .unwrap();
    let inter = a.intersection(&b).expect("non-empty overlap");
    // Every intersection vertex lies in (or on the boundary of) both parents.
    for &v in inter.vertices() {
        assert!(in_or_on(&a, v), "vertex escapes a");
        assert!(in_or_on(&b, v), "vertex escapes b");
    }
    // The shared interior point (1.5,1.5) is inside the intersection.
    assert!(inter.contains_point(sp(1.5, 1.5)));
    // Area is positive and no bigger than either parent.
    assert!(inter.area() > 0.0);
    assert!(inter.area() <= a.area() + 1e-9 && inter.area() <= b.area() + 1e-9);

    // Containment: inner ⊂ a ⟹ intersection ≈ inner (same area).
    let inner = ConvexLoop::from_vertices(square(1.0)).unwrap();
    let clipped = a.intersection(&inner).expect("non-empty");
    assert!((clipped.area() - inner.area()).abs() < 1e-6);

    // Disjoint ⟹ None.
    let far = ConvexLoop::from_vertices(vec![
        sp(-3.0, 40.0),
        sp(-3.0, 46.0),
        sp(3.0, 46.0),
        sp(3.0, 40.0),
    ])
    .unwrap();
    assert!(a.intersection(&far).is_none());
}
