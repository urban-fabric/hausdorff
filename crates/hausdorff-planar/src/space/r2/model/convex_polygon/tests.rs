use super::*;
use crate::space::r2::util::point_in_ring;
use hausdorff_traits::Orientation;

fn c(x: f64, y: f64) -> Coord {
    Coord::new(x, y)
}

fn square(s: f64) -> Vec<Coord> {
    vec![c(-s, -s), c(s, -s), c(s, s), c(-s, s)]
}

#[test]
fn validation_and_normalization() {
    let cp = ConvexPolygon::from_vertices(square(2.0)).expect("convex");
    assert_eq!(cp.num_vertices(), 4);
    assert!((cp.area() - 16.0).abs() < 1e-9);

    // CW input normalized to CCW (same region, same area).
    let mut rev = square(2.0);
    rev.reverse();
    let cp2 = ConvexPolygon::from_vertices(rev).expect("convex after normalize");
    assert!((cp2.area() - 16.0).abs() < 1e-9);
    assert!(cp2.contains(c(0.0, 0.0)));

    // A closed ring (trailing dup) is accepted (dup dropped).
    let mut closed = square(2.0);
    closed.push(closed[0]);
    assert!(ConvexPolygon::from_vertices(closed).is_some());

    // Concave (an arrow/dart) rejected.
    assert!(
        ConvexPolygon::from_vertices(vec![c(0.0, 0.0), c(4.0, 2.0), c(0.0, 4.0), c(1.0, 2.0)])
            .is_none()
    );

    // Degenerate / collinear rejected.
    assert!(ConvexPolygon::from_vertices(vec![c(0.0, 0.0), c(1.0, 1.0)]).is_none());
    assert!(ConvexPolygon::from_vertices(vec![c(0.0, 0.0), c(1.0, 1.0), c(2.0, 2.0)]).is_none());
}

#[test]
fn from_polygon_rejects_holes() {
    let solid = Polygon::new(LineString::new(square(3.0)), Vec::new());
    assert!(ConvexPolygon::from_polygon(&solid).is_some());

    let holed = Polygon::new(
        LineString::new(square(3.0)),
        vec![LineString::new(square(1.0))],
    );
    assert!(ConvexPolygon::from_polygon(&holed).is_none());
}

#[test]
fn binary_search_agrees_with_even_odd() {
    // A regular 12-gon so the O(log n) search exercises the fan.
    let n = 12;
    let mut verts = Vec::new();
    for k in 0..n {
        let t = 2.0 * std::f64::consts::PI * f64::from(k) / f64::from(n);
        verts.push(c(5.0 * t.cos(), 5.0 * t.sin()));
    }
    let cp = ConvexPolygon::from_vertices(verts).expect("convex 12-gon");
    let ring = cp.vertices().to_vec();

    let tol = 1e-6;
    let (mut inside, mut outside) = (0, 0);
    for i in -80..=80 {
        for j in -80..=80 {
            let p = c(f64::from(i) * 0.1, f64::from(j) * 0.1);
            let m = ring.len();
            let margin = (0..m)
                .map(|e| Coord::orient(ring[e], ring[(e + 1) % m], p))
                .fold(f64::INFINITY, f64::min);
            if margin > tol {
                assert!(cp.contains(p), "log-n missed an interior point");
                assert!(point_in_ring(p, &ring), "even-odd disagrees inside");
                inside += 1;
            } else if margin < -tol {
                assert!(!cp.contains(p), "log-n kept an exterior point");
                assert!(!point_in_ring(p, &ring), "even-odd disagrees outside");
                outside += 1;
            }
        }
    }
    assert!(inside > 100 && outside > 100, "in={inside} out={outside}");
}

#[test]
fn convex_containment_of_polygons() {
    let outer = ConvexPolygon::from_vertices(square(5.0)).unwrap();
    let inner = ConvexPolygon::from_vertices(square(2.0)).unwrap();
    assert!(outer.contains_polygon(&inner));
    assert!(!inner.contains_polygon(&outer));

    // A box poking outside is not contained.
    let poke =
        ConvexPolygon::from_vertices(vec![c(4.0, 4.0), c(7.0, 4.0), c(7.0, 7.0), c(4.0, 7.0)])
            .unwrap();
    assert!(!outer.contains_polygon(&poke));
}

#[test]
fn round_trip_to_polygon_is_convex() {
    let cp = ConvexPolygon::from_vertices(square(3.0)).unwrap();
    assert!(cp.to_polygon().is_convex());
}

#[test]
fn intersection_of_overlapping_squares() {
    // [-2,2]² ∩ [0,4]² = [0,2]² (a 2×2 box, area 4).
    let a = ConvexPolygon::from_vertices(square(2.0)).unwrap();
    let b = ConvexPolygon::from_vertices(vec![c(0.0, 0.0), c(4.0, 0.0), c(4.0, 4.0), c(0.0, 4.0)])
        .unwrap();
    let inter = a.intersection(&b).expect("overlap is non-empty");
    assert!((inter.area() - 4.0).abs() < 1e-9, "area = {}", inter.area());
    assert!(inter.contains(c(1.0, 1.0)));
    assert!(!inter.contains(c(3.0, 3.0)));
    // Intersection is commutative (same area).
    let inter2 = b.intersection(&a).unwrap();
    assert!((inter.area() - inter2.area()).abs() < 1e-9);
}

#[test]
fn intersection_containment_and_disjoint() {
    let outer = ConvexPolygon::from_vertices(square(5.0)).unwrap();
    let inner = ConvexPolygon::from_vertices(square(2.0)).unwrap();
    // inner ⊂ outer ⟹ intersection is inner (area 16).
    let clipped = outer.intersection(&inner).unwrap();
    assert!((clipped.area() - inner.area()).abs() < 1e-9);

    // Disjoint ⟹ None.
    let far = ConvexPolygon::from_vertices(vec![
        c(100.0, 100.0),
        c(102.0, 100.0),
        c(102.0, 102.0),
        c(100.0, 102.0),
    ])
    .unwrap();
    assert!(outer.intersection(&far).is_none());

    // Edge-touching only ⟹ degenerate ⟹ None.
    let touch =
        ConvexPolygon::from_vertices(vec![c(5.0, -1.0), c(7.0, -1.0), c(7.0, 1.0), c(5.0, 1.0)])
            .unwrap();
    assert!(outer.intersection(&touch).is_none());
}

#[test]
fn intersection_triangle_clips_square() {
    // A triangle overlapping a square clips to a convex polygon inside both.
    let sq = ConvexPolygon::from_vertices(square(3.0)).unwrap();
    let tri = ConvexPolygon::from_vertices(vec![c(-4.0, -4.0), c(4.0, -4.0), c(0.0, 4.0)]).unwrap();
    let inter = sq.intersection(&tri).expect("non-empty");
    // Every vertex of the intersection lies in both parents.
    for &v in inter.vertices() {
        assert!(sq.contains(v) || orient_on_boundary(&sq, v));
        assert!(tri.contains(v) || orient_on_boundary(&tri, v));
    }
    assert!(inter.area() > 0.0 && inter.area() <= tri.area());
}

/// True iff `v` is on (within tolerance of) some edge of `cp` — clip outputs
/// land on parent boundaries, where strict `contains` is `false` by design.
fn orient_on_boundary(cp: &ConvexPolygon, v: Coord) -> bool {
    let r = cp.vertices();
    let n = r.len();
    (0..n).any(|i| orient(r[i], r[(i + 1) % n], v).abs() < 1e-9)
}

#[test]
fn sat_disjointness_agrees_with_intersection() {
    let a = ConvexPolygon::from_vertices(square(2.0)).unwrap();
    // Overlapping.
    let b = ConvexPolygon::from_vertices(vec![c(1.0, 1.0), c(5.0, 1.0), c(5.0, 5.0), c(1.0, 5.0)])
        .unwrap();
    assert!(a.intersects(&b));
    assert!(!a.disjoint(&b));
    // Disjoint (separated along x).
    let far =
        ConvexPolygon::from_vertices(vec![c(10.0, 0.0), c(12.0, 0.0), c(12.0, 2.0), c(10.0, 2.0)])
            .unwrap();
    assert!(a.disjoint(&far));
    assert!(!a.intersects(&far));

    // Disjoint along a diagonal separating axis (corner-to-corner gap).
    let diag =
        ConvexPolygon::from_vertices(vec![c(3.0, 3.0), c(5.0, 3.0), c(5.0, 5.0), c(3.0, 5.0)])
            .unwrap();
    assert!(a.disjoint(&diag));

    // Consistency with `intersection`: disjoint ⟺ no (non-degenerate) overlap.
    for other in [&b, &far, &diag] {
        assert_eq!(a.intersects(other), a.intersection(other).is_some());
    }
    // Symmetry.
    assert_eq!(a.disjoint(&far), far.disjoint(&a));
}

#[test]
fn sat_triangle_vs_square() {
    let sq = ConvexPolygon::from_vertices(square(1.0)).unwrap();
    // Triangle poking into the square.
    let tri = ConvexPolygon::from_vertices(vec![c(0.0, 0.0), c(3.0, 0.5), c(0.5, 3.0)]).unwrap();
    assert!(sq.intersects(&tri));
    // Triangle entirely clear of the square.
    let away = ConvexPolygon::from_vertices(vec![c(5.0, 5.0), c(7.0, 5.0), c(6.0, 7.0)]).unwrap();
    assert!(sq.disjoint(&away));
}
