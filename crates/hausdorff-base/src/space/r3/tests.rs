use super::*;

#[test]
fn vector_is_generic_over_the_real_scalar_f32_parity() {
    // The ℝ³ `Vector` is no longer `f64`-locked: the same operations run at
    // `f32` and agree with the `f64` reference within single-precision
    // tolerance. This is the proof the geometry is generic over the
    // `Real` scalar trait.
    let a64 = Vector::<f64>::new(1.0, 2.0, 3.0);
    let b64 = Vector::<f64>::new(-4.0, 0.5, 6.0);
    let c64 = Vector::<f64>::new(0.0, -1.0, 2.0);
    let a32 = Vector::<f32>::new(1.0, 2.0, 3.0);
    let b32 = Vector::<f32>::new(-4.0, 0.5, 6.0);
    let c32 = Vector::<f32>::new(0.0, -1.0, 2.0);

    // dot product
    assert!((f64::from(a32.dot(b32)) - a64.dot(b64)).abs() < 1e-5);
    // norm
    assert!((f64::from(a32.norm()) - a64.norm()).abs() < 1e-5);
    // cross product, componentwise
    let x32 = a32.cross(b32);
    let x64 = a64.cross(b64);
    assert!((f64::from(x32.x) - x64.x).abs() < 1e-4);
    assert!((f64::from(x32.y) - x64.y).abs() < 1e-4);
    assert!((f64::from(x32.z) - x64.z).abs() < 1e-4);
    // normalized vector is unit length at f32 too
    let u32 = a32.normalize();
    assert!((f64::from(u32.norm()) - 1.0).abs() < 1e-6);
    // the scalar triple product / orientation predicate
    let o32 = a32.cross(b32).dot(c32);
    let o64 = a64.cross(b64).dot(c64);
    assert!((f64::from(o32) - o64).abs() < 1e-3);
}

#[test]
fn project_onto_tangent_is_perpendicular_with_pythagorean_length() {
    let v = Vector::<f64>::new(2.0, -3.0, 5.0);
    let n = Vector::<f64>::new(1.0, 1.0, 1.0).normalize();
    let proj = v.project_onto_tangent(n);
    // The projection is orthogonal to the normal (foot_perp_normal).
    assert!(proj.dot(n).abs() < 1e-12, "{}", proj.dot(n));
    // Its squared length is ‖v‖² − (v·n)² (tangent_normSq).
    let expected = v.dot(v) - v.dot(n).powi(2);
    assert!((proj.norm2() - expected).abs() < 1e-12, "{}", proj.norm2());
}

#[test]
fn project_onto_tangent_matches_inline_subtraction() {
    // The two duplicated inline copies in s2/arc.rs equal this primitive.
    let vp = Vector::<f64>::new(0.3, 0.7, -0.4);
    let n = Vector::<f64>::new(-0.2, 0.5, 0.84).normalize();
    let inline = vp.sub(n.mul(vp.dot(n)));
    assert_eq!(vp.project_onto_tangent(n), inline);
    // The latitude-extreme vertex ẑ − n·n_z is ẑ.project_onto_tangent(n).
    let z = Vector::<f64>::new(0.0, 0.0, 1.0);
    let vertex = Vector::<f64>::new(-n.x * n.z, -n.y * n.z, 1.0 - n.z * n.z);
    let proj = z.project_onto_tangent(n);
    assert!((proj.x - vertex.x).abs() < 1e-15);
    assert!((proj.y - vertex.y).abs() < 1e-15);
    assert!((proj.z - vertex.z).abs() < 1e-15);
}

#[test]
fn cos_angle_known_configurations_and_range() {
    let x = Vector::<f64>::new(1.0, 0.0, 0.0);
    let y = Vector::<f64>::new(0.0, 1.0, 0.0);
    // Orthogonal → 0.
    assert!(x.cos_angle(y).abs() < 1e-15);
    // Same direction (any positive scale) → 1.
    assert!((x.cos_angle(x.mul(7.0)) - 1.0).abs() < 1e-15);
    // Opposite → −1.
    assert!((x.cos_angle(x.mul(-3.0)) + 1.0).abs() < 1e-15);
    // 45° between (1,0,0) and (1,1,0) → 1/√2.
    let d = Vector::<f64>::new(1.0, 1.0, 0.0);
    assert!((x.cos_angle(d) - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-15);
    // Zero vector is graceful.
    assert_eq!(x.cos_angle(Vector::<f64>::new(0.0, 0.0, 0.0)), 0.0);
    // Always in [-1, 1] (Cauchy–Schwarz), for assorted vectors.
    let vecs = [
        Vector::<f64>::new(3.0, -1.0, 2.0),
        Vector::<f64>::new(-4.0, 5.0, 0.5),
        Vector::<f64>::new(0.1, 0.2, -9.0),
        Vector::<f64>::new(-2.0, -2.0, -2.0),
    ];
    for &u in &vecs {
        for &v in &vecs {
            let c = u.cos_angle(v);
            assert!((-1.0..=1.0).contains(&c), "{c}");
            // It is the cosine of the angle the robust `angle` method reports.
            assert!((c - u.angle(v).radians().cos()).abs() < 1e-12);
        }
    }
}

#[test]
fn gram_det_equals_triple_product_squared() {
    let triples = [
        (
            Vector::<f64>::new(1.0, 0.0, 0.0),
            Vector::<f64>::new(0.0, 1.0, 0.0),
            Vector::<f64>::new(0.0, 0.0, 1.0),
        ),
        (
            Vector::<f64>::new(2.0, -1.0, 3.0),
            Vector::<f64>::new(0.5, 4.0, -2.0),
            Vector::<f64>::new(-3.0, 1.0, 1.0),
        ),
        (
            Vector::<f64>::new(1.0, 1.0, 0.0),
            Vector::<f64>::new(0.0, 1.0, 1.0),
            Vector::<f64>::new(1.0, 0.0, 1.0),
        ),
        // A coplanar (degenerate) triple → triple product 0 → gram_det 0.
        (
            Vector::<f64>::new(1.0, 2.0, 3.0),
            Vector::<f64>::new(2.0, 4.0, 6.0),
            Vector::<f64>::new(-1.0, 0.0, 1.0),
        ),
    ];
    for (a, b, c) in triples {
        let triple = a.dot(b.cross(c));
        let gram = Vector::gram_det(a, b, c);
        assert!(
            (triple * triple - gram).abs() < 1e-9,
            "triple²={} vs gram={gram}",
            triple * triple
        );
    }
}

#[test]
fn gram_det_unit_vectors_match_spherical_cosine_form() {
    // For unit vectors: gram_det = 1 − Σcos² + 2∏cos.
    let unit = [
        Vector::<f64>::new(0.3, -0.5, 0.8).normalize(),
        Vector::<f64>::new(-0.7, 0.2, 0.4).normalize(),
        Vector::<f64>::new(0.1, 0.9, -0.2).normalize(),
    ];
    let (a, b, c) = (unit[0], unit[1], unit[2]);
    let (ab, bc, ca) = (a.dot(b), b.dot(c), c.dot(a));
    let spherical = 1.0 - ab * ab - bc * bc - ca * ca + 2.0 * ab * bc * ca;
    assert!((Vector::gram_det(a, b, c) - spherical).abs() < 1e-12);
}

#[test]
fn aabb_contains_overlaps_and_nesting() {
    // Corners accepted in any order and normalized.
    let b = Aabb::<f64>::new(
        Vector::<f64>::new(4.0, 0.0, 6.0),
        Vector::<f64>::new(0.0, 4.0, 2.0),
    );
    assert_eq!(b.min(), Vector::<f64>::new(0.0, 0.0, 2.0));
    assert_eq!(b.max(), Vector::<f64>::new(4.0, 4.0, 6.0));
    assert_eq!(b.center(), Vector::<f64>::new(2.0, 2.0, 4.0));
    // Membership (AABB3.contains): inside, on a face, and outside (per axis).
    assert!(b.contains(Vector::<f64>::new(2.0, 2.0, 4.0)));
    assert!(b.contains(Vector::<f64>::new(0.0, 4.0, 6.0))); // corner, inclusive
    assert!(!b.contains(Vector::<f64>::new(2.0, 2.0, 7.0))); // past max z
    // Overlap (overlaps_of_mem3) vs. a disjoint box (prune_sound3).
    let touching = Aabb::<f64>::new(
        Vector::<f64>::new(4.0, 4.0, 6.0),
        Vector::<f64>::new(6.0, 6.0, 8.0),
    );
    let far = Aabb::<f64>::new(
        Vector::<f64>::new(10.0, 10.0, 10.0),
        Vector::<f64>::new(12.0, 12.0, 12.0),
    );
    assert!(b.overlaps(&touching)); // shares the corner
    assert!(!b.overlaps(&far));
    // Nesting (contains_of_subset3) and its transitive union.
    let inner = Aabb::<f64>::new(
        Vector::<f64>::new(1.0, 1.0, 3.0),
        Vector::<f64>::new(3.0, 3.0, 5.0),
    );
    assert!(b.contains_box(&inner));
    assert!(!inner.contains_box(&b));
    assert!(b.union(&far).contains_box(&b));
    assert!(b.union(&far).contains_box(&far));
    // Expansion grows every face.
    assert!(
        b.expanded(1.0)
            .contains(Vector::<f64>::new(-0.5, -0.5, 1.5))
    );
}

#[test]
fn plane_side_test_equals_orient4() {
    // The z = 0 plane (normal +z) via three points.
    let (a, b, c) = (
        Vector::<f64>::new(0.0, 0.0, 0.0),
        Vector::<f64>::new(1.0, 0.0, 0.0),
        Vector::<f64>::new(0.0, 1.0, 0.0),
    );
    let plane = Plane::from_points(a, b, c);
    assert_eq!(plane.normal(), Vector::<f64>::new(0.0, 0.0, 1.0));
    // Sides: above (+z), below (−z), and on the plane.
    assert_eq!(plane.side_of(Vector::<f64>::new(0.5, 0.5, 2.0)), 1);
    assert_eq!(plane.side_of(Vector::<f64>::new(0.5, 0.5, -2.0)), -1);
    assert_eq!(plane.side_of(Vector::<f64>::new(0.5, 0.5, 0.0)), 0);
    assert!(plane.contains(Vector::<f64>::new(9.0, -3.0, 0.0), 1e-12));
    // Perpendicular signed distance is the bare z (unit normal).
    assert!((plane.signed_distance(Vector::<f64>::new(0.0, 0.0, 3.0)) - 3.0).abs() < 1e-12);

    // `eval(d)` equals orient4(a, b, c, d) = ((b−a)×(c−a))·(d−a) exactly.
    let orient4 =
        |a: Vector, b: Vector, c: Vector, d: Vector| b.sub(a).cross(c.sub(a)).dot(d.sub(a));
    for d in [
        Vector::<f64>::new(0.3, 0.4, 1.7),
        Vector::<f64>::new(-2.0, 1.0, -0.5),
        Vector::<f64>::new(5.0, 5.0, 0.0),
    ] {
        assert!((plane.eval(d) - orient4(a, b, c, d)).abs() < 1e-12);
    }
}

#[test]
fn plane_from_point_normal_and_straddle() {
    let plane = Plane::from_point_normal(
        Vector::<f64>::new(0.0, 0.0, 1.0),
        Vector::<f64>::new(0.0, 0.0, 2.0),
    );
    // Offset = normal·point = 2; the plane is z = 1.
    assert_eq!(plane.offset(), 2.0);
    assert_eq!(plane.side_of(Vector::<f64>::new(0.0, 0.0, 5.0)), 1);
    assert_eq!(plane.side_of(Vector::<f64>::new(0.0, 0.0, 0.0)), -1);
    // Opposite sides ⟹ the segment crosses the plane (orient4_lerp_last):
    // eval is affine along the segment, so it has a zero in (0, 1).
    let p = Vector::<f64>::new(0.0, 0.0, 0.0);
    let q = Vector::<f64>::new(0.0, 0.0, 4.0);
    assert!(plane.side_of(p) * plane.side_of(q) < 0);
    let t = -plane.eval(p) / (plane.eval(q) - plane.eval(p));
    assert!((0.0..=1.0).contains(&t));
    let crossing = p.add(q.sub(p).mul(t));
    assert!(plane.contains(crossing, 1e-12));
}
