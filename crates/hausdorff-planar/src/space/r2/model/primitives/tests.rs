use super::*;

#[test]
fn ellipse_membership_boundary_and_convexity() {
    let e = Ellipse::new(Coord::new(1.0, 2.0), 3.0, 2.0);
    // Centre inside; the two semi-axis tips lie on the (inclusive) boundary.
    assert!(e.contains_coord(Coord::new(1.0, 2.0)));
    assert!(e.contains_coord(Coord::new(4.0, 2.0)));
    assert!(e.contains_coord(Coord::new(1.0, 4.0)));
    // Just outside the boundary is excluded.
    assert!(!e.contains_coord(Coord::new(4.001, 2.0)));
    assert!(!e.contains_coord(Coord::new(1.0, 4.001)));
    // Convexity (inEllipse_convex): the midpoint of two contained points is
    // contained.
    let mid = Coord::new((4.0 + 1.0) * 0.5, (2.0 + 4.0) * 0.5);
    assert!(e.contains_coord(mid));
}

#[test]
fn circle_is_a_degenerate_ellipse() {
    let c = Ellipse::new(Coord::new(0.0, 0.0), 5.0, 5.0);
    assert!(c.contains_coord(Coord::new(3.0, 4.0))); // on the radius-5 circle
    assert!(!c.contains_coord(Coord::new(3.0, 4.1)));
}

#[test]
fn ellipse_area_and_bounds() {
    let e = Ellipse::new(Coord::new(1.0, 2.0), 3.0, 2.0);
    assert!((e.area() - std::f64::consts::PI * 6.0).abs() < 1e-12);
    let r = e.bounding_rect();
    assert_eq!(r.min(), Coord::new(-2.0, 0.0));
    assert_eq!(r.max(), Coord::new(4.0, 4.0));
}

#[test]
fn confidence_ellipse_axis_aligned() {
    // Diagonal covariance → axis-aligned; eigenvalues are the variances.
    let e = ConfidenceEllipse::from_covariance(Coord::new(0.0, 0.0), 4.0, 0.0, 1.0, 1.0);
    let (lmax, lmin) = e.eigenvalues();
    assert!((lmax - 4.0).abs() < 1e-12 && (lmin - 1.0).abs() < 1e-12);
    assert!((e.semi_major() - 2.0).abs() < 1e-12);
    assert!((e.semi_minor() - 1.0).abs() < 1e-12);
    assert!(e.rotation_rad().abs() < 1e-12);
    assert!((e.area() - std::f64::consts::PI * 2.0).abs() < 1e-12);
    // Boundary at the major-axis tip (k·σ_x = 2 along x); centre inside.
    assert!(e.contains_coord(Coord::new(2.0, 0.0)));
    assert!(!e.contains_coord(Coord::new(2.001, 0.0)));
    assert!(e.contains_coord(Coord::new(0.0, 0.0)));
}

#[test]
fn confidence_ellipse_correlated_is_rotated() {
    // Equicorrelated covariance → 45° major axis, eigenvalues 3 and 1.
    let e = ConfidenceEllipse::from_covariance(Coord::new(0.0, 0.0), 2.0, 1.0, 2.0, 1.0);
    let (lmax, lmin) = e.eigenvalues();
    assert!((lmax - 3.0).abs() < 1e-12 && (lmin - 1.0).abs() < 1e-12);
    assert!((e.rotation_rad() - std::f64::consts::FRAC_PI_4).abs() < 1e-12);
    // Area = π·k²·√(det Σ), det = 3.
    assert!((e.area() - std::f64::consts::PI * 3.0_f64.sqrt()).abs() < 1e-12);
}

#[test]
fn confidence_ellipse_eigenstructure() {
    let e = ConfidenceEllipse::from_covariance(Coord::new(0.0, 0.0), 2.0, 1.0, 2.0, 1.0);
    let (lmax, lmin) = e.eigenvalues();
    // Each eigenvalue solves the characteristic equation det(M − λI) = 0
    // (riemann sym2_eigenvalue).
    let (sxx, sxy, syy) = e.covariance();
    for l in [lmax, lmin] {
        let chi = (sxx - l) * (syy - l) - sxy * sxy;
        assert!(chi.abs() < 1e-12, "char eq {chi} at λ={l}");
    }
    // The eigenvalue product is det Σ (riemann sym2_eig_prod_eq_det).
    assert!((lmax * lmin - e.det()).abs() < 1e-12);
    // The two semi-axis vectors are orthogonal with the right lengths.
    let (maj, min) = e.axes();
    assert!((maj[0] * min[0] + maj[1] * min[1]).abs() < 1e-12);
    assert!((maj[0].hypot(maj[1]) - e.semi_major()).abs() < 1e-12);
    assert!((min[0].hypot(min[1]) - e.semi_minor()).abs() < 1e-12);
}

#[test]
fn ellipse_foci_and_eccentricity() {
    let e = Ellipse::new(Coord::new(0.0, 0.0), 3.0, 2.0); // major along x, c = √5
    let c = 5.0_f64.sqrt();
    let (f1, f2) = e.foci();
    assert!((f1.x + c).abs() < 1e-12 && f1.y.abs() < 1e-12);
    assert!((f2.x - c).abs() < 1e-12 && f2.y.abs() < 1e-12);
    assert!((e.eccentricity() - c / 3.0).abs() < 1e-12);
    // A circle has zero eccentricity and coincident foci at the centre.
    let circle = Ellipse::new(Coord::new(1.0, 1.0), 5.0, 5.0);
    assert!(circle.eccentricity().abs() < 1e-12);
    let (g1, g2) = circle.foci();
    assert_eq!(g1, Coord::new(1.0, 1.0));
    assert_eq!(g2, Coord::new(1.0, 1.0));
}

#[test]
fn ellipse_reflection_property() {
    // The normal makes equal angles with the two focal radii: a ray from one
    // focus reflects to the other (riemann Hausdorff/Ellipse.ellipse_reflection).
    let e = Ellipse::new(Coord::new(0.0, 0.0), 3.0, 2.0);
    let (f1, f2) = e.foci();
    for &t in &[0.3, 1.0, 2.5, 4.0, 5.5] {
        let p = e.point_at_angle(t);
        let n = e.normal_at(t);
        let u1 = [p.x - f1.x, p.y - f1.y];
        let u2 = [p.x - f2.x, p.y - f2.y];
        let cos1 = (n[0] * u1[0] + n[1] * u1[1]) / u1[0].hypot(u1[1]);
        let cos2 = (n[0] * u2[0] + n[1] * u2[1]) / u2[0].hypot(u2[1]);
        assert!((cos1 - cos2).abs() < 1e-12, "t={t}: {cos1} vs {cos2}");
        // Tangent ⟂ normal.
        let tan = e.tangent_at(t);
        assert!((n[0] * tan[0] + n[1] * tan[1]).abs() < 1e-12);
    }
}

#[test]
fn ellipse_tangent_line() {
    let e = Ellipse::new(Coord::new(1.0, -2.0), 3.0, 2.0);
    for &t in &[0.4, 1.3, 2.7, 4.5] {
        let [a, b, c] = e.tangent_line_at(t);
        let p = e.point_at_angle(t);
        // The boundary point lies on its tangent line.
        assert!((a * p.x + b * p.y - c).abs() < 1e-9);
        // The line normal (a, b) is parallel to the unit surface normal.
        let n = e.normal_at(t);
        let len = a.hypot(b);
        assert!((a / len - n[0]).abs() < 1e-12 && (b / len - n[1]).abs() < 1e-12);
        // The tangent direction runs along the line (⟂ the line normal).
        let tan = e.tangent_at(t);
        assert!((a * tan[0] + b * tan[1]).abs() < 1e-12);
    }
}

#[test]
fn ellipse_focus_directrix_property() {
    // Focus-distance / directrix-distance = eccentricity (major-x ellipse).
    let e = Ellipse::new(Coord::new(0.0, 0.0), 3.0, 2.0);
    let (_, f2) = e.foci(); // (c, 0)
    let d = e.directrix_offset(); // matching directrix at x = d
    let ecc = e.eccentricity();
    for &t in &[0.2, 1.1, 2.0, 3.5, 5.0] {
        let p = e.point_at_angle(t);
        let focus_dist = (p.x - f2.x).hypot(p.y - f2.y);
        let directrix_dist = (d - p.x).abs();
        assert!(
            (focus_dist - ecc * directrix_dist).abs() < 1e-12,
            "t={t}: {focus_dist} vs {}",
            ecc * directrix_dist
        );
    }
    assert!((e.semi_latus_rectum() - 4.0 / 3.0).abs() < 1e-12);
}

#[test]
fn ellipse_to_polygon_inscribes_and_approximates_area() {
    use crate::space::r2::Area;
    let e = Ellipse::new(Coord::new(1.0, -2.0), 3.0, 2.0);
    // The inscribed n-gon's area converges to π·a·b = 6π as n grows.
    let exact = std::f64::consts::PI * 6.0;
    let area = e.to_polygon(4000).unsigned_area();
    assert!(
        (area - exact).abs() / exact < 1e-5,
        "area {area} vs {exact}"
    );
    // `n` is clamped to ≥ 3 — no panic, still a valid positive-area polygon.
    assert!(e.to_polygon(0).unsigned_area() > 0.0);
}

#[test]
fn triangle_circumcircle_is_equidistant() {
    // Right triangle (0,0),(4,0),(0,3): the circumcentre is the hypotenuse
    // midpoint (2, 1.5) and the radius half the hypotenuse, 2.5.
    let t = Triangle(
        Coord::new(0.0, 0.0),
        Coord::new(4.0, 0.0),
        Coord::new(0.0, 3.0),
    );
    let (c, r) = t.circumcircle().expect("non-degenerate");
    assert!((c.x - 2.0).abs() < 1e-12 && (c.y - 1.5).abs() < 1e-12);
    assert!((r - 2.5).abs() < 1e-12);
    // Equidistant from every vertex (circumcenter_exists).
    for v in t.to_array() {
        assert!(((c.x - v.x).hypot(c.y - v.y) - r).abs() < 1e-12);
    }
    // Collinear vertices have no circumcircle.
    let flat = Triangle(
        Coord::new(0.0, 0.0),
        Coord::new(1.0, 1.0),
        Coord::new(2.0, 2.0),
    );
    assert!(flat.circumcenter().is_none());
    assert!(flat.circumcircle().is_none());
}

#[test]
fn line_perpendicular_bisector_is_the_equidistant_locus() {
    let l = Line::new(Coord::new(0.0, 0.0), Coord::new(4.0, 2.0));
    let [a, b, c] = l.perpendicular_bisector();
    // Its normal (A, B) is the chord end − start; the midpoint lies on it.
    assert!((a - 4.0).abs() < 1e-12 && (b - 2.0).abs() < 1e-12);
    let m = Coord::new(2.0, 1.0);
    assert!((a * m.x + b * m.y - c).abs() < 1e-12);
    // Every point on the line (stepping along its direction (−B, A)) is
    // equidistant from the two endpoints — `equidistant_iff`.
    for t in [-2.0, 0.5, 3.0] {
        let p = Coord::new(m.x - b * t, m.y + a * t);
        let da = (p.x - 0.0).hypot(p.y - 0.0);
        let db = (p.x - 4.0).hypot(p.y - 2.0);
        assert!((da - db).abs() < 1e-9, "t={t}: {da} vs {db}");
    }
}
