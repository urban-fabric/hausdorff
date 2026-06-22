use super::*;

fn deg(d: f64) -> Angle {
    Angle::from_degrees(d)
}

#[test]
fn membership_wraps_around_zero() {
    // A 30°-radius sector centred at 350°: contains 340..=20 (mod 360).
    let s = Sector::from_center_radius(deg(350.0), deg(30.0));
    assert!(s.contains_angle(deg(350.0))); // centre
    assert!(s.contains_angle(deg(10.0))); // 20° away, across 0°
    assert!(s.contains_angle(deg(340.0))); // 10° away
    assert!(s.contains_angle(deg(20.0))); // on the boundary
    assert!(!s.contains_angle(deg(40.0))); // 50° away — outside
    assert!(!s.contains_angle(deg(180.0))); // antipodal-ish — outside
}

#[test]
fn empty_and_full() {
    let e = Sector::empty(deg(10.0));
    assert!(e.is_empty());
    assert!(!e.contains_angle(deg(10.0))); // even the "centre" is excluded
    let f = Sector::full(deg(10.0));
    assert!(f.is_full());
    assert!(f.contains_angle(deg(200.0)));
    // A point sector contains exactly its centre (to within rounding).
    let p = Sector::from_angle(deg(45.0));
    assert!(p.contains_angle(deg(45.0)));
    assert!(!p.contains_angle(deg(45.001)));
}

#[test]
fn angular_distance_is_the_shorter_arc() {
    let s = Sector::from_center_radius(deg(10.0), deg(5.0));
    // 350° is 20° from 10° (the short way), not 340°.
    assert!((s.angular_distance(deg(350.0)).degrees() - 20.0).abs() < 1e-9);
    assert!(s.angular_distance(deg(10.0)).radians().abs() < 1e-12);
    // Always in [0, π].
    for d in [0.0, 90.0, 179.0, 181.0, 270.0, 359.0] {
        let dist = s.angular_distance(deg(d)).radians();
        assert!((0.0..=PI + 1e-9).contains(&dist));
    }
}

#[test]
fn intersects_and_contains_with_pruning() {
    let a = Sector::from_center_radius(deg(0.0), deg(20.0));
    let b = Sector::from_center_radius(deg(30.0), deg(20.0)); // 30 ≤ 40 → overlap
    let c = Sector::from_center_radius(deg(90.0), deg(20.0)); // 90 > 40 → disjoint
    assert!(a.intersects(b));
    assert!(b.intersects(a));
    assert!(!a.intersects(c)); // arcs_disjoint_of_far prunes
    // Containment: a wide sector swallows a narrow nearby one.
    let wide = Sector::from_center_radius(deg(0.0), deg(60.0));
    let narrow = Sector::from_center_radius(deg(10.0), deg(20.0)); // 10 + 20 ≤ 60
    assert!(wide.contains_sector(narrow));
    assert!(!narrow.contains_sector(wide));
    // Containment ⟹ intersection (a genuine consistency check).
    assert!(wide.intersects(narrow));
}

#[test]
fn contained_sector_points_are_all_inside() {
    // The triangle inequality (circDist_triangle) guarantee: if `wide`
    // contains `narrow`, every angle of `narrow` is in `wide`.
    let wide = Sector::from_center_radius(deg(100.0), deg(50.0));
    let narrow = Sector::from_center_radius(deg(120.0), deg(20.0));
    assert!(wide.contains_sector(narrow));
    for off in [-20.0, -10.0, 0.0, 10.0, 20.0] {
        assert!(wide.contains_angle(deg(120.0 + off)));
    }
}

#[test]
fn expanded_and_complement() {
    let s = Sector::from_center_radius(deg(0.0), deg(20.0));
    let big = s.expanded(deg(30.0));
    assert!((big.radius().degrees() - 50.0).abs() < 1e-9);
    assert!(big.contains_angle(deg(45.0))); // grew to include it
    assert!(!s.contains_angle(deg(45.0)));
    // Expansion only adds angles (InArc.mono): everything in `s` stays in.
    for d in [-15.0, 0.0, 15.0] {
        assert!(big.contains_angle(deg(d)));
    }
    // The complement covers exactly what the sector excludes.
    let comp = s.complement();
    assert!((comp.radius().degrees() - 160.0).abs() < 1e-9);
    assert!(comp.contains_angle(deg(180.0)));
    assert!(!comp.contains_angle(deg(0.0)));
    // Expansion past π saturates to a full circle.
    assert!(s.expanded(Angle::PI).is_full());
}

#[test]
fn width_tracks_radius() {
    let s = Sector::from_center_radius(deg(0.0), deg(20.0));
    assert!((s.width().degrees() - 40.0).abs() < 1e-9);
    assert_eq!(Sector::empty(deg(0.0)).width(), Angle::ZERO);
    // Full sector spans the whole circle.
    assert!((Sector::full(deg(0.0)).width().radians() - Angle::TAU.radians()).abs() < 1e-9);
}
