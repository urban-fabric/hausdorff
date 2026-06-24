use super::*;

fn g(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

/// A CCW rectangular ring (SW, SE, NE, NW) of the lat/lon box.
fn box_ring(s: f64, w: f64, n: f64, e: f64) -> Vec<GeoPoint> {
    vec![g(s, w), g(s, e), g(n, e), g(n, w)]
}

/// Karney (geodesic) area/perimeter of a ring — the ground truth to validate the
/// native measures against.
fn karney(ring: &[GeoPoint]) -> (f64, f64) {
    Ellipsoid::WGS84.ring_perimeter_area_unsigned_m(ring)
}

#[test]
fn area_matches_the_geodesic_area() {
    let ring = box_ring(30.0, 0.0, 40.0, 10.0);
    let poly = GeoPolygon::from_exterior(ring.clone());
    let (_, karney_area) = karney(&ring);
    // The authalic-sphere area reproduces the geodesic area to < 0.01 %.
    assert!(
        (poly.area_m2() - karney_area).abs() / karney_area < 1e-4,
        "area {} vs karney {karney_area}",
        poly.area_m2()
    );
}

#[test]
fn perimeter_matches_the_geodesic_perimeter() {
    let ring = box_ring(30.0, 0.0, 40.0, 10.0);
    let poly = GeoPolygon::from_exterior(ring.clone());
    let (karney_perim, _) = karney(&ring);
    // Same Karney geodesic edges — agree to sub-millimetre.
    assert!(
        (poly.perimeter_m() - karney_perim).abs() < 1e-3,
        "perimeter {} vs karney {karney_perim}",
        poly.perimeter_m()
    );
    // Open vs already-closed ring give the same perimeter.
    let mut closed = ring.clone();
    closed.push(ring[0]);
    assert!((GeoPolygon::from_exterior(closed).perimeter_m() - poly.perimeter_m()).abs() < 1e-6);
}

#[test]
fn holes_subtract_from_area_and_containment() {
    let shell = box_ring(30.0, 0.0, 40.0, 10.0);
    let hole = box_ring(34.0, 4.0, 36.0, 6.0);
    let poly = GeoPolygon::new(shell.clone(), vec![hole.clone()]);
    let (_, shell_a) = karney(&shell);
    let (_, hole_a) = karney(&hole);
    assert!(
        (poly.area_m2() - (shell_a - hole_a)).abs() / (shell_a - hole_a) < 1e-4,
        "area {} vs {}",
        poly.area_m2(),
        shell_a - hole_a
    );
    // Containment respects the hole.
    assert!(poly.contains(g(31.0, 1.0)), "inside shell, outside hole");
    assert!(!poly.contains(g(35.0, 5.0)), "inside the hole");
    assert!(!poly.contains(g(50.0, 50.0)), "outside the shell");
}

#[test]
fn containment_is_orientation_robust() {
    // A clockwise-wound exterior still contains its interior (ccw normalization).
    let cw = vec![g(30.0, 0.0), g(40.0, 0.0), g(40.0, 10.0), g(30.0, 10.0)];
    let poly = GeoPolygon::from_exterior(cw);
    assert!(poly.contains(g(35.0, 5.0)));
    assert!(!poly.contains(g(35.0, 20.0)));
}

#[test]
fn multipolygon_sums_area_and_unions_containment() {
    let a = GeoPolygon::from_exterior(box_ring(0.0, 0.0, 5.0, 5.0));
    let b = GeoPolygon::from_exterior(box_ring(0.0, 50.0, 5.0, 55.0));
    let area_sum = a.area_m2() + b.area_m2();
    let mp = GeoMultiPolygon::new(vec![a, b]);
    assert!((mp.area_m2() - area_sum).abs() < 1e-3);
    assert!(
        (mp.perimeter_m()
            - (GeoPolygon::from_exterior(box_ring(0.0, 0.0, 5.0, 5.0)).perimeter_m()
                + GeoPolygon::from_exterior(box_ring(0.0, 50.0, 5.0, 55.0)).perimeter_m()))
        .abs()
            < 1e-3
    );
    // A point in either member is contained; one between them is not.
    assert!(mp.contains(g(2.5, 2.5)));
    assert!(mp.contains(g(2.5, 52.5)));
    assert!(!mp.contains(g(2.5, 25.0)));
}

#[test]
fn intersects_overlapping_but_not_disjoint() {
    let a = GeoPolygon::from_exterior(box_ring(0.0, 0.0, 10.0, 10.0));
    let overlapping = GeoPolygon::from_exterior(box_ring(5.0, 5.0, 15.0, 15.0));
    let nested = GeoPolygon::from_exterior(box_ring(2.0, 2.0, 4.0, 4.0));
    let disjoint = GeoPolygon::from_exterior(box_ring(0.0, 50.0, 10.0, 60.0));
    assert!(a.intersects(&overlapping), "edges cross");
    assert!(a.intersects(&nested), "one contains the other");
    assert!(!a.intersects(&disjoint), "far apart");
}

#[test]
fn contains_polygon_nested_crossing_and_holes() {
    let outer = GeoPolygon::from_exterior(box_ring(0.0, 0.0, 10.0, 10.0));
    let inner = GeoPolygon::from_exterior(box_ring(2.0, 2.0, 4.0, 4.0)); // strictly inside
    let crossing = GeoPolygon::from_exterior(box_ring(5.0, 5.0, 15.0, 15.0)); // straddles edge
    let outside = GeoPolygon::from_exterior(box_ring(20.0, 20.0, 30.0, 30.0)); // disjoint
    assert!(outer.contains_polygon(&inner));
    assert!(!outer.contains_polygon(&crossing));
    assert!(!outer.contains_polygon(&outside));
    assert!(!inner.contains_polygon(&outer)); // not the other way round
    // A polygon with a hole does not contain a shape sitting in the hole.
    let holed = GeoPolygon::new(
        box_ring(0.0, 0.0, 10.0, 10.0),
        vec![box_ring(3.0, 3.0, 7.0, 7.0)],
    );
    let in_hole = GeoPolygon::from_exterior(box_ring(4.0, 4.0, 6.0, 6.0));
    assert!(!holed.contains_polygon(&in_hole));
}

#[test]
fn conversion_helpers_round_trip() {
    let ring = box_ring(30.0, 0.0, 40.0, 10.0);
    // To unit sphere: every point is a unit vector and recovers its lat/lon.
    for (p, s) in ring.iter().zip(geo_ring_to_s2(&ring)) {
        assert!((s.dot(s) - 1.0).abs() < 1e-12, "unit vector");
        let back = s.to_geo_point();
        assert!((back.lat - p.lat).abs() < 1e-9 && (back.lon - p.lon).abs() < 1e-9);
    }
    // To ECEF on the surface: ecef → geodetic recovers the lat/lon.
    for (p, e) in ring.iter().zip(geo_ring_to_ecef(&ring, Ellipsoid::WGS84)) {
        let (back, h) = Ellipsoid::WGS84.ecef_to_geodetic(e);
        assert!((back.lat - p.lat).abs() < 1e-9 && (back.lon - p.lon).abs() < 1e-9);
        assert!(h.abs() < 1e-6, "on the surface");
    }
}
