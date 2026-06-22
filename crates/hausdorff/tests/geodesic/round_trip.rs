use super::*;

#[test]
fn round_trip_inverse_then_direct() {
    // Self-contained breadth check across ALL pairs (incl. the 60-pair grid):
    // direct(inverse(a→b)) must return to b.
    for (a, b) in point_pairs() {
        let inv = Ellipsoid::WGS84.inverse(a, b);
        if inv.distance_m < 1e-9 {
            continue;
        }
        let fwd = Ellipsoid::WGS84.direct(a, inv.initial_bearing, inv.distance_m);
        assert!(
            (fwd.point.lat - b.lat).abs() < 1e-9,
            "lat {} vs {} for {a:?} → {b:?}",
            fwd.point.lat,
            b.lat
        );
        // At a pole the longitude (and arrival bearing) are degenerate.
        if b.lat.abs() >= 90.0 - 1e-9 {
            continue;
        }
        assert!(
            ang_close(fwd.point.lon, b.lon, 1e-9),
            "lon {} vs {} for {a:?} → {b:?}",
            fwd.point.lon,
            b.lon
        );
        assert!(
            ang_close(
                fwd.final_bearing.degrees(),
                inv.final_bearing.degrees(),
                1e-8
            ),
            "final bearing mismatch for {a:?} → {b:?}"
        );
    }
}

#[test]
fn r2_geodesic_metric_matches_ellipsoid() {
    // The r2 metric facade (x = lon, y = lat) must agree with the Ellipsoid API
    // it wraps — a self-contained facade-vs-core consistency check.
    use r2::{Bearing as _, Destination as _, Distance as _, Length as _};

    for (a, b) in point_pairs().into_iter().step_by(5) {
        let ours =
            r2::Geodesic.distance(r2::Point::new(a.lon, a.lat), r2::Point::new(b.lon, b.lat));
        let core = Ellipsoid::WGS84.distance_m(a, b);
        assert!(
            (ours - core).abs() < 1e-6,
            "{ours} vs {core} for {a:?} → {b:?}"
        );
    }

    let (a, b) = (g(40.6398, -73.7787), g(51.4775, -0.4614));
    let ours = r2::Geodesic.bearing(r2::Point::new(a.lon, a.lat), r2::Point::new(b.lon, b.lat));
    let core = Ellipsoid::WGS84.inverse(a, b).initial_bearing.degrees();
    assert!(ang_close(ours, core, 1e-8), "{ours} vs {core}");

    let ours = r2::Geodesic.destination(r2::Point::new(a.lon, a.lat), 51.0, 5_551_759.4);
    let core = Ellipsoid::WGS84
        .direct(a, Bearing::from_degrees(51.0), 5_551_759.4)
        .point;
    assert!((ours.x() - core.lon).abs() < 1e-9 && (ours.y() - core.lat).abs() < 1e-9);

    // Length over a polyline = sum of the segment inverses.
    let ls = r2::LineString::from(vec![(2.35, 48.85), (13.40, 52.52), (37.62, 55.75)]);
    let expect = r2::Geodesic.distance(r2::Point::new(2.35, 48.85), r2::Point::new(13.40, 52.52))
        + r2::Geodesic.distance(r2::Point::new(13.40, 52.52), r2::Point::new(37.62, 55.75));
    assert!((r2::Geodesic.length(&ls) - expect).abs() < 1e-9);
}

#[test]
fn geodesic_line_waypoints_match_direct() {
    // Sampling a line must agree with independent direct solves at every
    // distance, at a fraction of the cost.
    let origin = g(40.6398, -73.7787);
    let bearing = Bearing::from_degrees(51.0);
    let line = Ellipsoid::WGS84.line(origin, bearing);
    for i in 0..50 {
        let s = -2_000_000.0 + (i as f64) * 300_000.0; // backwards too
        let from_line = line.position_m(s);
        let from_direct = Ellipsoid::WGS84.direct(origin, bearing, s);
        assert!(
            (from_line.point.lat - from_direct.point.lat).abs() < 1e-13
                && (from_line.point.lon - from_direct.point.lon).abs() < 1e-13,
            "waypoint {s} diverges: {:?} vs {:?}",
            from_line.point,
            from_direct.point
        );
        assert_eq!(from_line.final_bearing, from_direct.final_bearing);
    }
}
