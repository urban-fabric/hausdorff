//! Ellipsoidal (Karney) geodesics: differential conformance of
//! `Ellipsoid::{inverse, direct}` against `geographiclib-rs` (the porting
//! source and `geo`'s geodesic backend), plus consistency between the
//! spherical and ellipsoidal layers.

use geographiclib_rs::{DirectGeodesic, Geodesic as OracleGeodesic, InverseGeodesic};
use hausdorff::earth::{Bearing, Ellipsoid, GeoPoint, haversine_m};
use hausdorff::space::r2;

fn g(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

/// A deterministic sweep of point pairs covering the algorithm's regimes:
/// short lines, long lines, meridional, equatorial, polar, sign quadrants,
/// antimeridian-straddling, and the nearly-antipodal Newton/astroid path.
fn point_pairs() -> Vec<(GeoPoint, GeoPoint)> {
    let mut pairs = vec![
        // Short and city-scale lines.
        (g(48.8566, 2.3522), g(48.8570, 2.3530)),
        (g(40.7128, -74.0060), g(40.7614, -73.9776)),
        // Long lines across quadrants.
        (g(40.6398, -73.7787), g(51.4775, -0.4614)), // JFK → LHR
        (g(-33.8688, 151.2093), g(35.6762, 139.6503)), // SYD → TYO
        (g(37.7749, -122.4194), g(-33.9249, 18.4241)), // SFO → CPT
        // Meridional, including pole-crossing.
        (g(10.0, 20.0), g(50.0, 20.0)),
        (g(80.0, 30.0), g(80.0, 210.0 - 360.0)),
        (g(-90.0, 0.0), g(90.0, 0.0)),
        (g(-90.0, 45.0), g(78.0, 45.0)),
        // Equatorial.
        (g(0.0, -10.0), g(0.0, 130.0)),
        (g(0.0, 170.0), g(0.0, -170.0)),
        // Nearly antipodal (Vincenty-failure regime).
        (g(0.5, 0.0), g(-0.4999, 179.9999)),
        (g(30.0, 0.0), g(-30.0001, 179.99)),
        (g(5.0, 10.0), g(-5.0, -169.99)),
        // Exactly antipodal-ish on the equator.
        (g(0.0, 0.0), g(0.0, 179.5)),
        // Antimeridian straddles.
        (g(60.0, 179.5), g(59.0, -179.5)),
        (g(-45.0, 178.0), g(-44.0, -178.0)),
        // Coincident points.
        (g(12.34, 56.78), g(12.34, 56.78)),
        // Poles to elsewhere.
        (g(90.0, 0.0), g(45.0, 90.0)),
        (g(-90.0, 12.0), g(-45.0, -120.0)),
    ];
    // A deterministic pseudo-grid for breadth (no RNG: fixed multipliers).
    for i in 0..60 {
        let lat1 = -80.0 + (i as f64) * 2.7;
        let lon1 = -170.0 + (i as f64) * 5.3;
        let lat2 = 80.0 - (i as f64) * 2.3;
        let lon2 = 170.0 - (i as f64) * 6.1;
        pairs.push((
            g(lat1.clamp(-90.0, 90.0), lon1),
            g(lat2.clamp(-90.0, 90.0), lon2),
        ));
    }
    pairs
}

/// |a − b| with degree wrap-around at 360.
fn ang_close(a: f64, b: f64, tol: f64) -> bool {
    let d = (a - b).rem_euclid(360.0);
    d.min(360.0 - d) <= tol
}

#[test]
fn inverse_matches_geographiclib_on_wgs84() {
    let oracle = OracleGeodesic::wgs84();
    for (a, b) in point_pairs() {
        let ours = Ellipsoid::WGS84.inverse(a, b);
        let (s12, azi1, azi2, _a12): (f64, f64, f64, f64) =
            oracle.inverse(a.lat, a.lon, b.lat, b.lon);
        assert!(
            (ours.distance_m - s12).abs() < 1e-8,
            "distance {} vs {} for {a:?} → {b:?}",
            ours.distance_m,
            s12
        );
        // Skip azimuth comparison for coincident points (it is arbitrary).
        if s12 > 1e-9 {
            assert!(
                ang_close(ours.initial_bearing.degrees(), azi1, 1e-9),
                "azi1 {} vs {azi1} for {a:?} → {b:?}",
                ours.initial_bearing.degrees()
            );
            assert!(
                ang_close(ours.final_bearing.degrees(), azi2, 1e-9),
                "azi2 {} vs {azi2} for {a:?} → {b:?}",
                ours.final_bearing.degrees()
            );
        }
    }
}

#[test]
fn direct_matches_geographiclib_on_wgs84() {
    let oracle = OracleGeodesic::wgs84();
    for (i, (origin, _)) in point_pairs().into_iter().enumerate() {
        let bearing = (i as f64) * 23.7 - 180.0;
        let distance = 1.0 + (i as f64) * 330_777.0; // 1 m … ~20,000 km
        let ours = Ellipsoid::WGS84.direct(origin, Bearing::from_degrees(bearing), distance);
        let (lat2, lon2, azi2): (f64, f64, f64) =
            oracle.direct(origin.lat, origin.lon, bearing, distance);
        assert!(
            (ours.point.lat - lat2).abs() < 1e-12,
            "lat {} vs {lat2} from {origin:?} az {bearing} d {distance}",
            ours.point.lat
        );
        assert!(
            ang_close(ours.point.lon, lon2, 1e-11),
            "lon {} vs {lon2} from {origin:?} az {bearing} d {distance}",
            ours.point.lon
        );
        assert!(
            ang_close(ours.final_bearing.degrees(), azi2, 1e-9),
            "azi2 {} vs {azi2} from {origin:?} az {bearing} d {distance}",
            ours.final_bearing.degrees()
        );
    }
}

#[test]
fn non_wgs84_ellipsoids_match_geographiclib() {
    // The uncached path: GRS80, a sphere, and Mars' much larger flattening.
    for e in [Ellipsoid::GRS80, Ellipsoid::SPHERE, Ellipsoid::MARS] {
        let oracle = OracleGeodesic::new(e.a, e.f);
        for (a, b) in point_pairs().into_iter().step_by(7) {
            let ours = e.inverse(a, b);
            let s12: f64 = oracle.inverse(a.lat, a.lon, b.lat, b.lon);
            assert!(
                (ours.distance_m - s12).abs() < 1e-8,
                "distance {} vs {s12} on {e:?} for {a:?} → {b:?}",
                ours.distance_m
            );
        }
    }
}

#[test]
fn round_trip_inverse_then_direct() {
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
fn r2_geodesic_metric_matches_geo() {
    // The r2 metric struct (x = lon, y = lat) against geo's Geodesic, which
    // wraps the same geographiclib algorithms.
    use geo::{Bearing as _, Destination as _, Distance as _};
    use r2::{Bearing as _, Destination as _, Distance as _, Length as _};

    for (a, b) in point_pairs().into_iter().step_by(5) {
        let ours =
            r2::Geodesic.distance(r2::Point::new(a.lon, a.lat), r2::Point::new(b.lon, b.lat));
        let geos =
            geo::Geodesic.distance(geo::Point::new(a.lon, a.lat), geo::Point::new(b.lon, b.lat));
        assert!(
            (ours - geos).abs() < 1e-8,
            "{ours} vs {geos} for {a:?} → {b:?}"
        );
    }

    let (a, b) = (g(40.6398, -73.7787), g(51.4775, -0.4614));
    let ours = r2::Geodesic.bearing(r2::Point::new(a.lon, a.lat), r2::Point::new(b.lon, b.lat));
    let geos = geo::Geodesic.bearing(geo::Point::new(a.lon, a.lat), geo::Point::new(b.lon, b.lat));
    assert!(
        ang_close(ours, geos.rem_euclid(360.0), 1e-9),
        "{ours} vs {geos}"
    );

    let ours = r2::Geodesic.destination(r2::Point::new(a.lon, a.lat), 51.0, 5_551_759.4);
    let geos = geo::Geodesic.destination(geo::Point::new(a.lon, a.lat), 51.0, 5_551_759.4);
    assert!((ours.x() - geos.x()).abs() < 1e-11 && (ours.y() - geos.y()).abs() < 1e-11);

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

#[test]
fn geodesic_area_matches_geo() {
    use geo::GeodesicArea as _;
    use r2::GeodesicArea as _;

    // CCW square with a CW hole, a CW-wound exterior, and a high-latitude
    // sliver crossing the antimeridian.
    let cases: Vec<(r2::Polygon, geo::Polygon)> = vec![
        (
            r2::Polygon::new(
                r2::LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
                vec![r2::LineString::from(vec![
                    (0.25, 0.25),
                    (0.25, 0.75),
                    (0.75, 0.75),
                    (0.75, 0.25),
                ])],
            ),
            geo::Polygon::new(
                geo::LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
                vec![geo::LineString::from(vec![
                    (0.25, 0.25),
                    (0.25, 0.75),
                    (0.75, 0.75),
                    (0.75, 0.25),
                ])],
            ),
        ),
        (
            // Clockwise exterior: signed area must come back negative.
            r2::Polygon::new(
                r2::LineString::from(vec![(2.0, 48.0), (2.0, 49.0), (3.0, 49.0), (3.0, 48.0)]),
                vec![],
            ),
            geo::Polygon::new(
                geo::LineString::from(vec![(2.0, 48.0), (2.0, 49.0), (3.0, 49.0), (3.0, 48.0)]),
                vec![],
            ),
        ),
        (
            r2::Polygon::new(
                r2::LineString::from(vec![
                    (179.5, 70.0),
                    (-179.5, 70.0),
                    (-179.5, 71.0),
                    (179.5, 71.0),
                ]),
                vec![],
            ),
            geo::Polygon::new(
                geo::LineString::from(vec![
                    (179.5, 70.0),
                    (-179.5, 70.0),
                    (-179.5, 71.0),
                    (179.5, 71.0),
                ]),
                vec![],
            ),
        ),
    ];

    for (ours, theirs) in &cases {
        let (p_ours, a_ours) = ours.geodesic_perimeter_area_signed();
        let (p_geo, a_geo) = theirs.geodesic_perimeter_area_signed();
        assert!((p_ours - p_geo).abs() < 1e-6, "{p_ours} vs {p_geo}");
        assert!((a_ours - a_geo).abs() < 1e-2, "{a_ours} vs {a_geo}");

        let (pu_ours, au_ours) = ours.geodesic_perimeter_area_unsigned();
        let (pu_geo, au_geo) = theirs.geodesic_perimeter_area_unsigned();
        assert!((pu_ours - pu_geo).abs() < 1e-6, "{pu_ours} vs {pu_geo}");
        assert!((au_ours - au_geo).abs() < 1e-2, "{au_ours} vs {au_geo}");
    }
}

#[test]
fn geodesic_area_improves_on_spherical() {
    use r2::{ChamberlainDuquetteArea as _, GeodesicArea as _};

    // A 1°×1° square at mid-latitude: the spherical approximation should sit
    // within ~0.6% of the exact ellipsoidal area, and they must disagree (the
    // whole point of the Karney area).
    let poly = r2::Polygon::new(
        r2::LineString::from(vec![(2.0, 48.0), (3.0, 48.0), (3.0, 49.0), (2.0, 49.0)]),
        vec![],
    );
    let exact = poly.geodesic_area_unsigned();
    let spherical = poly.chamberlain_duquette_unsigned_area();
    let rel = (exact - spherical).abs() / exact;
    assert!(rel > 1e-5 && rel < 6e-3, "relative gap {rel}");

    // Ellipsoid-level API agrees with the r2 trait.
    let ring: Vec<GeoPoint> = [(48.0, 2.0), (48.0, 3.0), (49.0, 3.0), (49.0, 2.0)]
        .iter()
        .map(|&(lat, lon)| g(lat, lon))
        .collect();
    let (_, area) = Ellipsoid::WGS84.ring_perimeter_area_m(&ring);
    assert!((area - exact).abs() < 1e-3, "{area} vs {exact}");

    // The whole planet, both hemispheres: total area sanity.
    let total = Ellipsoid::WGS84.total_area_m2();
    assert!((total - 510.0656e12).abs() / total < 1e-4, "{total}");
}

#[test]
fn ellipsoidal_vs_spherical_sanity() {
    // The ellipsoidal distance differs from haversine by less than ~0.6%
    // everywhere, and equals it (to round-off scale) on the f = 0 ellipsoid.
    for (a, b) in point_pairs() {
        let ell = Ellipsoid::WGS84.distance_m(a, b);
        let sph = haversine_m(a, b);
        if sph > 1_000.0 {
            assert!(
                (ell - sph).abs() / sph < 6e-3,
                "ellipsoidal {ell} vs spherical {sph} for {a:?} → {b:?}"
            );
        }
    }
}
