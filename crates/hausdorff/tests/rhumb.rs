//! Conformance tests for the rhumb-line (loxodrome) metric: cross-check the
//! in-house kernels and the r2 `Rhumb` metric space against georust `geo`'s
//! `Rhumb` (a dev-dependency oracle). The same inputs are fed to both
//! libraries and the outputs must agree to tight tolerance.

use geo::{
    Bearing as _, Densify as _, Destination as _, Distance as _, InterpolatePoint as _,
    Length as _, Rhumb as GeoRhumb,
};

use euclid::LengthMeters;
use hausdorff::earth::{
    Bearing as CompassBearing, GeoPoint, rhumb_bearing, rhumb_destination, rhumb_destination_m,
    rhumb_distance_m, rhumb_distance_meters, rhumb_point_at_ratio,
};
use hausdorff::{Bearing as _, Densify as _, Destination as _, Distance as _, Length as _};

const REL_EPS: f64 = 1e-9;
const DEG_EPS: f64 = 1e-9;

/// Point pairs as `((lon, lat), (lon, lat))`: a coarse global grid plus the
/// singular cases — antimeridian crossings, near-pole runs, east-west
/// (constant latitude) courses where the isometric-latitude delta vanishes,
/// meridians, and coincident points.
fn sample_pairs() -> Vec<((f64, f64), (f64, f64))> {
    let mut pairs = vec![
        ((179.5, 10.0), (-179.5, 12.0)),
        ((170.0, -20.0), (-170.0, -25.0)),
        ((-179.9, 60.0), (179.9, 60.0)),
        ((150.0, 0.0), (-160.0, 5.0)),
        ((30.0, 89.5), (60.0, 89.7)),
        ((-120.0, -89.0), (100.0, -89.9)),
        ((0.0, 89.9), (180.0, 89.9)),
        ((10.0, 45.0), (50.0, 45.0)),
        ((-100.0, -33.0), (20.0, -33.0)),
        ((0.0, 0.0), (90.0, 0.0)),
        ((179.0, 70.0), (-178.0, 70.0)),
        ((25.0, 10.0), (25.0, 60.0)),
        ((5.0, 5.0), (5.0, 5.0)),
        ((-74.006, 40.7128), (-0.1278, 51.5074)),
        ((151.21, -33.87), (139.69, 35.68)),
    ];
    let grid = [-60.0, -30.0, 0.0, 30.0, 60.0];
    let lons = [-150.0, -60.0, 0.0, 60.0, 150.0];
    for (i, &lat1) in grid.iter().enumerate() {
        for (j, &lat2) in grid.iter().enumerate() {
            pairs.push(((lons[i], lat1), (lons[j], lat2)));
        }
    }
    pairs
}

fn assert_close_rel(ours: f64, theirs: f64, what: &str) {
    let scale = theirs.abs().max(1.0);
    assert!(
        (ours - theirs).abs() / scale < REL_EPS,
        "{what}: {ours} vs {theirs}"
    );
}

fn assert_close_deg(ours: f64, theirs: f64, what: &str) {
    let diff = ((ours - theirs + 540.0).rem_euclid(360.0) - 180.0).abs();
    assert!(diff < DEG_EPS, "{what}: {ours} vs {theirs}");
}

#[test]
fn distance_matches_geo() {
    for &((x1, y1), (x2, y2)) in &sample_pairs() {
        let ours = rhumb_distance_m(
            GeoPoint::from_lon_lat(x1, y1),
            GeoPoint::from_lon_lat(x2, y2),
        );
        let theirs = GeoRhumb.distance(geo::Point::new(x1, y1), geo::Point::new(x2, y2));
        assert_close_rel(ours, theirs, &format!("distance ({x1},{y1})->({x2},{y2})"));

        let metric =
            hausdorff::Rhumb.distance(hausdorff::Point::new(x1, y1), hausdorff::Point::new(x2, y2));
        assert_close_rel(metric, theirs, "metric-space distance");
        assert_close_rel(
            rhumb_distance_meters(
                GeoPoint::from_lon_lat(x1, y1),
                GeoPoint::from_lon_lat(x2, y2),
            )
            .as_meters(),
            theirs,
            "typed distance",
        );
    }
}

#[test]
fn bearing_matches_geo() {
    for &((x1, y1), (x2, y2)) in &sample_pairs() {
        let ours = rhumb_bearing(
            GeoPoint::from_lon_lat(x1, y1),
            GeoPoint::from_lon_lat(x2, y2),
        );
        let theirs = GeoRhumb.bearing(geo::Point::new(x1, y1), geo::Point::new(x2, y2));
        assert_close_deg(
            ours.degrees(),
            theirs,
            &format!("bearing ({x1},{y1})->({x2},{y2})"),
        );

        let metric =
            hausdorff::Rhumb.bearing(hausdorff::Point::new(x1, y1), hausdorff::Point::new(x2, y2));
        assert_close_deg(metric, theirs, "metric-space bearing");
    }
}

#[test]
fn destination_matches_geo() {
    let origins = [
        (13.405, 52.52),
        (179.8, 10.0),
        (-179.9, -45.0),
        (0.0, 0.0),
        (30.0, 89.5),
        (-60.0, -89.5),
    ];
    let courses = [
        (0.0, 100_000.0),
        (45.0, 2_500_000.0),
        (90.0, 5_000_000.0),
        (90.0, 50_000.0),
        (180.0, 1_000_000.0),
        (270.0, 3_000_000.0),
        (215.0, 750.0),
        (359.0, 12_000.0),
    ];
    for &(x, y) in &origins {
        for &(bearing, dist) in &courses {
            let ours = rhumb_destination_m(
                GeoPoint::from_lon_lat(x, y),
                CompassBearing::from_degrees(bearing),
                dist,
            );
            let theirs = GeoRhumb.destination(geo::Point::new(x, y), bearing, dist);
            let what = format!("destination from ({x},{y}) bearing {bearing} dist {dist}");
            assert_close_deg(ours.lon, theirs.x(), &format!("{what}: lon"));
            assert!(
                (ours.lat - theirs.y()).abs() < DEG_EPS,
                "{what}: lat {} vs {}",
                ours.lat,
                theirs.y()
            );

            let metric = hausdorff::Rhumb.destination(hausdorff::Point::new(x, y), bearing, dist);
            assert_close_deg(metric.x(), theirs.x(), &format!("{what}: metric lon"));
            assert!(
                (metric.y() - theirs.y()).abs() < DEG_EPS,
                "{what}: metric lat"
            );

            let typed = rhumb_destination(
                GeoPoint::from_lon_lat(x, y),
                CompassBearing::from_degrees(bearing),
                LengthMeters::new(dist),
            );
            assert_eq!(
                (typed.lat, typed.lon),
                (ours.lat, ours.lon),
                "{what}: typed"
            );
        }
    }
}

#[test]
fn destination_is_inverse_of_distance_and_bearing() {
    let origin = GeoPoint::from_lon_lat(13.405, 52.52);
    for deg in [0.0, 37.0, 90.0, 180.0, 271.5, 359.0] {
        let dest = rhumb_destination_m(origin, CompassBearing::from_degrees(deg), 1_500.0);
        assert!(
            (rhumb_distance_m(origin, dest) - 1_500.0).abs() < 1e-6,
            "deg {deg}"
        );
        assert_close_deg(
            rhumb_bearing(origin, dest).degrees(),
            deg,
            &format!("round-trip bearing {deg}"),
        );
    }
}

#[test]
fn point_at_ratio_matches_geo() {
    for &((x1, y1), (x2, y2)) in &sample_pairs() {
        for ratio in [0.0, 0.1, 0.5, 0.9, 1.0] {
            let ours = rhumb_point_at_ratio(
                GeoPoint::from_lon_lat(x1, y1),
                GeoPoint::from_lon_lat(x2, y2),
                ratio,
            );
            let theirs = GeoRhumb.point_at_ratio_between(
                geo::Point::new(x1, y1),
                geo::Point::new(x2, y2),
                ratio,
            );
            let what = format!("ratio {ratio} on ({x1},{y1})->({x2},{y2})");
            assert_close_deg(ours.lon, theirs.x(), &format!("{what}: lon"));
            assert!(
                (ours.lat - theirs.y()).abs() < DEG_EPS,
                "{what}: lat {} vs {}",
                ours.lat,
                theirs.y()
            );
        }
    }
}

fn fixture_lines() -> Vec<Vec<(f64, f64)>> {
    vec![
        vec![(13.405, 52.52), (2.3522, 48.8566), (-0.1278, 51.5074)],
        vec![
            (170.0, -20.0),
            (179.5, -22.0),
            (-178.0, -23.5),
            (-160.0, -30.0),
        ],
        vec![(10.0, 45.0), (50.0, 45.0), (90.0, 45.0)],
        vec![(30.0, 88.0), (60.0, 89.0), (120.0, 89.5)],
        vec![(0.0, 0.0), (90.0, 0.0), (90.0, 45.0)],
    ]
}

fn h_line(coords: &[(f64, f64)]) -> hausdorff::LineString {
    hausdorff::LineString::from(coords.to_vec())
}

fn g_line(coords: &[(f64, f64)]) -> geo::LineString<f64> {
    geo::LineString(coords.iter().map(|&(x, y)| geo::Coord { x, y }).collect())
}

#[test]
fn line_string_length_matches_geo() {
    for coords in fixture_lines() {
        let ours = hausdorff::Rhumb.length(&h_line(&coords));
        let theirs = GeoRhumb.length(&g_line(&coords));
        assert_close_rel(ours, theirs, &format!("length of {coords:?}"));
    }
}

#[test]
fn densify_matches_geo_vertex_by_vertex() {
    for coords in fixture_lines() {
        for max_dist in [250_000.0, 1_000_000.0, 50_000_000.0] {
            let ours = hausdorff::Rhumb.densify(&h_line(&coords), max_dist);
            let theirs = GeoRhumb.densify(&g_line(&coords), max_dist);
            assert_eq!(
                ours.coords_count(),
                theirs.0.len(),
                "vertex count on {coords:?} max {max_dist}"
            );
            for (i, (h, g)) in ours.coords().zip(theirs.0.iter()).enumerate() {
                let what = format!("densify vertex {i} of {coords:?} max {max_dist}");
                assert_close_deg(h.x, g.x, &format!("{what}: lon"));
                assert!(
                    (h.y - g.y).abs() < DEG_EPS,
                    "{what}: lat {} vs {}",
                    h.y,
                    g.y
                );
            }
        }
    }
}

#[test]
fn euclidean_densify_matches_geo() {
    let coords = vec![(0.0, 0.0), (0.0, 6.0), (1.0, 7.0)];
    let ours = hausdorff::Euclidean.densify(&h_line(&coords), 2.0);
    let theirs = geo::Euclidean.densify(&g_line(&coords), 2.0);
    assert_eq!(ours.coords_count(), theirs.0.len());
    for (h, g) in ours.coords().zip(theirs.0.iter()) {
        assert!((h.x - g.x).abs() < 1e-12 && (h.y - g.y).abs() < 1e-12);
    }
}
