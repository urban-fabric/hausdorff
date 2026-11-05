//! Rhumb-line (loxodrome) metric conformance, oracle-free. Correctness is pinned
//! by: round-trip consistency (destination↔distance↔bearing), facade-vs-kernel
//! agreement (the r2 `Rhumb` metric, the typed and `_m` kernels), analytic
//! relations (equator == meridian quarter; east-west distance scales as cos φ;
//! linear in Δlon at the equator; cardinal bearings), proportional
//! point-at-ratio, length = segment sum, and densify path-length preservation,
//! plus a few golden absolute-scale anchors. No third-party geometry crate.

use crate::space::r2::{Bearing as _, Densify as _, Destination as _, Distance as _, Length as _};
use hausdorff_base::earth::{
    Bearing as CompassBearing, GeoPoint, LengthMeters, rhumb_bearing, rhumb_destination,
    rhumb_destination_m, rhumb_distance_m, rhumb_distance_meters, rhumb_point_at_ratio,
};

const REL_EPS: f64 = 1e-6;
const DEG_EPS: f64 = 1e-7;

// Golden absolute-scale anchors (captured from the validated kernel == geo).
const EQ_QUARTER: f64 = 10007557.221017962; // equator/meridian quarter circumference
const NYC: (f64, f64) = (-74.006, 40.7128);
const LON: (f64, f64) = (-0.1278, 51.5074);
const NYC_LON_DIST: f64 = 5794128.753915687;
const NYC_LON_BRG: f64 = 78.0440808681045;

fn gp(lon: f64, lat: f64) -> GeoPoint {
    GeoPoint::from_lon_lat(lon, lat)
}
fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    rhumb_distance_m(gp(a.0, a.1), gp(b.0, b.1))
}
fn brg(a: (f64, f64), b: (f64, f64)) -> f64 {
    rhumb_bearing(gp(a.0, a.1), gp(b.0, b.1)).degrees()
}

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
        ((-74.006, 40.7128), (-0.1278, 51.5074)),
        ((151.21, -33.87), (139.69, 35.68)),
    ];
    let grid = [-60.0, -30.0, 0.0, 30.0, 60.0];
    let lons = [-150.0, -60.0, 0.0, 60.0, 150.0];
    for (i, &lat1) in grid.iter().enumerate() {
        for &lat2 in grid.iter() {
            pairs.push(((lons[i], lat1), (lons[i], lat2)));
        }
    }
    pairs
}

fn rel_close(ours: f64, want: f64, what: &str) {
    let scale = want.abs().max(1.0);
    assert!(
        (ours - want).abs() / scale < REL_EPS,
        "{what}: {ours} vs {want}"
    );
}
fn deg_close(ours: f64, want: f64, what: &str) {
    let diff = ((ours - want + 540.0).rem_euclid(360.0) - 180.0).abs();
    assert!(diff < DEG_EPS, "{what}: {ours} vs {want}");
}

#[test]
fn distance_absolute_and_analytic() {
    // Golden absolute anchors.
    rel_close(dist((0.0, 0.0), (90.0, 0.0)), EQ_QUARTER, "equator quarter");
    rel_close(
        dist((0.0, 0.0), (0.0, 90.0)),
        EQ_QUARTER,
        "meridian quarter",
    );
    rel_close(dist(NYC, LON), NYC_LON_DIST, "NYC→London");

    // Linear in Δlon along the equator.
    for d in [10.0, 30.0, 45.0, 120.0] {
        rel_close(
            dist((0.0, 0.0), (d, 0.0)),
            EQ_QUARTER * d / 90.0,
            "equator linear",
        );
    }
    // East-west distance at latitude φ scales as cos φ.
    for lat in [15.0, 30.0, 45.0, 60.0] {
        rel_close(
            dist((0.0, lat), (90.0, lat)),
            EQ_QUARTER * lat.to_radians().cos(),
            "east-west cos(lat)",
        );
    }
}

#[test]
fn bearing_cardinals_and_anchor() {
    deg_close(brg((0.0, 0.0), (10.0, 0.0)), 90.0, "due east");
    deg_close(brg((0.0, 0.0), (0.0, 10.0)), 0.0, "due north");
    deg_close(brg((0.0, 0.0), (-10.0, 0.0)), 270.0, "due west");
    deg_close(brg((0.0, 10.0), (0.0, 0.0)), 180.0, "due south");
    deg_close(brg(NYC, LON), NYC_LON_BRG, "NYC→London bearing");
}

#[test]
fn destination_round_trips_distance_and_bearing() {
    // Mid-latitude origins only: a due-north course from near a pole overshoots
    // it (a rhumb singularity), which is out of scope for the round-trip.
    let origins = [
        (13.405, 52.52),
        (179.8, 10.0),
        (-179.9, -45.0),
        (0.0, 0.0),
        (-100.0, -33.0),
    ];
    let courses = [
        (0.0, 100_000.0),
        (45.0, 2_500_000.0),
        (90.0, 5_000_000.0),
        (215.0, 750.0),
        (359.0, 12_000.0),
    ];
    for &(x, y) in &origins {
        for &(bearing, d) in &courses {
            let dest = rhumb_destination_m(gp(x, y), CompassBearing::from_degrees(bearing), d);
            rel_close(dist((x, y), (dest.lon, dest.lat)), d, "round-trip distance");
            deg_close(
                brg((x, y), (dest.lon, dest.lat)),
                bearing,
                "round-trip bearing",
            );
        }
    }
    // Golden destination anchor.
    let dest = rhumb_destination_m(
        gp(13.405, 52.52),
        CompassBearing::from_degrees(45.0),
        2_500_000.0,
    );
    deg_close(dest.lon, 46.44157217795885, "dest lon");
    assert!(
        (dest.lat - 68.41788819121933).abs() < DEG_EPS,
        "dest lat {}",
        dest.lat
    );
}

#[test]
fn facade_and_typed_kernels_agree() {
    for &((x1, y1), (x2, y2)) in &sample_pairs() {
        let kernel = rhumb_distance_m(gp(x1, y1), gp(x2, y2));
        let metric = crate::space::r2::Rhumb.distance(
            crate::space::r2::Point::new(x1, y1),
            crate::space::r2::Point::new(x2, y2),
        );
        rel_close(metric, kernel, "metric distance == kernel");
        rel_close(
            rhumb_distance_meters(gp(x1, y1), gp(x2, y2)).value(),
            kernel,
            "typed distance == kernel",
        );

        let metric_brg = crate::space::r2::Rhumb.bearing(
            crate::space::r2::Point::new(x1, y1),
            crate::space::r2::Point::new(x2, y2),
        );
        deg_close(
            metric_brg,
            brg((x1, y1), (x2, y2)),
            "metric bearing == kernel",
        );
    }
    // destination: metric == _m == typed.
    let (o, bearing, d) = ((13.405, 52.52), 215.0, 750.0);
    let km = rhumb_destination_m(gp(o.0, o.1), CompassBearing::from_degrees(bearing), d);
    let metric =
        crate::space::r2::Rhumb.destination(crate::space::r2::Point::new(o.0, o.1), bearing, d);
    deg_close(metric.x(), km.lon, "dest metric lon");
    assert!((metric.y() - km.lat).abs() < DEG_EPS, "dest metric lat");
    let typed = rhumb_destination(
        gp(o.0, o.1),
        CompassBearing::from_degrees(bearing),
        LengthMeters::new(d),
    );
    assert_eq!((typed.lat, typed.lon), (km.lat, km.lon), "typed == _m");
}

#[test]
fn point_at_ratio_is_proportional_and_collinear() {
    for &((x1, y1), (x2, y2)) in &sample_pairs() {
        let total = dist((x1, y1), (x2, y2));
        let a = rhumb_point_at_ratio(gp(x1, y1), gp(x2, y2), 0.0);
        let b = rhumb_point_at_ratio(gp(x1, y1), gp(x2, y2), 1.0);
        deg_close(a.lon, x1, "ratio 0 lon");
        assert!((a.lat - y1).abs() < DEG_EPS, "ratio 0 lat");
        deg_close(b.lon, x2, "ratio 1 lon");
        assert!((b.lat - y2).abs() < DEG_EPS, "ratio 1 lat");
        if total < 1.0 {
            continue;
        }
        for ratio in [0.1, 0.5, 0.9] {
            let p = rhumb_point_at_ratio(gp(x1, y1), gp(x2, y2), ratio);
            rel_close(
                dist((x1, y1), (p.lon, p.lat)),
                ratio * total,
                "proportional distance",
            );
            deg_close(
                brg((x1, y1), (p.lon, p.lat)),
                brg((x1, y1), (x2, y2)),
                "collinear bearing",
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

fn h_line(coords: &[(f64, f64)]) -> crate::space::r2::LineString {
    crate::space::r2::LineString::from(coords.to_vec())
}

#[test]
fn length_equals_segment_sum() {
    for coords in fixture_lines() {
        let total = crate::space::r2::Rhumb.length(&h_line(&coords));
        let sum: f64 = coords.windows(2).map(|w| dist(w[0], w[1])).sum();
        rel_close(total, sum, &format!("length of {coords:?}"));
    }
}

#[test]
fn densify_preserves_path_endpoints_and_caps_gaps() {
    for coords in fixture_lines() {
        let original = h_line(&coords);
        let orig_len = crate::space::r2::Rhumb.length(&original);
        for max_dist in [250_000.0, 1_000_000.0, 50_000_000.0] {
            let dense = crate::space::r2::Rhumb.densify(&original, max_dist);
            let pts: Vec<(f64, f64)> = dense.coords().map(|c| (c.x, c.y)).collect();
            assert!(pts.len() >= coords.len(), "densify never drops vertices");
            // Endpoints preserved.
            deg_close(pts[0].0, coords[0].0, "densify start lon");
            deg_close(
                pts.last().unwrap().0,
                coords.last().unwrap().0,
                "densify end lon",
            );
            // No gap exceeds max_dist, and the rhumb path length is unchanged.
            for w in pts.windows(2) {
                assert!(dist(w[0], w[1]) <= max_dist + 1e-3, "gap exceeds max");
            }
            rel_close(
                crate::space::r2::Rhumb.length(&dense),
                orig_len,
                "densify preserves length",
            );
        }
    }
}

#[test]
fn euclidean_densify_caps_gaps_and_stays_collinear() {
    let coords = vec![(0.0, 0.0), (0.0, 6.0), (1.0, 7.0)];
    let line = h_line(&coords);
    let dense = crate::space::r2::Euclidean.densify(&line, 2.0);
    let pts: Vec<(f64, f64)> = dense.coords().map(|c| (c.x, c.y)).collect();
    assert_eq!(pts[0], (0.0, 0.0));
    assert_eq!(*pts.last().unwrap(), (1.0, 7.0));
    for w in pts.windows(2) {
        let gap = ((w[1].0 - w[0].0).powi(2) + (w[1].1 - w[0].1).powi(2)).sqrt();
        assert!(gap <= 2.0 + 1e-9, "euclidean gap exceeds max");
    }
    // Euclidean length is preserved (densified points lie on the segments).
    let seg_len = |a: (f64, f64), b: (f64, f64)| ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
    let orig: f64 = coords.windows(2).map(|w| seg_len(w[0], w[1])).sum();
    let dense_len: f64 = pts.windows(2).map(|w| seg_len(w[0], w[1])).sum();
    assert!(
        (orig - dense_len).abs() < 1e-9,
        "euclidean densify preserves length"
    );
}
