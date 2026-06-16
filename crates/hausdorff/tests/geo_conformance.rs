//! Conformance for the in-house planar/spherical algorithms, oracle-free.
//! Structural results are checked against independent dependency-free references
//! (shoelace area, min/max bounding box, ray-cast point-in-polygon, RDP
//! post-conditions); spherical/curve results without a simple closed form are
//! anchored by golden values captured from the validated implementation. No
//! third-party geometry crate.

use hausdorff::earth::{GeoPoint, frechet_distance_m, haversine_m};
use hausdorff::{
    Area, BoundingRect, Centroid, ChamberlainDuquetteArea, Contains, ConvexHull, CoordsIter,
    Simplify,
};
use hausdorff::{Bearing as _, Destination as _, Distance as _, InterpolatePoint as _, Length as _};

const EPS: f64 = 1e-9;

fn h_poly(coords: &[(f64, f64)]) -> hausdorff::Polygon {
    hausdorff::Polygon::new(
        hausdorff::LineString::new(coords.iter().map(|&(x, y)| hausdorff::Coord::new(x, y)).collect()),
        vec![],
    )
}

fn fixtures() -> Vec<Vec<(f64, f64)>> {
    vec![
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        vec![(0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0)],
        vec![(0.0, 0.0), (6.0, 0.0), (0.0, 3.0)],
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 4.0), (0.0, 4.0)],
        vec![(1.0, -3.1), (3.7, 2.7), (0.9, 7.6), (-4.8, 6.7), (-7.5, 0.9)],
    ]
}

/// Independent signed area (shoelace) over the closed ring.
fn shoelace(coords: &[(f64, f64)]) -> f64 {
    let n = coords.len();
    let mut s = 0.0;
    for i in 0..n {
        let (x0, y0) = coords[i];
        let (x1, y1) = coords[(i + 1) % n];
        s += x0 * y1 - x1 * y0;
    }
    s * 0.5
}

/// Independent ray-casting point-in-polygon (interior; boundary is ambiguous).
fn ray_cast(coords: &[(f64, f64)], x: f64, y: f64) -> bool {
    let n = coords.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = coords[i];
        let (xj, yj) = coords[j];
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

#[test]
fn area_matches_shoelace() {
    for f in fixtures() {
        let hp = h_poly(&f);
        let signed = shoelace(&f);
        assert!((hp.signed_area() - signed).abs() < EPS, "signed_area on {f:?}");
        assert!((hp.unsigned_area() - signed.abs()).abs() < EPS, "unsigned_area on {f:?}");
    }
}

#[test]
fn centroid_matches_golden() {
    const GOLDEN: [(f64, f64); 5] = [
        (2.0, 2.0),
        (2.0, 2.0),
        (2.0, 1.0),
        (1.6666666666666667, 1.6666666666666667),
        (-1.5294390972452703, 2.58154663126452),
    ];
    for (f, &(gx, gy)) in fixtures().iter().zip(GOLDEN.iter()) {
        let c = h_poly(f).centroid().unwrap();
        assert!((c.x() - gx).abs() < EPS && (c.y() - gy).abs() < EPS, "centroid on {f:?}");
    }
}

#[test]
fn bounding_rect_matches_minmax() {
    for f in fixtures() {
        let r = h_poly(&f).bounding_rect().unwrap();
        let min_x = f.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
        let min_y = f.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let max_x = f.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
        let max_y = f.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
        assert!((r.min().x - min_x).abs() < EPS && (r.min().y - min_y).abs() < EPS);
        assert!((r.max().x - max_x).abs() < EPS && (r.max().y - max_y).abs() < EPS);
    }
}

#[test]
fn contains_matches_ray_cast() {
    for f in fixtures() {
        let hp = h_poly(&f);
        for i in -10..=10 {
            for j in -10..=10 {
                let (x, y) = (i as f64 * 0.9, j as f64 * 0.9);
                let hc = hp.contains(&hausdorff::Coord::new(x, y));
                let rc = ray_cast(&f, x, y);
                if hc != rc {
                    // Boundary points are ambiguous for ray casting; only tolerate there.
                    assert!(
                        is_near_boundary(&f, x, y, 1e-6),
                        "contains mismatch at ({x},{y}) on {f:?}: hausdorff={hc} raycast={rc}"
                    );
                }
            }
        }
    }
}

#[test]
fn convex_hull_area_matches_golden() {
    let point_sets = vec![
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (2.0, 2.0)],
        vec![(0.0, 0.0), (1.0, 0.5), (2.0, 0.0), (2.0, 2.0), (1.0, 1.5), (0.0, 2.0), (1.0, 1.0)],
        vec![(-3.0, -1.0), (5.0, 0.0), (3.0, 6.0), (-2.0, 4.0), (0.0, 0.0), (1.0, 2.0)],
    ];
    const GOLDEN: [f64; 3] = [16.0, 4.0, 36.5];
    for (pts, &want) in point_sets.iter().zip(GOLDEN.iter()) {
        let a = hausdorff::MultiPoint::from(pts.clone()).convex_hull().unsigned_area();
        assert!((a - want).abs() < 1e-9, "convex hull area on {pts:?}: {a} vs {want}");
    }
}

#[test]
fn haversine_distance_matches_golden() {
    // samples as (lat, lon); golden meters captured from the validated kernel.
    let samples = [
        ((52.52, 13.405), (48.8566, 2.3522)),
        ((40.0, -74.0), (34.05, -118.24)),
        ((0.0, 0.0), (0.0, 1.0)),
        ((-33.87, 151.21), (35.68, 139.69)),
    ];
    const GOLDEN: [f64; 4] = [877464.5379215094, 3941829.06891843, 111195.0802335329, 7825755.604165367];
    for (((lat1, lon1), (lat2, lon2)), &want) in samples.iter().zip(GOLDEN.iter()) {
        let ours = haversine_m(GeoPoint::from_lat_lon(*lat1, *lon1), GeoPoint::from_lat_lon(*lat2, *lon2));
        assert!((ours - want).abs() <= 1e-9 * want, "haversine: {ours} vs {want}");
    }
}

#[test]
fn frechet_distance_matches_golden() {
    let path_a = [(52.52, 13.40), (52.40, 13.20), (52.30, 13.00)];
    let path_b = [(52.50, 13.42), (52.42, 13.25), (52.20, 13.10)];
    let ours = frechet_distance_m(
        &path_a.map(|(la, lo)| GeoPoint::from_lat_lon(la, lo)),
        &path_b.map(|(la, lo)| GeoPoint::from_lat_lon(la, lo)),
    )
    .expect("non-empty paths");
    assert!((ours - 13037.875025959464).abs() <= 1e-6, "frechet {ours}");
}

/// Distance from `(x, y)` to the nearest polygon edge is within `tol`.
fn is_near_boundary(coords: &[(f64, f64)], x: f64, y: f64, tol: f64) -> bool {
    let n = coords.len();
    for i in 0..n {
        let (ax, ay) = coords[i];
        let (bx, by) = coords[(i + 1) % n];
        let (dx, dy) = (bx - ax, by - ay);
        let len2 = dx * dx + dy * dy;
        let t = if len2 == 0.0 {
            0.0
        } else {
            (((x - ax) * dx + (y - ay) * dy) / len2).clamp(0.0, 1.0)
        };
        let (px, py) = (ax + t * dx, ay + t * dy);
        if ((x - px).powi(2) + (y - py).powi(2)).sqrt() < tol {
            return true;
        }
    }
    false
}

fn h_line(coords: &[(f64, f64)]) -> hausdorff::LineString {
    hausdorff::LineString::new(coords.iter().map(|&(x, y)| hausdorff::Coord::new(x, y)).collect())
}

#[test]
fn coords_count_closes_the_ring() {
    // hausdorff's polygon iterator yields the closing coordinate (len + 1).
    for f in fixtures() {
        assert_eq!(h_poly(&f).coords_count(), f.len() + 1, "coords_count on {f:?}");
    }
}

#[test]
fn simplify_satisfies_rdp_postconditions() {
    let wiggly = vec![
        (0.0, 0.0), (1.0, 0.1), (2.0, -0.1), (3.0, 5.0), (4.0, 0.05),
        (5.0, -0.02), (6.0, 0.0), (7.0, 8.0), (8.0, 0.0),
    ];
    let mut prev_len = wiggly.len() + 1;
    for eps in [0.2_f64, 1.0, 3.0] {
        let s: Vec<(f64, f64)> = h_line(&wiggly).simplify(eps).0.iter().map(|c| (c.x, c.y)).collect();
        // Endpoints preserved, subset of the original, monotone in eps.
        assert_eq!(s[0], wiggly[0], "first endpoint kept at eps={eps}");
        assert_eq!(*s.last().unwrap(), *wiggly.last().unwrap(), "last endpoint kept at eps={eps}");
        assert!(s.len() >= 2 && s.len() <= prev_len, "monotone vertex count at eps={eps}");
        prev_len = s.len();
        for &c in &s {
            assert!(wiggly.contains(&c), "kept point {c:?} is an original at eps={eps}");
        }
        // Every original point is within eps of the simplified polyline.
        let max_dev = wiggly
            .iter()
            .map(|&p| s.windows(2).map(|w| seg_dist(p, w[0], w[1])).fold(f64::INFINITY, f64::min))
            .fold(0.0_f64, f64::max);
        assert!(max_dev <= eps + 1e-9, "deviation {max_dev} > eps {eps}");
    }
}

#[test]
fn chamberlain_duquette_area_matches_golden() {
    const GOLDEN: [f64; 5] = [
        198111444408.08057,
        198111444408.08057,
        111477308082.54442,
        148613772731.42798,
        931358449518.4138,
    ];
    for (f, &want) in fixtures().iter().zip(GOLDEN.iter()) {
        let h = h_poly(f).chamberlain_duquette_unsigned_area();
        assert!((h - want).abs() <= 1e-9 * want.max(1.0), "CD area on {f:?}: {h} vs {want}");
    }

    let o2 = [
        (0.00388383, 51.501574), (0.00538587, 51.502278), (0.00553607, 51.503299),
        (0.00467777, 51.504181), (0.00327229, 51.504435), (0.00187754, 51.504168),
        (0.00087976, 51.503380), (0.00107288, 51.502324), (0.00185608, 51.501770),
    ];
    assert_eq!(h_poly(&o2).chamberlain_duquette_unsigned_area().round(), 78_478.0);
}

#[test]
fn haversine_distance_and_bearing_match_golden() {
    // samples as (lon, lat); golden (distance_m, bearing_deg).
    let samples = [
        ((13.405, 52.52), (2.3522, 48.8566)),
        ((-74.0, 40.0), (-118.24, 34.05)),
        ((0.0, 0.0), (1.0, 0.0)),
        ((151.21, -33.87), (139.69, 35.68)),
    ];
    const GOLDEN: [(f64, f64); 4] = [
        (877464.5379215094, 246.74342315057362),
        (3941829.06891843, 274.6844913710564),
        (111195.0802335329, 90.0),
        (7825755.604165367, 350.0829385305175),
    ];
    for (&((x1, y1), (x2, y2)), &(wd, wb)) in samples.iter().zip(GOLDEN.iter()) {
        let (a, b) = (hausdorff::Point::new(x1, y1), hausdorff::Point::new(x2, y2));
        let d = hausdorff::Haversine.distance(a, b);
        assert!((d - wd).abs() <= 1e-9 * wd, "distance {d} vs {wd}");
        let br = hausdorff::Haversine.bearing(a, b);
        assert!(((br - wb + 540.0) % 360.0 - 180.0).abs() < 1e-9, "bearing {br} vs {wb}");
    }
}

#[test]
fn lengths_match_reference() {
    let coords: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 0.5), (2.0, 2.0), (3.0, 1.0), (4.0, 3.0)];
    // Euclidean length: independent segment sum.
    let euclid_sum: f64 = coords
        .windows(2)
        .map(|w| ((w[1].0 - w[0].0).powi(2) + (w[1].1 - w[0].1).powi(2)).sqrt())
        .sum();
    assert!((hausdorff::Euclidean.length(&h_line(&coords)) - euclid_sum).abs() < EPS, "euclidean length");
    // Haversine length: golden.
    assert!((hausdorff::Haversine.length(&h_line(&coords)) - 730594.5476435976).abs() <= 1e-6, "haversine length");
}

#[test]
fn haversine_destination_matches_golden() {
    let origin = hausdorff::Point::new(13.405, 52.52);
    let courses = [(0.0, 1000.0), (90.0, 5000.0), (215.0, 250.0), (359.0, 12000.0)];
    const GOLDEN: [(f64, f64); 4] = [
        (13.404999999999973, 52.52899320363725),
        (13.478898358112474, 52.519976988255365),
        (13.402880770254797, 52.51815828078832),
        (13.401897084247594, 52.62790196665409),
    ];
    for (&(bearing, dist), &(gx, gy)) in courses.iter().zip(GOLDEN.iter()) {
        let d = hausdorff::Haversine.destination(origin, bearing, dist);
        assert!((d.x() - gx).abs() < 1e-9 && (d.y() - gy).abs() < 1e-9, "destination {bearing}/{dist}");
    }
}

#[test]
fn haversine_interpolation_follows_the_great_circle() {
    let start = hausdorff::Point::new(-122.4194, 37.7749);
    let end = hausdorff::Point::new(-73.9857, 40.7484);
    let total = hausdorff::Haversine.distance(start, end);

    let half = hausdorff::Haversine.point_at_ratio_between(start, end, 0.5);
    let first_half = hausdorff::Haversine.distance(start, half);
    let second_half = hausdorff::Haversine.distance(half, end);
    assert!((first_half - total * 0.5).abs() <= total * 1e-9);
    assert!((second_half - total * 0.5).abs() <= total * 1e-9);

    let waypoint = hausdorff::Haversine.point_at_distance_between(start, end, 100_000.0);
    let travelled = hausdorff::Haversine.distance(start, waypoint);
    assert!((travelled - 100_000.0).abs() <= 1e-4);
}

/// Distance from `p` to the bounded segment `a`–`b` (RDP's metric).
fn seg_dist(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (abx, aby) = (b.0 - a.0, b.1 - a.1);
    let len2 = abx * abx + aby * aby;
    let t = if len2 == 0.0 {
        0.0
    } else {
        (((p.0 - a.0) * abx + (p.1 - a.1) * aby) / len2).clamp(0.0, 1.0)
    };
    let (qx, qy) = (a.0 + abx * t, a.1 + aby * t);
    ((p.0 - qx).powi(2) + (p.1 - qy).powi(2)).sqrt()
}

#[test]
fn real_wkt_fixtures_conformance() {
    use hausdorff::format::wkt::{from_wkt, to_wkt};

    // Golden convex-hull areas captured from the validated implementation.
    for (name, wkt, hull_golden) in [
        ("poly1", include_str!("fixtures/poly1.wkt"), 154873.47276695012_f64),
        ("poly2", include_str!("fixtures/poly2.wkt"), 481365.4849054594_f64),
    ] {
        let g = from_wkt(wkt).unwrap_or_else(|e| panic!("parse {name}: {e:?}"));
        let coords: Vec<(f64, f64)> = match &g {
            hausdorff::Geometry::LineString(ls) => ls.0.iter().map(|c| (c.x, c.y)).collect(),
            other => panic!("{name}: expected LineString, got {other:?}"),
        };
        assert!(coords.len() > 50, "{name}: real shape ({} pts)", coords.len());

        match from_wkt(&to_wkt(&g)).expect("re-parse") {
            hausdorff::Geometry::LineString(ls) => assert_eq!(ls.0.len(), coords.len()),
            other => panic!("{name}: round-trip changed type to {other:?}"),
        }

        let ha = hausdorff::MultiPoint::from(coords.clone()).convex_hull().unsigned_area();
        assert!((ha - hull_golden).abs() <= hull_golden.abs() * 1e-9, "{name}: hull area {ha} vs {hull_golden}");

        let hls = hausdorff::LineString::from(coords.clone());
        let mut prev_len = coords.len();
        for eps in [0.5_f64, 2.0, 10.0] {
            let s: Vec<(f64, f64)> = hls.simplify(eps).0.iter().map(|c| (c.x, c.y)).collect();
            assert!(s.len() >= 2 && s.len() <= prev_len, "{name}: monotone at eps={eps}");
            prev_len = s.len();
            assert_eq!(s[0], coords[0], "{name}: first endpoint kept");
            assert_eq!(*s.last().unwrap(), *coords.last().unwrap(), "{name}: last endpoint kept");
            for &c in &s {
                assert!(coords.contains(&c), "{name}: kept point {c:?} is an original");
            }
            let max_dev = coords
                .iter()
                .map(|&p| s.windows(2).map(|w| seg_dist(p, w[0], w[1])).fold(f64::INFINITY, f64::min))
                .fold(0.0_f64, f64::max);
            assert!(max_dev <= eps + 1e-9, "{name}: simplify deviation {max_dev} > eps {eps}");
        }
    }
}
