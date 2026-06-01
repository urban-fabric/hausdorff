//! Conformance tests: cross-check the in-house planar/spherical algorithms
//! against georust `geo` (a dev-dependency oracle, exactly the role `h3o` plays
//! for the H3 core). The same geometry is built in both libraries and the
//! algorithm outputs must agree.

use geo::line_measures::FrechetDistance;
// Bring geo's algorithm traits into scope for method resolution without binding
// names that would collide with hausdorff's same-named traits.
use geo::{
    Area as _, BoundingRect as _, Centroid as _, Contains as _, ConvexHull as _, Distance as _,
    Haversine,
};

use hausdorff::geo::{GeoPoint, frechet_distance_m, haversine_m};
use hausdorff::{Area, BoundingRect, Centroid, Contains, ConvexHull};

const EPS: f64 = 1e-9;

// ---- builders: the same coordinates in each library ------------------------

fn h_poly(coords: &[(f64, f64)]) -> hausdorff::Polygon {
    hausdorff::Polygon::new(
        hausdorff::LineString::new(
            coords
                .iter()
                .map(|&(x, y)| hausdorff::Coord::new(x, y))
                .collect(),
        ),
        vec![],
    )
}

fn g_poly(coords: &[(f64, f64)]) -> geo::Polygon<f64> {
    geo::Polygon::new(
        geo::LineString(coords.iter().map(|&(x, y)| geo::Coord { x, y }).collect()),
        vec![],
    )
}

fn fixtures() -> Vec<Vec<(f64, f64)>> {
    vec![
        // unit-ish square (CCW)
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        // square wound CW
        vec![(0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0)],
        // right triangle
        vec![(0.0, 0.0), (6.0, 0.0), (0.0, 3.0)],
        // concave L-shape
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 4.0), (0.0, 4.0)],
        // irregular pentagon
        vec![(1.0, -3.1), (3.7, 2.7), (0.9, 7.6), (-4.8, 6.7), (-7.5, 0.9)],
    ]
}

// ---- planar algorithm conformance ------------------------------------------

#[test]
fn area_matches_geo() {
    for f in fixtures() {
        let (hp, gp) = (h_poly(&f), g_poly(&f));
        assert!(
            (hp.unsigned_area() - gp.unsigned_area()).abs() < EPS,
            "unsigned_area mismatch on {f:?}: {} vs {}",
            hp.unsigned_area(),
            gp.unsigned_area()
        );
        assert!(
            (hp.signed_area() - gp.signed_area()).abs() < EPS,
            "signed_area mismatch on {f:?}: {} vs {}",
            hp.signed_area(),
            gp.signed_area()
        );
    }
}

#[test]
fn centroid_matches_geo() {
    for f in fixtures() {
        let hc = h_poly(&f).centroid().unwrap();
        let gc = g_poly(&f).centroid().unwrap();
        assert!(
            (hc.x() - gc.x()).abs() < EPS && (hc.y() - gc.y()).abs() < EPS,
            "centroid mismatch on {f:?}: ({},{}) vs ({},{})",
            hc.x(),
            hc.y(),
            gc.x(),
            gc.y()
        );
    }
}

#[test]
fn bounding_rect_matches_geo() {
    for f in fixtures() {
        let hr = h_poly(&f).bounding_rect().unwrap();
        let gr = g_poly(&f).bounding_rect().unwrap();
        assert!((hr.min().x - gr.min().x).abs() < EPS);
        assert!((hr.min().y - gr.min().y).abs() < EPS);
        assert!((hr.max().x - gr.max().x).abs() < EPS);
        assert!((hr.max().y - gr.max().y).abs() < EPS);
    }
}

#[test]
fn contains_matches_geo() {
    for f in fixtures() {
        let (hp, gp) = (h_poly(&f), g_poly(&f));
        // sample a grid spanning the fixture's extent; skip near-boundary points
        // where boundary tie-breaking legitimately differs between libraries.
        for i in -10..=10 {
            for j in -10..=10 {
                let (x, y) = (i as f64 * 0.9, j as f64 * 0.9);
                let hc = hp.contains(&hausdorff::Coord::new(x, y));
                let gc = gp.contains(&geo::Coord { x, y });
                if hc != gc {
                    // tolerate disagreement only very close to an edge
                    let near_edge = is_near_boundary(&f, x, y, 1e-6);
                    assert!(
                        near_edge,
                        "contains mismatch at ({x},{y}) on {f:?}: hausdorff={hc} geo={gc}"
                    );
                }
            }
        }
    }
}

#[test]
fn convex_hull_area_matches_geo() {
    let point_sets = vec![
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (2.0, 2.0)],
        vec![(0.0, 0.0), (1.0, 0.5), (2.0, 0.0), (2.0, 2.0), (1.0, 1.5), (0.0, 2.0), (1.0, 1.0)],
        vec![(-3.0, -1.0), (5.0, 0.0), (3.0, 6.0), (-2.0, 4.0), (0.0, 0.0), (1.0, 2.0)],
    ];
    for pts in point_sets {
        let hmp = hausdorff::MultiPoint::from(pts.clone());
        let gmp: geo::MultiPoint<f64> =
            pts.iter().map(|&(x, y)| geo::Point::new(x, y)).collect();
        let ha = hmp.convex_hull().unsigned_area();
        let ga = gmp.convex_hull().unsigned_area();
        assert!(
            (ha - ga).abs() < 1e-9,
            "convex hull area mismatch on {pts:?}: {ha} vs {ga}"
        );
    }
}

// ---- spherical conformance --------------------------------------------------

#[test]
fn haversine_distance_matches_geo() {
    // geo's Haversine mean radius (6_371_008.7714 m) vs ours (6_371_008.8 m)
    // differ by ~0.03 m, so allow a small absolute tolerance.
    let samples = [
        ((52.52, 13.405), (48.8566, 2.3522)),  // Berlin -> Paris
        ((40.0, -74.0), (34.05, -118.24)),     // NYC -> LA
        ((0.0, 0.0), (0.0, 1.0)),
        ((-33.87, 151.21), (35.68, 139.69)),   // Sydney -> Tokyo
    ];
    for ((lat1, lon1), (lat2, lon2)) in samples {
        let ours = haversine_m(
            GeoPoint::from_lat_lon(lat1, lon1),
            GeoPoint::from_lat_lon(lat2, lon2),
        );
        let theirs = Haversine.distance(geo::Point::new(lon1, lat1), geo::Point::new(lon2, lat2));
        let rel = (ours - theirs).abs() / theirs.max(1.0);
        assert!(rel < 1e-5, "haversine mismatch: ours={ours} geo={theirs} rel={rel}");
    }
}

#[test]
fn frechet_distance_matches_geo() {
    let path_a = [(52.52, 13.40), (52.40, 13.20), (52.30, 13.00)];
    let path_b = [(52.50, 13.42), (52.42, 13.25), (52.20, 13.10)];

    let ours = frechet_distance_m(
        &path_a.map(|(la, lo)| GeoPoint::from_lat_lon(la, lo)),
        &path_b.map(|(la, lo)| GeoPoint::from_lat_lon(la, lo)),
    );
    // geo coords are (x = lon, y = lat).
    let ga = geo::LineString(path_a.iter().map(|&(la, lo)| geo::Coord { x: lo, y: la }).collect());
    let gb = geo::LineString(path_b.iter().map(|&(la, lo)| geo::Coord { x: lo, y: la }).collect());
    let theirs = Haversine.frechet_distance(&ga, &gb);

    let rel = (ours - theirs).abs() / theirs.max(1.0);
    assert!(rel < 1e-4, "frechet mismatch: ours={ours} geo={theirs} rel={rel}");
}

// ---- helpers ----------------------------------------------------------------

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
