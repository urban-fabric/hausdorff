use super::simplify::is_near_boundary;
use super::*;

#[test]
fn area_matches_shoelace() {
    for f in fixtures() {
        let hp = h_poly(&f);
        let signed = shoelace(&f);
        assert!(
            (hp.signed_area() - signed).abs() < EPS,
            "signed_area on {f:?}"
        );
        assert!(
            (hp.unsigned_area() - signed.abs()).abs() < EPS,
            "unsigned_area on {f:?}"
        );
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
        assert!(
            (c.x() - gx).abs() < EPS && (c.y() - gy).abs() < EPS,
            "centroid on {f:?}"
        );
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
                let hc = hp.contains(&hausdorff_planar::space::r2::Coord::new(x, y));
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
    let point_sets = [
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (2.0, 2.0)],
        vec![
            (0.0, 0.0),
            (1.0, 0.5),
            (2.0, 0.0),
            (2.0, 2.0),
            (1.0, 1.5),
            (0.0, 2.0),
            (1.0, 1.0),
        ],
        vec![
            (-3.0, -1.0),
            (5.0, 0.0),
            (3.0, 6.0),
            (-2.0, 4.0),
            (0.0, 0.0),
            (1.0, 2.0),
        ],
    ];
    const GOLDEN: [f64; 3] = [16.0, 4.0, 36.5];
    for (pts, &want) in point_sets.iter().zip(GOLDEN.iter()) {
        let a = hausdorff_planar::space::r2::MultiPoint::from(pts.clone())
            .convex_hull()
            .unsigned_area();
        assert!(
            (a - want).abs() < 1e-9,
            "convex hull area on {pts:?}: {a} vs {want}"
        );
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
    const GOLDEN: [f64; 4] = [
        877464.5379215094,
        3941829.06891843,
        111195.0802335329,
        7825755.604165367,
    ];
    for (((lat1, lon1), (lat2, lon2)), &want) in samples.iter().zip(GOLDEN.iter()) {
        let ours = haversine_m(
            GeoPoint::from_lat_lon(*lat1, *lon1),
            GeoPoint::from_lat_lon(*lat2, *lon2),
        );
        assert!(
            (ours - want).abs() <= 1e-9 * want,
            "haversine: {ours} vs {want}"
        );
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
