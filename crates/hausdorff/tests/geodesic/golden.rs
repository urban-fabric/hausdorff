use super::*;

#[test]
fn inverse_matches_golden() {
    for ((a, b), &(s12, azi1, azi2)) in point_pairs().into_iter().zip(INV_GOLDEN.iter()) {
        let ours = Ellipsoid::WGS84.inverse(a, b);
        assert!(
            (ours.distance_m - s12).abs() < 1e-8,
            "distance {} vs {s12} for {a:?} → {b:?}",
            ours.distance_m
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
fn direct_matches_golden() {
    for (i, ((origin, _), &(lat2, lon2, azi2))) in
        point_pairs().into_iter().zip(DIR_GOLDEN.iter()).enumerate()
    {
        let bearing = (i as f64) * 23.7 - 180.0;
        let distance = 1.0 + (i as f64) * 330_777.0; // 1 m … ~6,000 km
        let ours = Ellipsoid::WGS84.direct(origin, Bearing::from_degrees(bearing), distance);
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
fn non_wgs84_ellipsoids_match_golden() {
    // The uncached path: GRS80, a sphere, and Mars' much larger flattening.
    for (e, golden) in [Ellipsoid::GRS80, Ellipsoid::SPHERE, Ellipsoid::MARS]
        .into_iter()
        .zip(NONWGS_GOLDEN.iter())
    {
        for ((a, b), &s12) in point_pairs().into_iter().step_by(7).zip(golden.iter()) {
            let ours = e.inverse(a, b);
            assert!(
                (ours.distance_m - s12).abs() < 1e-8,
                "distance {} vs {s12} on {e:?} for {a:?} → {b:?}",
                ours.distance_m
            );
        }
    }
}
