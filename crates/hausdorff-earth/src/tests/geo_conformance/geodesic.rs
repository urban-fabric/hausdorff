use super::simplify::h_line;
use super::*;

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
        assert!(
            (h - want).abs() <= 1e-9 * want.max(1.0),
            "CD area on {f:?}: {h} vs {want}"
        );
    }

    let o2 = [
        (0.00388383, 51.501574),
        (0.00538587, 51.502278),
        (0.00553607, 51.503299),
        (0.00467777, 51.504181),
        (0.00327229, 51.504435),
        (0.00187754, 51.504168),
        (0.00087976, 51.503380),
        (0.00107288, 51.502324),
        (0.00185608, 51.501770),
    ];
    assert_eq!(
        h_poly(&o2).chamberlain_duquette_unsigned_area().round(),
        78_478.0
    );
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
        let (a, b) = (
            hausdorff_planar::space::r2::Point::new(x1, y1),
            hausdorff_planar::space::r2::Point::new(x2, y2),
        );
        let d = hausdorff_planar::space::r2::Haversine.distance(a, b);
        assert!((d - wd).abs() <= 1e-9 * wd, "distance {d} vs {wd}");
        let br = hausdorff_planar::space::r2::Haversine.bearing(a, b);
        assert!(
            ((br - wb + 540.0) % 360.0 - 180.0).abs() < 1e-9,
            "bearing {br} vs {wb}"
        );
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
    assert!(
        (hausdorff_planar::space::r2::Euclidean.length(&h_line(&coords)) - euclid_sum).abs() < EPS,
        "euclidean length"
    );
    // Haversine length: golden.
    assert!(
        (hausdorff_planar::space::r2::Haversine.length(&h_line(&coords)) - 730594.5476435976).abs()
            <= 1e-6,
        "haversine length"
    );
}

#[test]
fn haversine_destination_matches_golden() {
    let origin = hausdorff_planar::space::r2::Point::new(13.405, 52.52);
    let courses = [
        (0.0, 1000.0),
        (90.0, 5000.0),
        (215.0, 250.0),
        (359.0, 12000.0),
    ];
    const GOLDEN: [(f64, f64); 4] = [
        (13.404999999999973, 52.52899320363725),
        (13.478898358112474, 52.519976988255365),
        (13.402880770254797, 52.51815828078832),
        (13.401897084247594, 52.62790196665409),
    ];
    for (&(bearing, dist), &(gx, gy)) in courses.iter().zip(GOLDEN.iter()) {
        let d = hausdorff_planar::space::r2::Haversine.destination(origin, bearing, dist);
        assert!(
            (d.x() - gx).abs() < 1e-9 && (d.y() - gy).abs() < 1e-9,
            "destination {bearing}/{dist}"
        );
    }
}

#[test]
fn haversine_interpolation_follows_the_great_circle() {
    let start = hausdorff_planar::space::r2::Point::new(-122.4194, 37.7749);
    let end = hausdorff_planar::space::r2::Point::new(-73.9857, 40.7484);
    let total = hausdorff_planar::space::r2::Haversine.distance(start, end);

    let half = hausdorff_planar::space::r2::Haversine.point_at_ratio_between(start, end, 0.5);
    let first_half = hausdorff_planar::space::r2::Haversine.distance(start, half);
    let second_half = hausdorff_planar::space::r2::Haversine.distance(half, end);
    assert!((first_half - total * 0.5).abs() <= total * 1e-9);
    assert!((second_half - total * 0.5).abs() <= total * 1e-9);

    let waypoint =
        hausdorff_planar::space::r2::Haversine.point_at_distance_between(start, end, 100_000.0);
    let travelled = hausdorff_planar::space::r2::Haversine.distance(start, waypoint);
    assert!((travelled - 100_000.0).abs() <= 1e-4);
}
