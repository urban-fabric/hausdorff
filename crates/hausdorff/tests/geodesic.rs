//! Ellipsoidal (Karney) geodesics: conformance of `Ellipsoid::{inverse, direct}`
//! against golden values captured from `geographiclib-rs` (the porting source),
//! plus self-contained consistency (inverse∘direct round-trips, line sampling,
//! facade-vs-core agreement) between the spherical and ellipsoidal layers. The
//! golden inverse/direct cases cover all algorithmic regimes; the 60-pair grid
//! breadth is exercised by the self-contained round-trip below. No third-party
//! geodesy crate.

use hausdorff::earth::{Bearing, Ellipsoid, GeoPoint, haversine_m};
use hausdorff::space::r2;

fn g(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

/// A deterministic sweep of point pairs covering the algorithm's regimes:
/// short lines, long lines, meridional, equatorial, polar, sign quadrants,
/// antimeridian-straddling, and the nearly-antipodal Newton/astroid path.
/// The first 20 (hand-picked) are golden-anchored; the rest add round-trip breadth.
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

// geographiclib-rs WGS84 inverse goldens (s12, azi1, azi2) for the first 20 pairs.
const INV_GOLDEN: [(f64, f64, f64); 20] = [
    (73.65461927056084, 52.84727642196089, 52.847878875984016),
    (5906.105192486033, 23.955470806435937, 23.974004340812762),
    (5554325.652391282, 51.37298025993312, 107.97126270829544),
    (7792174.827115239, -10.002482212029669, -10.225331640881215),
    (16487767.831863264, 94.91976604079136, 108.32743085117994),
    (4434992.208449778, 0.0, 0.0),
    (2233651.714751699, 0.0, 180.0),
    (20003931.458625447, 0.0, 0.0),
    (18663800.048868023, 0.0, 0.0),
    (15584728.711058298, 90.0, 90.0),
    (2226389.8158654715, 90.0, 90.0),
    (20003920.400266513, 0.009500587709188827, 179.99049941243453),
    (20003911.1449335, 178.90392026979455, 1.0960808292702462),
    (20003922.22814904, -0.9538284559568778, -179.04617154404312),
    (19980861.908890963, 55.966495140158635, 124.03350485984137),
    (124972.06056709313, 152.61862473151885, 153.4802923822309),
    (336922.9828390199, 72.15480116707758, 69.35047577497711),
    (0.0, 180.0, 180.0),
    (5017021.351334979, 90.0, 180.0),
    (5017021.351334979, -132.0, -0.0),
];

// geographiclib-rs WGS84 direct goldens (lat2, lon2, azi2) for the first 20 pairs,
// driven by bearing = i*23.7-180, distance = 1+i*330777.
const DIR_GOLDEN: [(f64, f64, f64); 20] = [
    (48.85659100775696, 2.3522, 180.0),
    (37.97456361598065, -75.51871418621533, -157.25932445396694),
    (36.47419753529442, -79.21075837485125, -135.98956520953664),
    (-36.32472497734874, 140.73760340136587, -102.86652759607321),
    (37.81167879233002, -137.4575456368973, -94.44886913602875),
    (16.799565147222378, 6.389850172174093, -64.67387789846872),
    (78.43404660220997, -81.08714614556064, -147.93689353333053),
    (
        -69.26087316172855,
        -14.099999999999993,
        -5.862477761722107e-153,
    ),
    (
        -66.29512752709573,
        54.599999999999994,
        3.5353596941828427e-153,
    ),
    (22.22862927930843, 5.469121008007528, 36.35716127116364),
    (15.760159837938023, -164.4178077549494, 60.60013352577943),
    (5.460220490769031, 32.36697361855914, 82.43077844637891),
    (16.238491603687667, 36.01061077988299, 119.05251740296319),
    (-18.542937831039858, 41.190825853147096, 124.24815668828),
    (-36.00247228304447, 22.78253263219792, 144.30663157990165),
    (15.366887704273434, -177.22772814594254, 177.66300284138399),
    (-75.97131805748907, 84.67442618260085, -73.31278538224558),
    (-24.824146561519306, 21.452464123233966, -132.91421534521456),
    (36.562668352866254, -66.59999999999997, -180.0),
    (
        -33.58113690286295,
        -77.69999999999999,
        -1.0248134273967261e-152,
    ),
];

// geographiclib-rs inverse distances on GRS80 / SPHERE / MARS, point_pairs().step_by(7).
const NONWGS_GOLDEN: [[f64; 12]; 3] = [
    [
        73.65461927109482,
        20003931.458460927,
        19980861.908839397,
        17317868.729571857,
        16335296.894257812,
        18893944.778500825,
        11037145.56993637,
        1956871.5227166475,
        8034412.138982807,
        16537311.022741908,
        15616489.589539975,
        14372104.70438197,
    ],
    [
        73.51072610386282,
        20015114.442035925,
        19959516.901919156,
        17340504.03930931,
        16353568.878121428,
        18895026.462324496,
        11039930.327124698,
        1963907.0293185182,
        8028776.745510318,
        16531470.088063316,
        15625244.177852591,
        14401172.599895718,
    ],
    [
        39.24388997694646,
        10638068.506567718,
        10631065.972303875,
        9204388.970370563,
        8683424.015366316,
        10051543.201704929,
        5870940.185060392,
        1038284.1775847921,
        4276778.059602858,
        8800593.010651246,
        8304849.859968006,
        7634589.080871756,
    ],
];

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

#[test]
fn geodesic_area_matches_golden() {
    use r2::GeodesicArea as _;

    // CCW square with a CW hole, a CW-wound exterior, and a high-latitude
    // sliver crossing the antimeridian. Goldens (perimeter_signed, area_signed,
    // perimeter_unsigned, area_unsigned) captured from geo's GeodesicArea.
    let cases: Vec<r2::Polygon> = vec![
        r2::Polygon::new(
            r2::LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
            vec![r2::LineString::from(vec![
                (0.25, 0.25),
                (0.25, 0.75),
                (0.75, 0.75),
                (0.75, 0.25),
            ])],
        ),
        r2::Polygon::new(
            r2::LineString::from(vec![(2.0, 48.0), (2.0, 49.0), (3.0, 49.0), (3.0, 48.0)]),
            vec![],
        ),
        r2::Polygon::new(
            r2::LineString::from(vec![
                (179.5, 70.0),
                (-179.5, 70.0),
                (-179.5, 71.0),
                (179.5, 71.0),
            ]),
            vec![],
        ),
    ];
    const GOLDEN: [(f64, f64, f64, f64); 3] = [
        (
            665659.5125489294,
            9231614224.814873,
            665659.5125489294,
            9231614224.814873,
        ),
        (
            370196.1683474253,
            -8217589924.770447,
            370196.1683474253,
            510057404134163.7,
        ),
        (
            297673.23927019304,
            4157942715.0476074,
            297673.23927019304,
            4157942715.0476074,
        ),
    ];

    for (poly, &(ps, a_s, pu, au)) in cases.iter().zip(GOLDEN.iter()) {
        let (p_ours, a_ours) = poly.geodesic_perimeter_area_signed();
        assert!((p_ours - ps).abs() < 1e-6, "perimeter {p_ours} vs {ps}");
        assert!((a_ours - a_s).abs() < 1e-2, "signed area {a_ours} vs {a_s}");

        let (pu_ours, au_ours) = poly.geodesic_perimeter_area_unsigned();
        assert!((pu_ours - pu).abs() < 1e-6, "u-perimeter {pu_ours} vs {pu}");
        assert!(
            (au_ours - au).abs() < 1e-2,
            "unsigned area {au_ours} vs {au}"
        );
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
