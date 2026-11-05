use super::*;

#[test]
fn web_mercator_matches_proj() {
    let refs: &[Ref] = &[
        (0.0, 0.0, 0.0, 0.0),
        (51.5074, -0.1278, -14226.6309, 6_711_542.4756),
        (-33.8688, 151.2093, 16_832_542.2792, -4_011_198.6473),
        (40.7128, -74.0060, -8_238_310.2356, 4_970_071.5791),
        (85.0, 179.0, 19_926_188.8520, 19_971_868.8804),
    ];
    check("WebMercator", &WebMercator::EPSG_3857, refs, 5e-3);
}

#[test]
fn utm_matches_proj() {
    let utm30n: &[Ref] = &[
        (51.5074, -0.1278, 699_316.2343, 5_710_163.7581),
        (0.0, -3.0, 500_000.0, 0.0),
        (60.0, -1.0, 611_544.0420, 6_653_097.4353),
        (40.0, -6.0, 243_900.3520, 4_432_069.0569),
    ];
    check("UTM30N", &TransverseMercator::utm(30, true), utm30n, 1e-3);
    let utm33n: &[Ref] = &[
        (48.2, 16.37, 601_799.1458, 5_339_437.1538),
        (52.52, 13.405, 391_779.2593, 5_820_072.1592),
        (45.0, 15.0, 500_000.0, 4_982_950.4002),
    ];
    check("UTM33N", &TransverseMercator::utm(33, true), utm33n, 1e-3);
    let utm1s: &[Ref] = &[
        (-40.0, -177.0, 500_000.0, 5_572_242.7813),
        (-10.0, -179.5, 225_928.9466, 8_893_548.7217),
        (-60.0, -176.0, 555_776.2668, 3_348_167.2646),
    ];
    check("UTM1S", &TransverseMercator::utm(1, false), utm1s, 1e-3);
}

/// Transverse Mercator **far from the central meridian** — out to ±20° (≈2000 km
/// easting), where the low-order Snyder series degrades to centimetres. The
/// Karney–Krüger order-n⁶ series stays sub-millimetre against PROJ's accurate
/// `tmerc` across the whole span. CM = 3°E, k0 = 0.9996 (UTM-like).
#[test]
fn transverse_mercator_far_from_central_meridian_matches_proj() {
    let tm = TransverseMercator {
        ellipsoid: Ellipsoid::WGS84,
        central_meridian: 3.0,
        lat_origin: 0.0,
        k0: 0.9996,
        false_easting: 500_000.0,
        false_northing: 0.0,
    };
    // (lat, lon, easting, northing) — PROJ 9.4.0 accurate tmerc, CM = 3°E.
    let far: &[Ref] = &[
        (45.0, 3.0, 500_000.000000, 4_982_950.400227),
        (45.0, 13.0, 1_288_141.060230, 5_031_833.622250), // +10°
        (45.0, 23.0, 2_075_750.235110, 5_181_473.169330), // +20°
        (45.0, -7.0, -288_141.060230, 5_031_833.622250),  // −10°
        (60.0, 18.0, 1_331_808.333412, 6_746_522.323960), // +15°, high lat
        (20.0, -12.0, -1_083_030.922104, 2_283_512.921526), // −15°
        (-35.0, 25.0, 2_523_906.171801, -4_102_603.370651), // +22°, southern
    ];
    check("TM far-CM", &tm, far, 1e-4);
}

/// **3-D geographic height vs PROJ.** WGS84 geographic-3D (EPSG:4979) → geocentric
/// ECEF (EPSG:4978): the ellipsoidal height rides through the transform pivot.
/// Validates the 3-D path against PROJ 9.4.0 (not just the internal round-trip).
#[test]
fn geographic_3d_height_matches_proj() {
    let g3d = Crs::from_epsg(4979).unwrap();
    let ecef = Crs::from_epsg(4978).unwrap();
    // (lat, lon, height_m, X, Y, Z) — pyproj 4979→4978.
    let refs: &[(f64, f64, f64, f64, f64, f64)] = &[
        (
            45.0,
            45.0,
            1000.0,
            3_194_919.1451,
            3_194_919.1451,
            4_488_055.5156,
        ),
        (
            51.5,
            -0.1,
            100.0,
            3_978_704.7222,
            -6944.1568,
            4_968_440.7181,
        ),
        (
            -33.87,
            151.21,
            50.0,
            -4_646_053.8387,
            2_553_133.8488,
            -3_534_510.7711,
        ),
        (0.0, 0.0, 0.0, 6_378_137.0, 0.0, 0.0),
        (
            80.0,
            -120.0,
            8848.0,
            -556_350.6549,
            -963_627.6012,
            6_268_256.5400,
        ),
    ];
    let mut max_m = 0.0_f64;
    for &(lat, lon, h, x, y, z) in refs {
        let p = CrsCoord::geographic_with_height(GeoPoint::from_lat_lon(lat, lon), h);
        let CrsCoord::Geocentric(got) = transform(g3d, ecef, p) else {
            panic!("expected geocentric")
        };
        max_m = max_m
            .max((got.x - x).abs())
            .max((got.y - y).abs())
            .max((got.z - z).abs());
        // Inverse of our own forward recovers lat/lon AND ellipsoidal height
        // exactly (PROJ's reference is only transcribed to 0.1 mm).
        let CrsCoord::Geographic {
            point: back,
            height: hb,
        } = transform(ecef, g3d, CrsCoord::Geocentric(got))
        else {
            panic!("expected geographic")
        };
        assert!(
            (back.lat - lat).abs() < 1e-9 && (hb - h).abs() < 1e-6,
            "inverse h {h} → {hb}"
        );
    }
    eprintln!("4979→4978 (3-D height) max Δ = {max_m:.6} m vs PROJ");
    assert!(
        max_m < 1e-3,
        "3-D height transform off by {max_m} m vs PROJ"
    );
}
