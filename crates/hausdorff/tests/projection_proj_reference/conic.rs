use super::*;

#[test]
fn albers_conus_matches_proj() {
    let refs: &[Ref] = &[
        (40.0, -100.0, -338_390.5876, 1_894_100.1400),
        (30.0, -80.0, 1_534_849.0388, 898_886.0884),
        (48.0, -120.0, -1_786_745.6096, 3_004_984.2834),
        (25.0, -97.0, -102_181.5305, 218_657.5760),
    ];
    let albers = AlbersEqualArea {
        ellipsoid: Ellipsoid::GRS80,
        lat1: 29.5,
        lat2: 45.5,
        lat0: 23.0,
        lon0: -96.0,
        false_easting: 0.0,
        false_northing: 0.0,
    };
    check("Albers5070", &albers, refs, 5e-3);
}

#[test]
fn lambert93_matches_proj() {
    let refs: &[Ref] = &[
        (46.5, 3.0, 700_000.0, 6_600_000.0),
        (48.85, 2.35, 652_301.5648, 6_861_302.7259),
        (43.3, 5.4, 894_829.4091, 6_247_498.2043),
        (50.6, 3.05, 703_545.8177, 7_055_882.5829),
    ];
    let lcc = LambertConformalConic {
        ellipsoid: Ellipsoid::GRS80,
        lat1: 44.0,
        lat2: 49.0,
        lat0: 46.5,
        lon0: 3.0,
        false_easting: 700_000.0,
        false_northing: 6_600_000.0,
    };
    check("Lambert93", &lcc, refs, 5e-3);
}

#[test]
fn polar_stereographic_matches_proj() {
    let arctic: &[Ref] = &[
        (84.0, 40.0, 419_410.5411, -499_834.0191),
        (72.0, -120.0, -1_707_230.1002, 985_669.7579),
        (89.0, 175.0, 9_469.9150, 108_241.6235),
        (75.0, 0.0, 0.0, -1_638_783.2384),
    ];
    check(
        "Arctic3995",
        &PolarStereographic {
            ellipsoid: Ellipsoid::WGS84,
            lat_ts: 71.0,
            lon0: 0.0,
            north: true,
            false_easting: 0.0,
            false_northing: 0.0,
        },
        arctic,
        5e-3,
    );
    let antarctic: &[Ref] = &[
        (-84.0, 40.0, 419_410.5411, 499_834.0191),
        (-72.0, -120.0, -1_707_230.1002, -985_669.7579),
        (-89.0, 175.0, 9_469.9150, -108_241.6235),
        (-75.0, 0.0, 0.0, 1_638_783.2384),
    ];
    check(
        "Antarc3031",
        &PolarStereographic {
            ellipsoid: Ellipsoid::WGS84,
            lat_ts: -71.0,
            lon0: 0.0,
            north: false,
            false_easting: 0.0,
            false_northing: 0.0,
        },
        antarctic,
        5e-3,
    );
}

#[test]
fn british_national_grid_matches_proj() {
    let refs: &[Ref] = &[
        (52.6575703, 1.7179216, 651_409.9041, 313_177.2698),
        (51.5, -0.1, 531_866.1304, 179_660.9048),
        (55.0, -3.0, 336_038.7622, 567_734.8572),
        (50.1, -5.5, 149_737.3470, 28_150.4517),
    ];
    let bng = TransverseMercator {
        ellipsoid: airy1830(),
        central_meridian: -2.0,
        lat_origin: 49.0,
        k0: 0.999_601_271_7,
        false_easting: 400_000.0,
        false_northing: -100_000.0,
    };
    check("BNG27700", &bng, refs, 5e-3);
}

/// The **full datum-shift pipeline** (Helmert + ellipsoid + projection) for the
/// safety-critical "GPS WGS84 → national datum/grid" operation, vs PROJ's
/// 7-parameter transform. Each reference is PROJ transforming `4326 → target`.
#[test]
fn datum_shift_pipeline_matches_proj() {
    // WGS84 (4326) → OSGB36 geographic (4277): lat/lon after the datum shift.
    // The shift is ~tens of metres — getting it wrong mislocates by that much.
    let to_osgb36: &[(f64, f64, f64, f64)] = &[
        (52.6575703, 1.7179216, 52.657_161_919, 1.719_791_404),
        (51.5, -0.1, 51.499_488_290, -0.098_390_664),
        (55.0, -3.0, 54.999_941_910, -2.998_591_795),
        (50.1, -5.5, 50.099_397_697, -5.499_049_126),
    ];
    let osgb = Crs::from_epsg(4277).unwrap();
    let mut max_deg = 0.0_f64;
    for &(lat, lon, rlat, rlon) in to_osgb36 {
        let out = transform(
            Crs::wgs84(),
            osgb,
            CrsCoord::geographic(GeoPoint::from_lat_lon(lat, lon)),
        );
        let CrsCoord::Geographic { point: g, .. } = out else {
            panic!("expected geographic")
        };
        max_deg = max_deg.max((g.lat - rlat).abs()).max((g.lon - rlon).abs());
    }
    eprintln!("WGS84→OSGB36 (4277) max Δ = {max_deg:.3e}° vs PROJ Helmert");
    assert!(max_deg < 1e-7, "OSGB36 datum shift off by {max_deg}°");

    // WGS84 (4326) → British National Grid (27700): the full GPS→grid pipeline.
    // Note the easting differs from the same-datum projection by ~129 m — the
    // datum shift — so this exercises Helmert + Airy ellipsoid + TM together.
    let to_bng: &[(f64, f64, f64, f64)] = &[
        (52.6575703, 1.7179216, 651_538.6418, 313_138.4039),
        (51.5, -0.1, 531_979.2927, 179_606.9085),
        (55.0, -3.0, 336_128.7373, 567_727.1070),
        (50.1, -5.5, 149_802.1828, 28_080.3225),
    ];
    let bng = Crs::from_epsg(27_700).unwrap();
    let mut max_m = 0.0_f64;
    for &(lat, lon, re, rn) in to_bng {
        let out = transform(
            Crs::wgs84(),
            bng,
            CrsCoord::geographic(GeoPoint::from_lat_lon(lat, lon)),
        );
        let CrsCoord::Projected(c) = out else {
            panic!("expected projected")
        };
        max_m = max_m.max((c.x - re).abs()).max((c.y - rn).abs());
    }
    eprintln!("WGS84→BNG (27700) full pipeline max Δ = {max_m:.6} m vs PROJ Helmert");
    assert!(max_m < 5e-3, "BNG datum pipeline off by {max_m} m");
}
