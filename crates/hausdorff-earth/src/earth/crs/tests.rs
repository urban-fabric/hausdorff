use super::*;

fn approx(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

#[test]
fn epsg_registry_recognises_the_core_systems() {
    assert_eq!(Crs::from_epsg(4326), Some(Crs::wgs84()));
    assert!(matches!(Crs::from_epsg(4978), Some(Crs::Geocentric { .. })));
    assert!(matches!(Crs::from_epsg(3857), Some(Crs::Projected { .. })));
    // UTM zone 30N / 30S.
    assert_eq!(Crs::from_epsg(32_630).unwrap().epsg(), Some(32_630));
    assert!(matches!(Crs::from_epsg(32_730), Some(Crs::Projected { .. })));
    assert_eq!(Crs::from_epsg(9999), None);
}

#[test]
fn geographic_to_webmercator_and_back() {
    let wgs84 = Crs::from_epsg(4326).unwrap();
    let webmerc = Crs::from_epsg(3857).unwrap();
    // Null Island maps to the planar origin.
    let origin = transform(
        wgs84,
        webmerc,
        CrsCoord::Geographic(GeoPoint::from_lat_lon(0.0, 0.0)),
    );
    assert!(matches!(origin, CrsCoord::Projected(c) if approx(c.x, 0.0, 1e-6) && approx(c.y, 0.0, 1e-6)));
    // London (51.5074, -0.1278): EPSG:3857 ≈ (-14_226, 6_711_544) m.
    let london = GeoPoint::from_lat_lon(51.5074, -0.1278);
    let proj = transform(wgs84, webmerc, CrsCoord::Geographic(london)).projected();
    assert!(approx(proj.x, -14_226.0, 5.0), "x = {}", proj.x);
    assert!(approx(proj.y, 6_711_544.0, 50.0), "y = {}", proj.y);
    // Round-trips back to the same lat/lon.
    let back = transform(webmerc, wgs84, CrsCoord::Projected(proj)).geographic();
    assert!(approx(back.lat, london.lat, 1e-9) && approx(back.lon, london.lon, 1e-9));
}

#[test]
fn geographic_to_utm_and_back() {
    let wgs84 = Crs::from_epsg(4326).unwrap();
    let utm30n = Crs::from_epsg(32_630).unwrap();
    // London is in UTM zone 30N; easting near 500 km (central meridian −3°),
    // northing ≈ 5.7 Mm.
    let london = GeoPoint::from_lat_lon(51.5074, -0.1278);
    let utm = transform(wgs84, utm30n, CrsCoord::Geographic(london)).projected();
    assert!(approx(utm.x, 699_330.0, 200.0), "easting = {}", utm.x);
    assert!(approx(utm.y, 5_710_158.0, 200.0), "northing = {}", utm.y);
    // Round-trip.
    let back = transform(utm30n, wgs84, CrsCoord::Projected(utm)).geographic();
    assert!(approx(back.lat, london.lat, 1e-7) && approx(back.lon, london.lon, 1e-7));
}

#[test]
fn geographic_to_geocentric_ecef_and_back() {
    let wgs84 = Crs::from_epsg(4326).unwrap();
    let ecef = Crs::from_epsg(4978).unwrap();
    let p = GeoPoint::from_lat_lon(45.0, 45.0);
    let xyz = transform(wgs84, ecef, CrsCoord::Geographic(p)).geocentric();
    // On the WGS84 ellipsoid at 45/45: a known ECEF (~3.19e6, 3.19e6, 4.49e6).
    assert!(approx(xyz.x, 3_194_419.0, 50.0), "x = {}", xyz.x);
    assert!(approx(xyz.z, 4_487_348.0, 50.0), "z = {}", xyz.z);
    let back = transform(ecef, wgs84, CrsCoord::Geocentric(xyz)).geographic();
    assert!(approx(back.lat, 45.0, 1e-9) && approx(back.lon, 45.0, 1e-9));
}

#[test]
fn projected_to_projected_pivots_through_ecef() {
    // UTM 30N → Web Mercator directly (projected→projected).
    let utm30n = Crs::from_epsg(32_630).unwrap();
    let webmerc = Crs::from_epsg(3857).unwrap();
    let london = GeoPoint::from_lat_lon(51.5074, -0.1278);
    let via_utm = transform(Crs::wgs84(), utm30n, CrsCoord::Geographic(london));
    let webmerc_coord = transform(utm30n, webmerc, via_utm).projected();
    let direct = transform(Crs::wgs84(), webmerc, CrsCoord::Geographic(london)).projected();
    assert!(approx(webmerc_coord.x, direct.x, 1e-3));
    assert!(approx(webmerc_coord.y, direct.y, 1e-3));
}

#[test]
fn phase2_projected_crses_round_trip_through_the_registry() {
    let wgs84 = Crs::wgs84();
    // Each new projected CRS: WGS84 geographic → CRS → back is the identity.
    for (code, lat, lon) in [
        (5070_u32, 40.0, -100.0), // Albers CONUS
        (2154, 46.5, 2.5),        // Lambert-93 (France)
        (3995, 84.0, 40.0),       // Arctic Polar Stereographic
        (3031, -82.0, -60.0),     // Antarctic Polar Stereographic
    ] {
        let crs = Crs::from_epsg(code).unwrap();
        let p = GeoPoint::from_lat_lon(lat, lon);
        let proj = transform(wgs84, crs, CrsCoord::Geographic(p));
        assert!(matches!(proj, CrsCoord::Projected(_)), "EPSG:{code} is projected");
        let back = transform(crs, wgs84, proj).geographic();
        assert!(
            approx(back.lat, lat, 1e-6) && approx(back.lon, lon, 1e-6),
            "EPSG:{code} round-trip: back=({},{}) dlat={:e} dlon={:e}",
            back.lat,
            back.lon,
            back.lat - lat,
            back.lon - lon
        );
    }
}

#[test]
fn british_national_grid_matches_os_worked_example() {
    // Within the OSGB36 datum (no shift), EPSG:4277 → EPSG:27700 is the pure
    // BNG projection. OS "Guide to Coordinate Systems" worked example
    // (Caister water tower): OSGB36 (52.6575703°, 1.7179216°) → BNG
    // (651409.903, 313177.270) m.
    let osgb_geo = Crs::from_epsg(4277).unwrap();
    let bng = Crs::from_epsg(27_700).unwrap();
    let p = GeoPoint::from_lat_lon(52.657_570_3, 1.717_921_6);
    let en = transform(osgb_geo, bng, CrsCoord::Geographic(p)).projected();
    assert!(approx(en.x, 651_409.903, 0.5), "easting = {}", en.x);
    assert!(approx(en.y, 313_177.270, 0.5), "northing = {}", en.y);
    // Round-trips within the datum.
    let back = transform(bng, osgb_geo, CrsCoord::Projected(en)).geographic();
    assert!(approx(back.lat, p.lat, 1e-7) && approx(back.lon, p.lon, 1e-7));
}

#[test]
fn cross_datum_shift_is_nonzero_and_reversible() {
    let wgs84 = Crs::wgs84();
    // The same physical location has a *different* lat/lon on OSGB36 — the
    // datum shift is real (order tens of metres), not zero.
    for crs in [
        Crs::from_epsg(4277).unwrap(), // OSGB36
        Crs::from_epsg(4230).unwrap(), // ED50
        Crs::from_epsg(4267).unwrap(), // NAD27
    ] {
        let p = GeoPoint::from_lat_lon(51.5, -0.1);
        let shifted = transform(wgs84, crs, CrsCoord::Geographic(p)).geographic();
        // Non-trivial shift (more than ~1 m ≈ 1e-5°, less than ~1 km).
        let dlat = (shifted.lat - p.lat).abs();
        let dlon = (shifted.lon - p.lon).abs();
        assert!(dlat > 1e-5 || dlon > 1e-5, "{:?}: shift too small", crs.epsg());
        assert!(dlat < 1e-2 && dlon < 1e-2, "{:?}: shift too large", crs.epsg());
        // And the shift reverses exactly (helmert_round_trip).
        let back = transform(crs, wgs84, CrsCoord::Geographic(shifted)).geographic();
        assert!(
            approx(back.lat, p.lat, 1e-7) && approx(back.lon, p.lon, 1e-7),
            "{:?}: round-trip",
            crs.epsg()
        );
    }
}

#[test]
fn ntv2_grid_datum_matches_proj_end_to_end() {
    // Grid-shifted geographic CRS → WGS84 through the full registry pipeline,
    // which uses the embedded NTv2 grids. PROJ 9.4.0 ground truth.
    let wgs84 = Crs::wgs84();
    let cases: &[(u32, &str, &[(f64, f64, f64, f64)])] = &[
        (
            4314, // DHDN — BETA2007
            "DHDN",
            &[
                (52.0, 10.0, 51.998_628_928, 9.998_786_383),
                (48.5, 9.2, 48.499_021_746, 9.198_951_791),
                (49.0, 7.0, 48.998_951_376, 6.999_247_781),
            ],
        ),
        (
            4275, // NTF — ntf_r93
            "NTF",
            &[
                (48.85, 2.35, 48.849_933_563, 2.349_295_594),
                (43.6, 1.44, 43.599_980_695, 1.439_310_962),
                (45.75, 4.85, 45.749_971_904, 4.849_417_525),
            ],
        ),
        (
            4272, // NZGD49 — nzgd2kgrid0005
            "NZGD49",
            &[
                (-36.85, 174.76, -36.848_196_691, 174.760_191_647),
                (-41.29, 174.78, -41.288_275_516, 174.780_190_614),
                (-43.53, 172.64, -43.528_327_299, 172.640_130_644),
            ],
        ),
    ];
    for &(code, name, refs) in cases {
        let crs = Crs::from_epsg(code).unwrap();
        let mut max_deg = 0.0_f64;
        for &(lat, lon, rlat, rlon) in refs {
            let p = CrsCoord::Geographic(GeoPoint::from_lat_lon(lat, lon));
            let CrsCoord::Geographic(g) = transform(crs, wgs84, p) else { panic!() };
            max_deg = max_deg.max((g.lat - rlat).abs()).max((g.lon - rlon).abs());
            // Round-trips: source → WGS84 → source → WGS84.
            let back = transform(wgs84, crs, transform(crs, wgs84, p));
            let CrsCoord::Geographic(b) = transform(crs, wgs84, back) else { panic!() };
            assert!((b.lat - rlat).abs() < 1e-7 && (b.lon - rlon).abs() < 1e-7);
        }
        eprintln!("{name}→WGS84 (grid) end-to-end max Δ = {max_deg:.3e}° vs PROJ");
        assert!(max_deg < 1e-7, "{name} grid datum pipeline off by {max_deg}°");
    }
}

#[test]
fn local_tangent_plane_crs_gives_site_local_metres() {
    // A local CRS anchored at a site origin; WGS84 → local East/North metres.
    let origin = GeoPoint::from_lat_lon(52.0, 10.0);
    let site = Crs::local_tangent_plane(origin);
    let wgs84 = Crs::wgs84();
    // The origin is the local (0, 0).
    let o = transform(wgs84, site, CrsCoord::Geographic(origin)).projected();
    assert!(o.x.abs() < 1e-6 && o.y.abs() < 1e-6);
    // PROJ topocentric: (52.05, 10.07) → E≈4802.10, N≈5565.70 m.
    let p = GeoPoint::from_lat_lon(52.05, 10.07);
    let local = transform(wgs84, site, CrsCoord::Geographic(p)).projected();
    assert!((local.x - 4802.102).abs() < 1e-2, "E={}", local.x);
    assert!((local.y - 5565.702).abs() < 1e-2, "N={}", local.y);
    // Round-trips back to the same lat/lon.
    let back = transform(site, wgs84, CrsCoord::Projected(local)).geographic();
    assert!((back.lat - p.lat).abs() < 1e-9 && (back.lon - p.lon).abs() < 1e-9);
}

#[test]
fn helmert_apply_inverse_round_trips() {
    // A non-trivial Bursa–Wolf (OSGB36→WGS84 params) round-trips exactly.
    let h = Helmert::bursa_wolf(446.448, -125.157, 542.06, 0.15, 0.247, 0.842, -20.489);
    assert!(!h.is_identity());
    let p = EcefPoint {
        x: 3_980_000.0,
        y: -10_000.0,
        z: 4_970_000.0,
    };
    let back = h.inverse_apply(h.apply(p));
    assert!(approx(back.x, p.x, 1e-6), "x {} vs {}", back.x, p.x);
    assert!(approx(back.y, p.y, 1e-6), "y {} vs {}", back.y, p.y);
    assert!(approx(back.z, p.z, 1e-6), "z {} vs {}", back.z, p.z);
    // The identity is a no-op.
    assert_eq!(Helmert::IDENTITY.apply(p), p);
}
