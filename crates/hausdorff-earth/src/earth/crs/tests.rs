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
        CrsCoord::geographic(GeoPoint::from_lat_lon(0.0, 0.0)),
    );
    assert!(matches!(origin, CrsCoord::Projected(c) if approx(c.x, 0.0, 1e-6) && approx(c.y, 0.0, 1e-6)));
    // London (51.5074, -0.1278): EPSG:3857 ≈ (-14_226, 6_711_544) m.
    let london = GeoPoint::from_lat_lon(51.5074, -0.1278);
    let proj = transform(wgs84, webmerc, CrsCoord::geographic(london)).projected();
    assert!(approx(proj.x, -14_226.0, 5.0), "x = {}", proj.x);
    assert!(approx(proj.y, 6_711_544.0, 50.0), "y = {}", proj.y);
    // Round-trips back to the same lat/lon.
    let back = transform(webmerc, wgs84, CrsCoord::Projected(proj)).point();
    assert!(approx(back.lat, london.lat, 1e-9) && approx(back.lon, london.lon, 1e-9));
}

#[test]
fn geographic_to_utm_and_back() {
    let wgs84 = Crs::from_epsg(4326).unwrap();
    let utm30n = Crs::from_epsg(32_630).unwrap();
    // London is in UTM zone 30N; easting near 500 km (central meridian −3°),
    // northing ≈ 5.7 Mm.
    let london = GeoPoint::from_lat_lon(51.5074, -0.1278);
    let utm = transform(wgs84, utm30n, CrsCoord::geographic(london)).projected();
    assert!(approx(utm.x, 699_330.0, 200.0), "easting = {}", utm.x);
    assert!(approx(utm.y, 5_710_158.0, 200.0), "northing = {}", utm.y);
    // Round-trip.
    let back = transform(utm30n, wgs84, CrsCoord::Projected(utm)).point();
    assert!(approx(back.lat, london.lat, 1e-7) && approx(back.lon, london.lon, 1e-7));
}

#[test]
fn geographic_to_geocentric_ecef_and_back() {
    let wgs84 = Crs::from_epsg(4326).unwrap();
    let ecef = Crs::from_epsg(4978).unwrap();
    let p = GeoPoint::from_lat_lon(45.0, 45.0);
    let xyz = transform(wgs84, ecef, CrsCoord::geographic(p)).geocentric();
    // On the WGS84 ellipsoid at 45/45: a known ECEF (~3.19e6, 3.19e6, 4.49e6).
    assert!(approx(xyz.x, 3_194_419.0, 50.0), "x = {}", xyz.x);
    assert!(approx(xyz.z, 4_487_348.0, 50.0), "z = {}", xyz.z);
    let back = transform(ecef, wgs84, CrsCoord::Geocentric(xyz)).point();
    assert!(approx(back.lat, 45.0, 1e-9) && approx(back.lon, 45.0, 1e-9));
}

#[test]
fn projected_to_projected_pivots_through_ecef() {
    // UTM 30N → Web Mercator directly (projected→projected).
    let utm30n = Crs::from_epsg(32_630).unwrap();
    let webmerc = Crs::from_epsg(3857).unwrap();
    let london = GeoPoint::from_lat_lon(51.5074, -0.1278);
    let via_utm = transform(Crs::wgs84(), utm30n, CrsCoord::geographic(london));
    let webmerc_coord = transform(utm30n, webmerc, via_utm).projected();
    let direct = transform(Crs::wgs84(), webmerc, CrsCoord::geographic(london)).projected();
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
        let proj = transform(wgs84, crs, CrsCoord::geographic(p));
        assert!(matches!(proj, CrsCoord::Projected(_)), "EPSG:{code} is projected");
        let back = transform(crs, wgs84, proj).point();
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
    let en = transform(osgb_geo, bng, CrsCoord::geographic(p)).projected();
    assert!(approx(en.x, 651_409.903, 0.5), "easting = {}", en.x);
    assert!(approx(en.y, 313_177.270, 0.5), "northing = {}", en.y);
    // Round-trips within the datum.
    let back = transform(bng, osgb_geo, CrsCoord::Projected(en)).point();
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
        let shifted = transform(wgs84, crs, CrsCoord::geographic(p)).point();
        // Non-trivial shift (more than ~1 m ≈ 1e-5°, less than ~1 km).
        let dlat = (shifted.lat - p.lat).abs();
        let dlon = (shifted.lon - p.lon).abs();
        assert!(dlat > 1e-5 || dlon > 1e-5, "{:?}: shift too small", crs.epsg());
        assert!(dlat < 1e-2 && dlon < 1e-2, "{:?}: shift too large", crs.epsg());
        // And the shift reverses exactly (helmert_round_trip).
        let back = transform(crs, wgs84, CrsCoord::geographic(shifted)).point();
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
            let p = CrsCoord::geographic(GeoPoint::from_lat_lon(lat, lon));
            let CrsCoord::Geographic { point: g, .. } = transform(crs, wgs84, p) else { panic!() };
            max_deg = max_deg.max((g.lat - rlat).abs()).max((g.lon - rlon).abs());
            // Round-trips: source → WGS84 → source → WGS84.
            let back = transform(wgs84, crs, transform(crs, wgs84, p));
            let CrsCoord::Geographic { point: b, .. } = transform(crs, wgs84, back) else { panic!() };
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
    let o = transform(wgs84, site, CrsCoord::geographic(origin)).projected();
    assert!(o.x.abs() < 1e-6 && o.y.abs() < 1e-6);
    // PROJ topocentric: (52.05, 10.07) → E≈4802.10, N≈5565.70 m.
    let p = GeoPoint::from_lat_lon(52.05, 10.07);
    let local = transform(wgs84, site, CrsCoord::geographic(p)).projected();
    assert!((local.x - 4802.102).abs() < 1e-2, "E={}", local.x);
    assert!((local.y - 5565.702).abs() < 1e-2, "N={}", local.y);
    // Round-trips back to the same lat/lon.
    let back = transform(site, wgs84, CrsCoord::Projected(local)).point();
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

#[test]
fn broadened_registry_recognises_the_new_codes() {
    // Geographic bases on GRS80 (NAD83 / ETRS89).
    for code in [4269_u32, 4258] {
        assert!(matches!(Crs::from_epsg(code), Some(Crs::Geographic { .. })));
        assert_eq!(Crs::from_epsg(code).unwrap().epsg(), Some(code));
        assert_eq!(Crs::from_epsg(code).unwrap().datum().ellipsoid, Ellipsoid::GRS80);
    }
    // A projected representative from every new family, by EPSG code.
    for code in [
        25_828_u32, 25_838, 26_901, 26_923, 28_348, 28_358, 2193, 2180, 3978, 3577, 3413, 3976,
    ] {
        let crs = Crs::from_epsg(code).unwrap_or_else(|| panic!("EPSG:{code} unrecognised"));
        assert!(matches!(crs, Crs::Projected { .. }), "EPSG:{code} should be projected");
        assert_eq!(crs.epsg(), Some(code));
    }
    // Zone ranges are bounded — one past each end is not recognised.
    assert_eq!(Crs::from_epsg(25_827), None);
    assert_eq!(Crs::from_epsg(25_839), None);
    assert_eq!(Crs::from_epsg(28_347), None);
    assert_eq!(Crs::from_epsg(28_359), None);
}

#[test]
fn new_projected_crses_round_trip_through_the_registry() {
    let wgs84 = Crs::wgs84();
    // (code, an in-area lat/lon): WGS84 geographic → CRS → back is the identity.
    for (code, lat, lon) in [
        (25_828_u32, 28.1, -15.4),  // ETRS89 / UTM 28N (Canaries)
        (26_923, 45.0, -45.0),      // NAD83 / UTM 23N
        (28_356, -33.87, 151.21),   // GDA94 / MGA 56 (Sydney)
        (2193, -41.29, 174.78),     // NZTM2000 (Wellington)
        (2180, 52.23, 21.01),       // ETRS89 / Poland CS92 (Warsaw)
        (3978, 60.0, -95.0),        // Canada Atlas Lambert
        (3577, -25.0, 133.0),       // Australian Albers
        (3413, 75.0, -40.0),        // NSIDC Polar Stereographic North
        (3976, -75.0, 0.0),         // NSIDC Polar Stereographic South
    ] {
        let crs = Crs::from_epsg(code).unwrap();
        let p = GeoPoint::from_lat_lon(lat, lon);
        let proj = transform(wgs84, crs, CrsCoord::geographic(p));
        assert!(matches!(proj, CrsCoord::Projected(_)), "EPSG:{code} is projected");
        let back = transform(crs, wgs84, proj).point();
        assert!(
            approx(back.lat, lat, 1e-6) && approx(back.lon, lon, 1e-6),
            "EPSG:{code} round-trip: back=({}, {}) dlat={:e} dlon={:e}",
            back.lat,
            back.lon,
            back.lat - lat,
            back.lon - lon
        );
    }
}

#[test]
fn area_of_use_bounds_and_containment() {
    // A UTM zone's area of use is its 6° band (zone 30N: lon −6..0, lat 0..84).
    let utm30n = Crs::from_epsg(32_630).unwrap();
    let band = utm30n.area_of_use().unwrap();
    assert_eq!((band.min_lon, band.max_lon, band.min_lat, band.max_lat), (-6.0, 0.0, 0.0, 84.0));
    assert!(utm30n.contains(GeoPoint::from_lat_lon(51.5, -0.1))); // London, in band
    assert!(!utm30n.contains(GeoPoint::from_lat_lon(51.5, 50.0))); // far east, out
    assert!(!utm30n.contains(GeoPoint::from_lat_lon(-5.0, -3.0))); // south of equator, out

    // A southern MGA zone bands into the southern hemisphere (lat −80..0).
    let mga56 = Crs::from_epsg(28_356).unwrap().area_of_use().unwrap();
    assert_eq!((mga56.min_lat, mga56.max_lat), (-80.0, 0.0));

    // An antimeridian-wrapping geographic extent (NAD83: lon 167.65 → −40.73).
    let nad83 = Crs::from_epsg(4269).unwrap();
    let aou = nad83.area_of_use().unwrap();
    assert!(aou.min_lon > aou.max_lon, "NAD83 area of use wraps the antimeridian");
    assert!(nad83.contains(GeoPoint::from_lat_lon(45.0, -100.0))); // CONUS (west side)
    assert!(nad83.contains(GeoPoint::from_lat_lon(52.0, 175.0))); // Aleutians (east side)
    assert!(!nad83.contains(GeoPoint::from_lat_lon(50.0, 10.0))); // Europe, outside the wrap

    // Polar CRSes guard on latitude.
    let nsidc_north = Crs::from_epsg(3413).unwrap();
    assert!(nsidc_north.contains(GeoPoint::from_lat_lon(75.0, -40.0)));
    assert!(!nsidc_north.contains(GeoPoint::from_lat_lon(40.0, -40.0)));

    // A codeless CRS has no declared area of use: it contains every finite point
    // but rejects non-finite ones.
    let local = Crs::local_tangent_plane(GeoPoint::from_lat_lon(52.0, 10.0));
    assert_eq!(local.area_of_use(), None);
    assert!(local.contains(GeoPoint::from_lat_lon(0.0, 0.0)));
    assert!(!local.contains(GeoPoint::from_lat_lon(f64::NAN, 0.0)));
}

#[test]
fn epsg_registry_recognises_3d_geographic() {
    // 4979 (WGS84 3D) and 4937 (ETRS89 3D) are geographic-3D CRSes.
    let wgs84_3d = Crs::from_epsg(4979).unwrap();
    assert_eq!(wgs84_3d, Crs::wgs84_3d());
    assert!(matches!(wgs84_3d, Crs::Geographic { dim: GeoDim::D3, .. }));
    assert_eq!(wgs84_3d.epsg(), Some(4979));
    let etrs89_3d = Crs::from_epsg(4937).unwrap();
    assert!(matches!(etrs89_3d, Crs::Geographic { dim: GeoDim::D3, .. }));
    assert_eq!(etrs89_3d.datum().ellipsoid, Ellipsoid::GRS80);
    // Both declare an area of use (global / ETRS89 extent).
    assert!(wgs84_3d.area_of_use().is_some() && etrs89_3d.area_of_use().is_some());
}

#[test]
fn geographic_3d_round_trips_height_and_2d_drops_it() {
    let wgs84_3d = Crs::wgs84_3d();
    let wgs84_2d = Crs::wgs84();
    let ecef = Crs::from_epsg(4978).unwrap();
    let p = GeoPoint::from_lat_lon(45.0, 45.0);
    // 3-D geographic → ECEF → 3-D geographic preserves lat/lon AND height.
    for &h in &[-200.0, 0.0, 1234.5, 8848.0] {
        let xyz = transform(wgs84_3d, ecef, CrsCoord::geographic_with_height(p, h));
        let CrsCoord::Geographic { point: back, height: hb } = transform(ecef, wgs84_3d, xyz)
        else {
            panic!()
        };
        assert!(approx(back.lat, p.lat, 1e-9) && approx(back.lon, p.lon, 1e-9));
        assert!(approx(hb, h, 1e-6), "height {h} → {hb}");
    }
    // A 2-D CRS (4326) drops the height: a 3-D input comes out at height 0.
    let CrsCoord::Geographic { height: h2d, .. } =
        transform(wgs84_3d, wgs84_2d, CrsCoord::geographic_with_height(p, 1000.0))
    else {
        panic!()
    };
    assert_eq!(h2d, 0.0);
    // The 2-D constructor is height 0.
    assert_eq!(CrsCoord::geographic(p).height(), 0.0);
}

#[test]
fn height_displaces_ecef_along_the_normal() {
    // Two coords at the same lat/lon, 1000 m apart in ellipsoidal height, map to
    // ECEF points exactly 1000 m apart — height is measured along the normal.
    let wgs84_3d = Crs::wgs84_3d();
    let ecef = Crs::from_epsg(4978).unwrap();
    let p = GeoPoint::from_lat_lon(45.0, 45.0);
    let a = transform(wgs84_3d, ecef, CrsCoord::geographic_with_height(p, 0.0)).geocentric();
    let b = transform(wgs84_3d, ecef, CrsCoord::geographic_with_height(p, 1000.0)).geocentric();
    let d = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
    assert!(approx(d, 1000.0, 1e-6), "‖Δecef‖ = {d}");
}

#[test]
fn datum_shift_carries_height_in_3d() {
    // A 3-D geographic CRS on the OSGB36 datum (Airy 1830 + Helmert). The Helmert
    // acts in 3-D (on ECEF), so a point's height survives the datum shift and
    // reverses exactly (the proven Helmert inverse).
    let osgb36_3d = Crs::Geographic {
        datum: Datum::osgb36(),
        epsg: None,
        dim: GeoDim::D3,
    };
    let wgs84_3d = Crs::wgs84_3d();
    let p = GeoPoint::from_lat_lon(51.5, -0.1);
    let h = 142.0;
    let CrsCoord::Geographic { point: shifted, height: hs } =
        transform(osgb36_3d, wgs84_3d, CrsCoord::geographic_with_height(p, h))
    else {
        panic!()
    };
    // The horizontal shift is real (datum change), and the height is carried
    // (not zeroed), differing only by the ellipsoid/Helmert (< ~100 m for the UK).
    assert!(
        (shifted.lat - p.lat).abs() > 1e-5 || (shifted.lon - p.lon).abs() > 1e-5,
        "datum shift should be nonzero"
    );
    assert!((hs - h).abs() < 200.0, "height {h} → {hs} across datum");
    // Round-trips exactly back to the OSGB36 lat/lon AND height.
    let CrsCoord::Geographic { point: back, height: hb } =
        transform(wgs84_3d, osgb36_3d, CrsCoord::geographic_with_height(shifted, hs))
    else {
        panic!()
    };
    assert!(approx(back.lat, p.lat, 1e-7) && approx(back.lon, p.lon, 1e-7));
    assert!(approx(hb, h, 1e-6), "height round-trip {h} → {hb}");
}

#[test]
fn compound_crs_collapses_to_geographic_3d() {
    // Geographic-2D horizontal + ellipsoidal vertical ≡ geographic-3D.
    let compound = CompoundCrs::new(Crs::wgs84(), VerticalCrs::Ellipsoidal { epsg: None });
    let equiv = compound.as_crs().expect("geographic + ellipsoidal collapses");
    assert!(matches!(equiv, Crs::Geographic { dim: GeoDim::D3, .. }));
    assert_eq!(equiv.datum(), Datum::WGS84);
    // It flows through `transform` like any 3-D geographic CRS — height retained.
    let ecef = Crs::from_epsg(4978).unwrap();
    let p = GeoPoint::from_lat_lon(40.0, -3.0);
    let c = CrsCoord::geographic_with_height(p, 655.0);
    let CrsCoord::Geographic { height: hb, .. } =
        transform(ecef, equiv, transform(equiv, ecef, c))
    else {
        panic!()
    };
    assert!(approx(hb, 655.0, 1e-6));
    // A projected-horizontal compound doesn't collapse yet (deferred).
    let projected =
        CompoundCrs::new(Crs::from_epsg(3857).unwrap(), VerticalCrs::Ellipsoidal { epsg: None });
    assert!(projected.as_crs().is_none());
}

#[test]
fn proj_string_parses_to_registry_equivalent() {
    // (EPSG code, its PROJ-string from pyproj, a test lat/lon). The parsed CRS
    // must project identically to the registry CRS. Both are exercised on their
    // *own* datum (PROJ-strings carry no datum shift), isolating projection +
    // ellipsoid.
    let cases: &[(u32, &str, f64, f64)] = &[
        (4326, "+proj=longlat +datum=WGS84 +no_defs +type=crs", 51.5, -0.1),
        (4978, "+proj=geocent +datum=WGS84 +units=m +no_defs", 45.0, 45.0),
        (
            3857,
            "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 \
             +k=1 +units=m +nadgrids=@null +wktext +no_defs",
            51.5,
            -0.1,
        ),
        (32_630, "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs", 51.5, -0.1),
        (
            27_700,
            "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 \
             +y_0=-100000 +ellps=airy +units=m +no_defs",
            52.0,
            -1.0,
        ),
        (
            2154,
            "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 \
             +y_0=6600000 +ellps=GRS80 +units=m +no_defs",
            46.5,
            2.5,
        ),
        (
            5070,
            "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 \
             +datum=NAD83 +units=m +no_defs",
            40.0,
            -100.0,
        ),
        (
            3031,
            "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 \
             +datum=WGS84 +units=m +no_defs",
            -75.0,
            30.0,
        ),
    ];
    for &(code, s, lat, lon) in cases {
        let from_code = Crs::from_epsg(code).unwrap();
        let from_str = Crs::from_proj_string(s).unwrap_or_else(|| panic!("parse failed: {s}"));
        assert_eq!(from_str.epsg(), None, "EPSG:{code}: a parsed CRS carries no code");
        let geo = |c: Crs| Crs::Geographic { datum: c.datum(), epsg: None, dim: GeoDim::D2 };
        let p = CrsCoord::geographic(GeoPoint::from_lat_lon(lat, lon));
        let a = transform(geo(from_code), from_code, p);
        let b = transform(geo(from_str), from_str, p);
        let d = match (a, b) {
            (CrsCoord::Projected(a), CrsCoord::Projected(b)) => (a.x - b.x).abs().max((a.y - b.y).abs()),
            (CrsCoord::Geocentric(a), CrsCoord::Geocentric(b)) => {
                (a.x - b.x).abs().max((a.y - b.y).abs()).max((a.z - b.z).abs())
            }
            (CrsCoord::Geographic { point: a, .. }, CrsCoord::Geographic { point: b, .. }) => {
                (a.lat - b.lat).abs().max((a.lon - b.lon).abs())
            }
            _ => panic!("EPSG:{code}: parsed CRS is a different kind"),
        };
        assert!(d < 1e-6, "EPSG:{code}: parsed CRS differs from registry by {d}");
    }
}

#[test]
fn proj_string_init_and_towgs84_and_rejects() {
    // `+init=epsg:NNNN` defers to the registry, keeping the code.
    assert_eq!(Crs::from_proj_string("+init=epsg:4326"), Some(Crs::wgs84()));
    assert_eq!(Crs::from_proj_string("+init=epsg:32630").unwrap().epsg(), Some(32_630));
    // `+towgs84=` builds a non-identity Helmert datum shift.
    let ed50 = Crs::from_proj_string(
        "+proj=longlat +ellps=intl +towgs84=-87,-98,-121,0,0,0,0 +no_defs",
    )
    .unwrap();
    assert!(matches!(
        ed50.datum().to_wgs84,
        DatumShift::Helmert(h) if !h.is_identity()
    ));
    // A bare `+ellps` with no shift is the identity (own-datum) Helmert.
    let airy = Crs::from_proj_string("+proj=longlat +ellps=airy +no_defs").unwrap();
    assert!(matches!(airy.datum().to_wgs84, DatumShift::Helmert(h) if h.is_identity()));
    // Unsupported projection / garbage / oblique stereographic → None.
    assert_eq!(Crs::from_proj_string("+proj=robin +datum=WGS84"), None);
    assert_eq!(Crs::from_proj_string("not a proj string"), None);
    assert_eq!(Crs::from_proj_string("+proj=stere +lat_0=45 +datum=WGS84"), None);
}
