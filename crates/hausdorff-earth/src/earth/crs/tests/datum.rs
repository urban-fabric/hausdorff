use super::*;

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
        assert!(
            dlat > 1e-5 || dlon > 1e-5,
            "{:?}: shift too small",
            crs.epsg()
        );
        assert!(
            dlat < 1e-2 && dlon < 1e-2,
            "{:?}: shift too large",
            crs.epsg()
        );
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
    // EPSG code, label, and `(lat, lon, expected_lat, expected_lon)` samples.
    type DatumCase<'a> = (u32, &'a str, &'a [(f64, f64, f64, f64)]);
    let cases: &[DatumCase] = &[
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
            let CrsCoord::Geographic { point: g, .. } = transform(crs, wgs84, p) else {
                panic!()
            };
            max_deg = max_deg.max((g.lat - rlat).abs()).max((g.lon - rlon).abs());
            // Round-trips: source → WGS84 → source → WGS84.
            let back = transform(wgs84, crs, transform(crs, wgs84, p));
            let CrsCoord::Geographic { point: b, .. } = transform(crs, wgs84, back) else {
                panic!()
            };
            assert!((b.lat - rlat).abs() < 1e-7 && (b.lon - rlon).abs() < 1e-7);
        }
        eprintln!("{name}→WGS84 (grid) end-to-end max Δ = {max_deg:.3e}° vs PROJ");
        assert!(
            max_deg < 1e-7,
            "{name} grid datum pipeline off by {max_deg}°"
        );
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
