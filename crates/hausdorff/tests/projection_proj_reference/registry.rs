use super::*;

/// Forward+inverse for a **registry** projected CRS, driving it the way callers
/// do — `Crs::from_epsg` then [`transform`] from the CRS's own-datum geographic
/// base (projection only, no datum shift). Validates the registry wiring and the
/// projection together vs PROJ.
fn check_registry(code: u32, refs: &[Ref], fwd_tol_m: f64) {
    let proj = Crs::from_epsg(code).unwrap_or_else(|| panic!("EPSG:{code} unrecognised"));
    let geo = Crs::Geographic {
        datum: proj.datum(),
        epsg: None,
        dim: GeoDim::D2,
    };
    let (mut max_e, mut max_n, mut max_inv) = (0.0_f64, 0.0_f64, 0.0_f64);
    for &(lat, lon, e, n) in refs {
        let CrsCoord::Projected(got) = transform(
            geo,
            proj,
            CrsCoord::geographic(GeoPoint::from_lat_lon(lat, lon)),
        ) else {
            panic!("EPSG:{code} should project")
        };
        max_e = max_e.max((got.x - e).abs());
        max_n = max_n.max((got.y - n).abs());
        let CrsCoord::Geographic { point: back, .. } =
            transform(proj, geo, CrsCoord::Projected(Coord::new(e, n)))
        else {
            panic!("EPSG:{code} inverse should be geographic")
        };
        max_inv = max_inv
            .max((back.lat - lat).abs())
            .max((back.lon - lon).abs());
    }
    eprintln!(
        "EPSG:{code:<6} forward max |ΔE|={max_e:.6} m |ΔN|={max_n:.6} m   inverse max Δ={max_inv:.3e}°"
    );
    assert!(
        max_e < fwd_tol_m && max_n < fwd_tol_m,
        "EPSG:{code}: forward ΔE={max_e} ΔN={max_n} exceeds {fwd_tol_m} m"
    );
    assert!(
        max_inv < 1e-6,
        "EPSG:{code}: inverse Δ {max_inv}° too large"
    );
}

/// The broadened registry's new projected CRSes vs PROJ 9.4.0 (each from its own
/// GRS80/WGS84 geographic base, in-zone points for the Transverse Mercator ones).
/// `(epsg, &[(lat, lon, easting, northing)])` — a dense reference table laid out
/// one CRS per row (kept verbatim with `rustfmt::skip`).
#[rustfmt::skip]
#[test]
fn broadened_registry_projections_match_proj() {
    let cases: &[(u32, &[Ref])] = &[
        // ETRS89 / NAD83 / GDA94 UTM-family zones (Transverse Mercator on GRS80).
        (25_828, &[(28.1, -15.4, 460_707.1786, 3_108_344.5672), (40.0, -16.0, 414_639.5382, 4_428_236.0645)]),
        (25_838, &[(40.0, 45.0, 500_000.0000, 4_427_757.2186), (38.0, 44.0, 412_201.5822, 4_206_286.7579)]),
        (26_901, &[(60.0, -176.0, 555_776.2668, 6_651_832.7353), (52.0, -174.5, 671_613.2238, 5_763_989.3666)]),
        (26_923, &[(45.0, -45.0, 500_000.0000, 4_982_950.4001), (50.0, -48.0, 285_015.7633, 5_542_944.0185)]),
        (28_348, &[(-30.0, 105.0, 500_000.0000, 6_681_214.6475), (-25.0, 103.5, 348_624.2414, 7_234_214.7049)]),
        (28_358, &[(-43.0, 165.0, 500_000.0000, 5_239_185.2040), (-30.0, 166.5, 644_679.8540, 6_680_267.5834)]),
        // NZGD2000 / NZTM2000 and ETRS89 / Poland CS92 (k0 = 0.9993).
        (2193, &[(-41.29, 174.78, 1_749_045.8029, 5_427_521.3834), (-36.85, 174.76, 1_756_911.9449, 5_920_321.8182), (-43.53, 172.64, 1_570_910.4458, 5_180_263.9268)]),
        (2180, &[(52.23, 21.01, 637_231.0903, 486_786.3937), (50.06, 19.94, 567_262.4529, 244_060.6145)]),
        // NAD83 / Canada Atlas Lambert (LCC) and GDA94 / Australian Albers.
        (3978, &[(60.0, -95.0, 0.0000, 1_204_290.9001), (45.42, -75.7, 1_510_453.2007, -170_038.9917)]),
        (3577, &[(-33.87, 151.21, 1_761_025.1670, -3_828_891.3656), (-25.0, 133.0, 99_795.5592, -2_702_961.5341)]),
        // WGS84 / NSIDC Sea Ice Polar Stereographic North and South.
        (3413, &[(75.0, -40.0, 142_401.9812, -1_627_662.0927), (80.0, -60.0, -281_056.8544, -1_048_918.4605)]),
        (3976, &[(-75.0, 0.0, 0.0000, 1_633_879.4975), (-80.0, 120.0, 940_434.5640, -542_960.1487)]),
    ];
    for &(code, refs) in cases {
        check_registry(code, refs, 5e-3);
    }
}
