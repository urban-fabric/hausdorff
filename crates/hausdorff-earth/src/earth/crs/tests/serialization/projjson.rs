use super::super::*;

#[test]
fn projjson_round_trips_and_parses() {
    // Self round-trip: emit PROJJSON, parse it back, project identically. Spans
    // geographic 2-D/3-D, geocentric, and every projection family.
    for &(code, lat, lon) in &[
        (4326u32, 51.5, -0.1),
        (4979, 51.5, -0.1),
        (4978, 45.0, 45.0),
        (3857, 51.5, -0.1),
        (32_630, 51.5, -0.1),
        (27_700, 52.0, -1.0),
        (2154, 46.5, 2.5),
        (5070, 40.0, -100.0),
        (3031, -75.0, 30.0),
    ] {
        let crs = Crs::from_epsg(code).unwrap();
        let j = crs.to_projjson();
        let back = Crs::from_projjson(&j)
            .unwrap_or_else(|| panic!("EPSG:{code}: PROJJSON re-parse failed: {j}"));
        assert!(
            super::projection_delta(crs, back, lat, lon) < 1e-6,
            "EPSG:{code}: PROJJSON round-trip off — {j}"
        );
        assert_eq!(
            back.epsg(),
            Some(code),
            "EPSG:{code}: PROJJSON keeps the code"
        );
    }
}

#[test]
fn projjson_from_epsg_carries_names_id_and_round_trips() {
    // NAD83 (TIGER's datum): the wrapper yields PROJJSON with the real names, the
    // correct GRS 1980 ellipsoid, and the EPSG id readers key on.
    let json = Crs::projjson_from_epsg(4269).expect("4269 is in the registry");
    assert!(json.contains(r#""name":"NAD83""#), "{json}");
    assert!(
        json.contains(r#""name":"North American Datum 1983""#),
        "{json}"
    );
    assert!(json.contains(r#""name":"GRS 1980""#), "{json}");
    assert!(
        json.contains(r#""inverse_flattening":298.257222101"#),
        "{json}"
    );
    assert!(
        json.contains(r#""id":{"authority":"EPSG","code":4269}"#),
        "{json}"
    );
    // Names are cosmetic — it still round-trips back to the same CRS by id.
    assert_eq!(
        Crs::from_projjson(&json),
        Some(Crs::from_epsg(4269).unwrap())
    );

    // WGS 84 likewise; an unknown code yields None.
    assert!(
        Crs::projjson_from_epsg(4326)
            .unwrap()
            .contains(r#""name":"WGS 84""#)
    );
    assert_eq!(Crs::projjson_from_epsg(999_999), None);
}

#[test]
fn projjson_parses_hand_written_and_rejects() {
    // Geographic by EPSG id (registry datum); 2 axes → 2-D.
    let g = Crs::from_projjson(
        r#"{"type":"GeographicCRS","name":"WGS 84","datum_ensemble":{"name":"e","ellipsoid":{"name":"WGS 84","semi_major_axis":6378137,"inverse_flattening":298.257223563}},"coordinate_system":{"subtype":"ellipsoidal","axis":[{"name":"lat","direction":"north"},{"name":"lon","direction":"east"}]},"id":{"authority":"EPSG","code":4326}}"#,
    )
    .unwrap();
    assert_eq!(g, Crs::wgs84());
    // 3 axes → 3-D geographic.
    let g3 = Crs::from_projjson(
        r#"{"type":"GeographicCRS","name":"x","datum":{"name":"d","ellipsoid":{"name":"WGS 84","semi_major_axis":6378137,"inverse_flattening":298.257223563}},"coordinate_system":{"subtype":"ellipsoidal","axis":[{},{},{}]},"id":{"authority":"EPSG","code":4979}}"#,
    )
    .unwrap();
    assert!(matches!(
        g3,
        Crs::Geographic {
            dim: GeoDim::D3,
            ..
        }
    ));
    // Geocentric (Cartesian CS).
    let geo = Crs::from_projjson(
        r#"{"type":"GeodeticCRS","name":"x","datum":{"name":"d","ellipsoid":{"name":"WGS 84","semi_major_axis":6378137,"inverse_flattening":298.257223563}},"coordinate_system":{"subtype":"Cartesian","axis":[{},{},{}]},"id":{"authority":"EPSG","code":4978}}"#,
    )
    .unwrap();
    assert!(matches!(geo, Crs::Geocentric { .. }));
    // Projected whose base id is NOT in our registry → ELLIPSOID fallback (GRS80),
    // and the projection still matches the registry Lambert-93.
    let lcc = Crs::from_projjson(
        r#"{"type":"ProjectedCRS","name":"x","base_crs":{"name":"RGF93","datum":{"name":"d","ellipsoid":{"name":"GRS 1980","semi_major_axis":6378137,"inverse_flattening":298.257222101}},"id":{"authority":"EPSG","code":4171}},"conversion":{"name":"c","method":{"name":"Lambert Conic Conformal (2SP)"},"parameters":[{"name":"Latitude of false origin","value":46.5},{"name":"Longitude of false origin","value":3},{"name":"Latitude of 1st standard parallel","value":49},{"name":"Latitude of 2nd standard parallel","value":44},{"name":"Easting at false origin","value":700000},{"name":"Northing at false origin","value":6600000}]},"coordinate_system":{"subtype":"Cartesian","axis":[{},{}]}}"#,
    )
    .unwrap();
    assert!(super::projection_delta(Crs::from_epsg(2154).unwrap(), lcc, 46.5, 2.5) < 1e-6);
    // Malformed / unsupported → None.
    assert_eq!(Crs::from_projjson("{not json"), None);
    assert_eq!(Crs::from_projjson("[]"), None);
    assert_eq!(Crs::from_projjson(r#"{"type":"VerticalCRS"}"#), None);
}

#[test]
fn serialized_datum_shift_round_trips() {
    // OSGB36 / BNG (27700) carries a real Helmert shift; both serializers preserve
    // it — proj-string via +towgs84, WKT2 via the BASEGEOGCRS EPSG id — so a full
    // WGS84→BNG transform through the re-parsed CRS matches the registry.
    let reg = Crs::from_epsg(27_700).unwrap();
    let p = CrsCoord::geographic(GeoPoint::from_lat_lon(52.0, -1.0));
    let want = transform(Crs::wgs84(), reg, p).projected();
    for parsed in [
        Crs::from_proj_string(&reg.to_proj_string()).unwrap(),
        Crs::from_wkt2(&reg.to_wkt2()).unwrap(),
    ] {
        // The shift survived (non-identity), not silently dropped to projection-only.
        assert!(matches!(parsed.datum().to_wgs84, DatumShift::Helmert(h) if !h.is_identity()));
        let got = transform(Crs::wgs84(), parsed, p).projected();
        assert!(
            (got.x - want.x).abs() < 1e-6 && (got.y - want.y).abs() < 1e-6,
            "datum shift not preserved: {got:?} vs {want:?}"
        );
    }
}
