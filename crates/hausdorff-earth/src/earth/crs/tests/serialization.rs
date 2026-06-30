use super::*;

#[test]
fn proj_string_parses_to_registry_equivalent() {
    // (EPSG code, its PROJ-string from pyproj, a test lat/lon). The parsed CRS
    // must project identically to the registry CRS. Both are exercised on their
    // *own* datum (PROJ-strings carry no datum shift), isolating projection +
    // ellipsoid.
    let cases: &[(u32, &str, f64, f64)] = &[
        (
            4326,
            "+proj=longlat +datum=WGS84 +no_defs +type=crs",
            51.5,
            -0.1,
        ),
        (
            4978,
            "+proj=geocent +datum=WGS84 +units=m +no_defs",
            45.0,
            45.0,
        ),
        (
            3857,
            "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 \
             +k=1 +units=m +nadgrids=@null +wktext +no_defs",
            51.5,
            -0.1,
        ),
        (
            32_630,
            "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs",
            51.5,
            -0.1,
        ),
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
        assert_eq!(
            from_str.epsg(),
            None,
            "EPSG:{code}: a parsed CRS carries no code"
        );
        let geo = |c: Crs| Crs::Geographic {
            datum: c.datum(),
            epsg: None,
            dim: GeoDim::D2,
        };
        let p = CrsCoord::geographic(GeoPoint::from_lat_lon(lat, lon));
        let a = transform(geo(from_code), from_code, p);
        let b = transform(geo(from_str), from_str, p);
        let d = match (a, b) {
            (CrsCoord::Projected(a), CrsCoord::Projected(b)) => {
                (a.x - b.x).abs().max((a.y - b.y).abs())
            }
            (CrsCoord::Geocentric(a), CrsCoord::Geocentric(b)) => (a.x - b.x)
                .abs()
                .max((a.y - b.y).abs())
                .max((a.z - b.z).abs()),
            (CrsCoord::Geographic { point: a, .. }, CrsCoord::Geographic { point: b, .. }) => {
                (a.lat - b.lat).abs().max((a.lon - b.lon).abs())
            }
            _ => panic!("EPSG:{code}: parsed CRS is a different kind"),
        };
        assert!(
            d < 1e-6,
            "EPSG:{code}: parsed CRS differs from registry by {d}"
        );
    }
}

#[test]
fn proj_string_init_and_towgs84_and_rejects() {
    // `+init=epsg:NNNN` defers to the registry, keeping the code.
    assert_eq!(Crs::from_proj_string("+init=epsg:4326"), Some(Crs::wgs84()));
    assert_eq!(
        Crs::from_proj_string("+init=epsg:32630").unwrap().epsg(),
        Some(32_630)
    );
    // `+towgs84=` builds a non-identity Helmert datum shift.
    let ed50 =
        Crs::from_proj_string("+proj=longlat +ellps=intl +towgs84=-87,-98,-121,0,0,0,0 +no_defs")
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
    assert_eq!(
        Crs::from_proj_string("+proj=stere +lat_0=45 +datum=WGS84"),
        None
    );
}

/// Compare a parsed CRS to the registry CRS by projecting a test point on each
/// CRS's *own* datum (isolating projection + ellipsoid), returning the worst
/// component difference. Panics if the two CRSes are different kinds.
fn projection_delta(from_code: Crs, from_text: Crs, lat: f64, lon: f64) -> f64 {
    let geo = |c: Crs| Crs::Geographic {
        datum: c.datum(),
        epsg: None,
        dim: GeoDim::D2,
    };
    let p = CrsCoord::geographic(GeoPoint::from_lat_lon(lat, lon));
    let a = transform(geo(from_code), from_code, p);
    let b = transform(geo(from_text), from_text, p);
    match (a, b) {
        (CrsCoord::Projected(a), CrsCoord::Projected(b)) => {
            (a.x - b.x).abs().max((a.y - b.y).abs())
        }
        (CrsCoord::Geocentric(a), CrsCoord::Geocentric(b)) => (a.x - b.x)
            .abs()
            .max((a.y - b.y).abs())
            .max((a.z - b.z).abs()),
        (CrsCoord::Geographic { point: a, .. }, CrsCoord::Geographic { point: b, .. }) => {
            (a.lat - b.lat).abs().max((a.lon - b.lon).abs())
        }
        _ => panic!("parsed CRS is a different kind than the registry CRS"),
    }
}

#[test]
fn wkt2_parses_to_registry_equivalent() {
    // (EPSG code, its WKT2:2019 string from pyproj, test lat/lon). Covers
    // GEOGCRS+ENSEMBLE, GEODCRS (geocentric), and PROJCRS with TM (known base
    // 4277 → registry datum), LCC (base 4171 not in our registry → ELLIPSOID
    // fallback), Albers, and Polar Stereographic (variant B).
    let cases: &[(u32, &str, f64, f64)] = &[
        (
            4326,
            r#"GEOGCRS["WGS 84",ENSEMBLE["World Geodetic System 1984 ensemble",MEMBER["World Geodetic System 1984 (Transit)"],ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ENSEMBLEACCURACY[2.0]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],CS[ellipsoidal,2],AXIS["geodetic latitude (Lat)",north,ORDER[1],ANGLEUNIT["degree",0.0174532925199433]],AXIS["geodetic longitude (Lon)",east,ORDER[2],ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4326]]"#,
            51.5,
            -0.1,
        ),
        (
            4978,
            r#"GEODCRS["WGS 84",ENSEMBLE["World Geodetic System 1984 ensemble",MEMBER["World Geodetic System 1984 (Transit)"],ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ENSEMBLEACCURACY[2.0]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],CS[Cartesian,3],AXIS["(X)",geocentricX,ORDER[1],LENGTHUNIT["metre",1]],AXIS["(Y)",geocentricY,ORDER[2],LENGTHUNIT["metre",1]],AXIS["(Z)",geocentricZ,ORDER[3],LENGTHUNIT["metre",1]],ID["EPSG",4978]]"#,
            45.0,
            45.0,
        ),
        (
            27_700,
            r#"PROJCRS["OSGB36 / British National Grid",BASEGEOGCRS["OSGB36",DATUM["Ordnance Survey of Great Britain 1936",ELLIPSOID["Airy 1830",6377563.396,299.3249646,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4277]],CONVERSION["British National Grid",METHOD["Transverse Mercator",ID["EPSG",9807]],PARAMETER["Latitude of natural origin",49,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8801]],PARAMETER["Longitude of natural origin",-2,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["Scale factor at natural origin",0.9996012717,SCALEUNIT["unity",1],ID["EPSG",8805]],PARAMETER["False easting",400000,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",-100000,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",27700]]"#,
            52.0,
            -1.0,
        ),
        (
            2154,
            r#"PROJCRS["RGF93 v1 / Lambert-93",BASEGEOGCRS["RGF93 v1",DATUM["Reseau Geodesique Francais 1993 v1",ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4171]],CONVERSION["Lambert-93",METHOD["Lambert Conic Conformal (2SP)",ID["EPSG",9802]],PARAMETER["Latitude of false origin",46.5,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8821]],PARAMETER["Longitude of false origin",3,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8822]],PARAMETER["Latitude of 1st standard parallel",49,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Latitude of 2nd standard parallel",44,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8824]],PARAMETER["Easting at false origin",700000,LENGTHUNIT["metre",1],ID["EPSG",8826]],PARAMETER["Northing at false origin",6600000,LENGTHUNIT["metre",1],ID["EPSG",8827]]],CS[Cartesian,2],AXIS["easting (X)",east,ORDER[1],LENGTHUNIT["metre",1]],AXIS["northing (Y)",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",2154]]"#,
            46.5,
            2.5,
        ),
        (
            5070,
            r#"PROJCRS["NAD83 / Conus Albers",BASEGEOGCRS["NAD83",DATUM["North American Datum 1983",ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4269]],CONVERSION["Conus Albers",METHOD["Albers Equal Area",ID["EPSG",9822]],PARAMETER["Latitude of false origin",23,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8821]],PARAMETER["Longitude of false origin",-96,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8822]],PARAMETER["Latitude of 1st standard parallel",29.5,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Latitude of 2nd standard parallel",45.5,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8824]],PARAMETER["Easting at false origin",0,LENGTHUNIT["metre",1],ID["EPSG",8826]],PARAMETER["Northing at false origin",0,LENGTHUNIT["metre",1],ID["EPSG",8827]]],CS[Cartesian,2],AXIS["easting (X)",east,ORDER[1],LENGTHUNIT["metre",1]],AXIS["northing (Y)",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",5070]]"#,
            40.0,
            -100.0,
        ),
        (
            3031,
            r#"PROJCRS["WGS 84 / Antarctic Polar Stereographic",BASEGEOGCRS["WGS 84",ENSEMBLE["World Geodetic System 1984 ensemble",MEMBER["World Geodetic System 1984 (Transit)"],ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ENSEMBLEACCURACY[2.0]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4326]],CONVERSION["Antarctic Polar Stereographic",METHOD["Polar Stereographic (variant B)",ID["EPSG",9829]],PARAMETER["Latitude of standard parallel",-71,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8832]],PARAMETER["Longitude of origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8833]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",north,MERIDIAN[90,ANGLEUNIT["degree",0.0174532925199433]],ORDER[1],LENGTHUNIT["metre",1]],AXIS["(N)",north,MERIDIAN[0,ANGLEUNIT["degree",0.0174532925199433]],ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",3031]]"#,
            -75.0,
            30.0,
        ),
    ];
    for &(code, wkt, lat, lon) in cases {
        let from_code = Crs::from_epsg(code).unwrap();
        let from_text =
            Crs::from_wkt2(wkt).unwrap_or_else(|| panic!("WKT2 parse failed for EPSG:{code}"));
        // WKT2 carries the ID, so the parsed CRS keeps its EPSG code.
        assert_eq!(
            from_text.epsg(),
            Some(code),
            "EPSG:{code}: code preserved from ID[]"
        );
        let d = projection_delta(from_code, from_text, lat, lon);
        assert!(
            d < 1e-6,
            "EPSG:{code}: parsed WKT2 CRS differs from registry by {d}"
        );
    }
}

#[test]
fn wkt2_dimension_datum_and_rejects() {
    // CS[ellipsoidal,3] → a 3-D geographic CRS (height retained).
    let g3d = Crs::from_wkt2(
        r#"GEOGCRS["WGS 84 (3D)",DATUM["WGS84",ELLIPSOID["WGS 84",6378137,298.257223563]],CS[ellipsoidal,3],AXIS["lat",north],AXIS["lon",east],AXIS["h",up],ID["EPSG",4979]]"#,
    )
    .unwrap();
    assert!(matches!(
        g3d,
        Crs::Geographic {
            dim: GeoDim::D3,
            ..
        }
    ));
    assert_eq!(g3d.epsg(), Some(4979));
    // A known base code pulls the registry datum (OSGB36's Helmert shift), so a
    // full WGS84→BNG transform matches the registry end to end (not just the
    // projection). Caister worked example, ~0.5 m.
    let bng = Crs::from_wkt2(
        r#"PROJCRS["x",BASEGEOGCRS["OSGB36",DATUM["d",ELLIPSOID["Airy 1830",6377563.396,299.3249646]],ID["EPSG",4277]],CONVERSION["c",METHOD["Transverse Mercator"],PARAMETER["Latitude of natural origin",49],PARAMETER["Longitude of natural origin",-2],PARAMETER["Scale factor at natural origin",0.9996012717],PARAMETER["False easting",400000],PARAMETER["False northing",-100000]],CS[Cartesian,2],AXIS["e",east],AXIS["n",north]]"#,
    )
    .unwrap();
    // The base's Helmert shift is non-identity (the registry OSGB36 datum, not a
    // bare ellipsoid), and a full WGS84→BNG transform matches the registry CRS
    // exactly — datum shift included.
    assert!(matches!(bng.datum().to_wgs84, DatumShift::Helmert(h) if !h.is_identity()));
    let p = CrsCoord::geographic(GeoPoint::from_lat_lon(52.657_570_3, 1.717_921_6));
    let a = transform(Crs::wgs84(), bng, p).projected();
    let b = transform(Crs::wgs84(), Crs::from_epsg(27_700).unwrap(), p).projected();
    assert!(
        (a.x - b.x).abs() < 1e-6 && (a.y - b.y).abs() < 1e-6,
        "WKT2 BNG vs registry: {a:?} vs {b:?}"
    );
    // Malformed / unsupported → None.
    assert_eq!(Crs::from_wkt2("GEOGCRS[\"x\","), None);
    assert_eq!(Crs::from_wkt2("not wkt"), None);
    assert_eq!(Crs::from_wkt2(r#"VERTCRS["height",VDATUM["d"]]"#), None);
}

#[test]
fn serializers_round_trip_through_the_parsers() {
    // Each registry CRS emitted to PROJ-string and WKT2, then re-parsed, must
    // project identically. (Spans geographic 2-D/3-D, geocentric, and every
    // projection family.)
    let cases: &[(u32, f64, f64)] = &[
        (4326, 51.5, -0.1),
        (4979, 51.5, -0.1),
        (4978, 45.0, 45.0),
        (4269, 40.0, -100.0),
        (3857, 51.5, -0.1),
        (32_630, 51.5, -0.1),
        (27_700, 52.0, -1.0),
        (2154, 46.5, 2.5),
        (5070, 40.0, -100.0),
        (3031, -75.0, 30.0),
    ];
    for &(code, lat, lon) in cases {
        let crs = Crs::from_epsg(code).unwrap();
        let ps = crs.to_proj_string();
        let from_ps = Crs::from_proj_string(&ps)
            .unwrap_or_else(|| panic!("EPSG:{code}: re-parse of proj-string failed: {ps}"));
        assert!(
            projection_delta(crs, from_ps, lat, lon) < 1e-6,
            "EPSG:{code}: proj-string round-trip off — {ps}"
        );
        let w = crs.to_wkt2();
        let from_w = Crs::from_wkt2(&w)
            .unwrap_or_else(|| panic!("EPSG:{code}: re-parse of WKT2 failed: {w}"));
        assert!(
            projection_delta(crs, from_w, lat, lon) < 1e-6,
            "EPSG:{code}: WKT2 round-trip off — {w}"
        );
        // WKT2 carries the ID, so the code survives the round-trip.
        assert_eq!(
            from_w.epsg(),
            Some(code),
            "EPSG:{code}: WKT2 keeps the code"
        );
    }
}

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
            projection_delta(crs, back, lat, lon) < 1e-6,
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
    assert!(json.contains(r#""name":"North American Datum 1983""#), "{json}");
    assert!(json.contains(r#""name":"GRS 1980""#), "{json}");
    assert!(json.contains(r#""inverse_flattening":298.257222101"#), "{json}");
    assert!(json.contains(r#""id":{"authority":"EPSG","code":4269}"#), "{json}");
    // Names are cosmetic — it still round-trips back to the same CRS by id.
    assert_eq!(Crs::from_projjson(&json), Some(Crs::from_epsg(4269).unwrap()));

    // WGS 84 likewise; an unknown code yields None.
    assert!(Crs::projjson_from_epsg(4326).unwrap().contains(r#""name":"WGS 84""#));
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
    assert!(projection_delta(Crs::from_epsg(2154).unwrap(), lcc, 46.5, 2.5) < 1e-6);
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
