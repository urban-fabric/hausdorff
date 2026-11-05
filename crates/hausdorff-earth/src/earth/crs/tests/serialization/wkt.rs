use super::super::*;

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
        let d = super::projection_delta(from_code, from_text, lat, lon);
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
