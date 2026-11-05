use super::super::*;

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
