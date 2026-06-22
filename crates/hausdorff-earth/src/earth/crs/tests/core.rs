use super::*;

#[test]
fn epsg_registry_recognises_the_core_systems() {
    assert_eq!(Crs::from_epsg(4326), Some(Crs::wgs84()));
    assert!(matches!(Crs::from_epsg(4978), Some(Crs::Geocentric { .. })));
    assert!(matches!(Crs::from_epsg(3857), Some(Crs::Projected { .. })));
    // UTM zone 30N / 30S.
    assert_eq!(Crs::from_epsg(32_630).unwrap().epsg(), Some(32_630));
    assert!(matches!(
        Crs::from_epsg(32_730),
        Some(Crs::Projected { .. })
    ));
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
    assert!(
        matches!(origin, CrsCoord::Projected(c) if approx(c.x, 0.0, 1e-6) && approx(c.y, 0.0, 1e-6))
    );
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
        assert!(
            matches!(proj, CrsCoord::Projected(_)),
            "EPSG:{code} is projected"
        );
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
