use super::*;

#[test]
fn broadened_registry_recognises_the_new_codes() {
    // Geographic bases on GRS80 (NAD83 / ETRS89).
    for code in [4269_u32, 4258] {
        assert!(matches!(Crs::from_epsg(code), Some(Crs::Geographic { .. })));
        assert_eq!(Crs::from_epsg(code).unwrap().epsg(), Some(code));
        assert_eq!(
            Crs::from_epsg(code).unwrap().datum().ellipsoid,
            Ellipsoid::GRS80
        );
    }
    // A projected representative from every new family, by EPSG code.
    for code in [
        25_828_u32, 25_838, 26_901, 26_923, 28_348, 28_358, 2193, 2180, 3978, 3577, 3413, 3976,
    ] {
        let crs = Crs::from_epsg(code).unwrap_or_else(|| panic!("EPSG:{code} unrecognised"));
        assert!(
            matches!(crs, Crs::Projected { .. }),
            "EPSG:{code} should be projected"
        );
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
        (25_828_u32, 28.1, -15.4), // ETRS89 / UTM 28N (Canaries)
        (26_923, 45.0, -45.0),     // NAD83 / UTM 23N
        (28_356, -33.87, 151.21),  // GDA94 / MGA 56 (Sydney)
        (2193, -41.29, 174.78),    // NZTM2000 (Wellington)
        (2180, 52.23, 21.01),      // ETRS89 / Poland CS92 (Warsaw)
        (3978, 60.0, -95.0),       // Canada Atlas Lambert
        (3577, -25.0, 133.0),      // Australian Albers
        (3413, 75.0, -40.0),       // NSIDC Polar Stereographic North
        (3976, -75.0, 0.0),        // NSIDC Polar Stereographic South
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
    assert_eq!(
        (band.min_lon, band.max_lon, band.min_lat, band.max_lat),
        (-6.0, 0.0, 0.0, 84.0)
    );
    assert!(utm30n.contains(GeoPoint::from_lat_lon(51.5, -0.1))); // London, in band
    assert!(!utm30n.contains(GeoPoint::from_lat_lon(51.5, 50.0))); // far east, out
    assert!(!utm30n.contains(GeoPoint::from_lat_lon(-5.0, -3.0))); // south of equator, out

    // A southern MGA zone bands into the southern hemisphere (lat −80..0).
    let mga56 = Crs::from_epsg(28_356).unwrap().area_of_use().unwrap();
    assert_eq!((mga56.min_lat, mga56.max_lat), (-80.0, 0.0));

    // An antimeridian-wrapping geographic extent (NAD83: lon 167.65 → −40.73).
    let nad83 = Crs::from_epsg(4269).unwrap();
    let aou = nad83.area_of_use().unwrap();
    assert!(
        aou.min_lon > aou.max_lon,
        "NAD83 area of use wraps the antimeridian"
    );
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
    assert!(matches!(
        wgs84_3d,
        Crs::Geographic {
            dim: GeoDim::D3,
            ..
        }
    ));
    assert_eq!(wgs84_3d.epsg(), Some(4979));
    let etrs89_3d = Crs::from_epsg(4937).unwrap();
    assert!(matches!(
        etrs89_3d,
        Crs::Geographic {
            dim: GeoDim::D3,
            ..
        }
    ));
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
        let CrsCoord::Geographic {
            point: back,
            height: hb,
        } = transform(ecef, wgs84_3d, xyz)
        else {
            panic!()
        };
        assert!(approx(back.lat, p.lat, 1e-9) && approx(back.lon, p.lon, 1e-9));
        assert!(approx(hb, h, 1e-6), "height {h} → {hb}");
    }
    // A 2-D CRS (4326) drops the height: a 3-D input comes out at height 0.
    let CrsCoord::Geographic { height: h2d, .. } = transform(
        wgs84_3d,
        wgs84_2d,
        CrsCoord::geographic_with_height(p, 1000.0),
    ) else {
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
    let CrsCoord::Geographic {
        point: shifted,
        height: hs,
    } = transform(osgb36_3d, wgs84_3d, CrsCoord::geographic_with_height(p, h))
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
    let CrsCoord::Geographic {
        point: back,
        height: hb,
    } = transform(
        wgs84_3d,
        osgb36_3d,
        CrsCoord::geographic_with_height(shifted, hs),
    )
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
    let equiv = compound
        .as_crs()
        .expect("geographic + ellipsoidal collapses");
    assert!(matches!(
        equiv,
        Crs::Geographic {
            dim: GeoDim::D3,
            ..
        }
    ));
    assert_eq!(equiv.datum(), Datum::WGS84);
    // It flows through `transform` like any 3-D geographic CRS — height retained.
    let ecef = Crs::from_epsg(4978).unwrap();
    let p = GeoPoint::from_lat_lon(40.0, -3.0);
    let c = CrsCoord::geographic_with_height(p, 655.0);
    let CrsCoord::Geographic { height: hb, .. } = transform(ecef, equiv, transform(equiv, ecef, c))
    else {
        panic!()
    };
    assert!(approx(hb, 655.0, 1e-6));
    // A projected-horizontal compound doesn't collapse yet (deferred).
    let projected = CompoundCrs::new(
        Crs::from_epsg(3857).unwrap(),
        VerticalCrs::Ellipsoidal { epsg: None },
    );
    assert!(projected.as_crs().is_none());
}
