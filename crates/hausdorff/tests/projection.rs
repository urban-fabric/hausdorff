//! Map-projection interchange: WGS84 geographic ↔ projected planar (Web
//! Mercator, UTM). Validated by analytic anchors, round-trips, and monotonicity.

use std::f64::consts::PI;

use hausdorff::earth::GeoPoint;
use hausdorff::earth::projection::{Ellipsoid, Projection, TransverseMercator, Utm, WebMercator};
use hausdorff::space::r2::{Area, Geometry, Length, LineString, Point, Polygon};

fn g(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

#[test]
fn web_mercator_anchors_and_round_trip() {
    let wm = WebMercator::EPSG_3857;
    let o = wm.project(g(0.0, 0.0));
    assert!(o.x.abs() < 1e-6 && o.y.abs() < 1e-6);
    let e = wm.project(g(0.0, 180.0));
    assert!((e.x - wm.radius * PI).abs() < 1e-3, "x = {}", e.x);
    for &(lat, lon) in &[(0.0, 0.0), (51.5, -0.12), (-33.87, 151.21), (60.0, 24.9)] {
        let back = wm.inverse(wm.project(g(lat, lon)));
        assert!((back.lat - lat).abs() < 1e-9 && (back.lon - lon).abs() < 1e-9);
    }
}

#[test]
fn try_utm_validates_zone() {
    assert!(TransverseMercator::try_utm(1, true).is_some());
    assert!(TransverseMercator::try_utm(60, false).is_some());
    assert!(TransverseMercator::try_utm(0, true).is_none());
    assert!(TransverseMercator::try_utm(61, true).is_none());
    assert_eq!(
        TransverseMercator::try_utm(33, true),
        Some(TransverseMercator::utm(33, true))
    );
}

#[test]
fn utm_central_meridian_anchor() {
    let u = Utm::from_geo(g(0.0, 3.0));
    assert_eq!(u.zone, 31);
    assert!(u.north);
    assert!(
        (u.easting - 500_000.0).abs() < 1e-3,
        "easting {}",
        u.easting
    );
    assert!(u.northing.abs() < 1e-3, "northing {}", u.northing);
}

#[test]
fn utm_round_trips_to_mm() {
    for &(lat, lon) in &[
        (40.7128, -74.0060),
        (48.8566, 2.3522),
        (-33.8688, 151.2093),
        (35.6762, 139.6503),
    ] {
        let back = Utm::from_geo(g(lat, lon)).to_geo();
        assert!((back.lat - lat).abs() < 1e-8, "lat {lat} -> {}", back.lat);
        assert!((back.lon - lon).abs() < 1e-8, "lon {lon} -> {}", back.lon);
    }
}

#[test]
fn utm_easting_increases_eastward() {
    let cm = Utm::from_geo(g(45.0, 3.0)).easting;
    let east = Utm::from_geo(g(45.0, 5.0)).easting;
    assert!(east > cm && (cm - 500_000.0).abs() < 1.0);
}

#[test]
fn geometry_projects_and_round_trips() {
    // A city-block polygon (with a hole) in lon/lat degrees around Paris.
    let poly = Polygon::new(
        LineString::from(vec![
            (2.350, 48.853),
            (2.360, 48.853),
            (2.360, 48.860),
            (2.350, 48.860),
        ]),
        vec![LineString::from(vec![
            (2.354, 48.855),
            (2.356, 48.855),
            (2.356, 48.857),
            (2.354, 48.857),
        ])],
    );
    let g = Geometry::Polygon(poly);
    let utm = TransverseMercator::utm(31, true);

    let projected = utm.project_geometry(&g);
    // In projected metres the planar area must be physically plausible:
    // ~0.01° lon × 0.007° lat at 48.86°N ≈ (733 m × 778 m) minus the hole.
    let Geometry::Polygon(p) = &projected else {
        panic!("polygon should project to polygon");
    };
    let area = p.unsigned_area();
    assert!(
        (500_000.0..600_000.0).contains(&area),
        "projected area {area} m²"
    );
    assert!(p.exterior().is_closed() && p.interiors().len() == 1);

    let back = utm.inverse_geometry(&projected);
    let Geometry::Polygon(b) = &back else {
        panic!("round trip should preserve the variant");
    };
    for (orig, rt) in g_coords(&g).zip(g_coords(&Geometry::Polygon(b.clone()))) {
        assert!((orig.0 - rt.0).abs() < 1e-8 && (orig.1 - rt.1).abs() < 1e-8);
    }
}

fn g_coords(g: &Geometry) -> impl Iterator<Item = (f64, f64)> + '_ {
    let Geometry::Polygon(p) = g else {
        panic!("test helper expects a polygon")
    };
    p.exterior()
        .coords()
        .chain(p.interiors().iter().flat_map(|r| r.coords()))
        .map(|c| (c.x, c.y))
}

#[test]
fn projected_length_approximates_haversine() {
    // A ~1.1 km segment: its planar length in UTM metres should approximate
    // the great-circle metres. The dominant error is haversine's own spherical
    // model (mean-radius sphere vs the ellipsoid, up to ~0.5%), not the
    // projection (UTM scale error ≤ 4e-4 near the CM).
    let ls = LineString::from(vec![(2.350, 48.853), (2.360, 48.860)]);
    let utm = TransverseMercator::utm(31, true);
    let projected = utm.project_geometry(&ls);

    let planar = hausdorff::space::r2::Euclidean.length(&projected);
    let spherical = hausdorff::space::r2::Haversine.length(&ls);
    assert!(
        (planar - spherical).abs() / spherical < 5e-3,
        "planar {planar} vs spherical {spherical}"
    );
}

#[test]
fn project_geometry_works_for_bare_types() {
    // The generic surface accepts any MapCoords geometry, not just the enum.
    let utm = TransverseMercator::utm(31, true);
    let pt = Point::new(3.0, 0.0);
    let projected = utm.project_geometry(&pt);
    assert!((projected.x() - 500_000.0).abs() < 1e-3);
    let back = utm.inverse_geometry(&projected);
    assert!((back.x() - 3.0).abs() < 1e-8 && back.y().abs() < 1e-8);
}

#[test]
fn projection_is_ellipsoid_flexible() {
    let wgs84 = TransverseMercator {
        ellipsoid: Ellipsoid::WGS84,
        central_meridian: 0.0,
        lat_origin: 0.0,
        k0: 0.9996,
        false_easting: 500_000.0,
        false_northing: 0.0,
    };
    let sphere = TransverseMercator {
        ellipsoid: Ellipsoid::SPHERE,
        ..wgs84
    };
    let p = g(45.0, 2.0);
    assert!((wgs84.project(p).y - sphere.project(p).y).abs() > 1.0);
    let back = sphere.inverse(sphere.project(p));
    assert!((back.lat - 45.0).abs() < 1e-8 && (back.lon - 2.0).abs() < 1e-8);
}
