//! Façade-surface guard: the workspace split must keep every original
//! `hausdorff::…` path working. Compile-time coverage for the paths the rest
//! of the suite doesn't already exercise — the transitional `geo` alias, the
//! root re-exports, the construction macros, and the assembled `earth::batch`
//! / `config` modules.

use hausdorff::geo::GeoPoint as AliasGeoPoint;
use hausdorff::{GeoPoint, MetricSpace as _, coord, point};

#[test]
fn transitional_geo_alias_still_resolves() {
    let p: AliasGeoPoint = GeoPoint::from_lat_lon(52.52, 13.405);
    assert_eq!(p.lat, 52.52);
    assert!(hausdorff::geo::haversine_m(p, p).abs() < 1e-12);
}

#[test]
fn root_reexports_and_macros() {
    let d = hausdorff::Euclidean.distance(coord! { x: 0.0, y: 0.0 }, coord! { x: 3.0, y: 4.0 });
    assert_eq!(d, 5.0);
    assert_eq!(point! { x: 1.0, y: 2.0 }, hausdorff::Point::new(1.0, 2.0));
}

#[test]
fn assembled_modules_reachable() {
    let p = GeoPoint::from_lat_lon(52.52, 13.405);
    let batch = hausdorff::earth::batch::LonLatBatch::from_points(std::slice::from_ref(&p));
    assert_eq!(batch.len(), 1);
    let cfg = hausdorff::config::config();
    assert!(cfg.index.node_capacity > 0);
    let sphere = hausdorff::kernels::earth::Sphere::<f64>::mean_earth();
    assert!(sphere.radius_m > 0.0);
}
