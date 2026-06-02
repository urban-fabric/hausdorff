//! Public-API tests for `hausdorff::convert`, relocated out of `src/`.

use hausdorff::convert::*;
use hausdorff::earth::GeoPoint;
use hausdorff::space::r2::Coord;

#[test]
fn geo_point_round_trips_through_coord() {
    let p = GeoPoint::from_lat_lon(40.0, -74.0);
    let c = geo_point_to_coord(p);
    assert_eq!(c, Coord::new(-74.0, 40.0)); // x = lon, y = lat
    let back = coord_to_geo_point(c);
    assert_eq!(back, p);
}

#[test]
fn line_string_round_trips() {
    let pts = [
        GeoPoint::from_lat_lon(0.0, 0.0),
        GeoPoint::from_lat_lon(1.0, 2.0),
    ];
    let ls = geo_points_to_line_string(&pts);
    assert_eq!(ls.0[1], Coord::new(2.0, 1.0));
    assert_eq!(line_string_to_geo_points(&ls), pts.to_vec());
}
