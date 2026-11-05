//! Public-API tests for the S2 region layer: `CellUnion` set algebra, the
//! adaptive `RegionCoverer`, the lat/lng `Rect` region, and `Polyline` line
//! covers — invariant/property checks, no external oracle.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use crate::space::s2::{
    CellUnion, Polyline, Rect, Region, RegionCoverer, S2CellId, S2Point, cover_r2_geometry,
    cover_r2_geometry_at_level,
};
use hausdorff_base::earth::GeoPoint;

fn gp(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

fn pt(lat: f64, lon: f64) -> S2Point {
    S2Point::from_geo_point(gp(lat, lon))
}

#[test]
fn cell_union_normalizes_and_answers_queries() {
    let cell = S2CellId::from_geo_point_at_level(gp(48.85, 2.35), 9);
    let mut ids = cell.children().to_vec();
    ids.push(cell.children()[0].children()[1]);
    let u = CellUnion::from(ids);
    assert_eq!(u.cell_ids(), &[cell][..]);
    assert!(u.contains_point(pt(48.85, 2.35)));
    assert!(!u.contains_point(pt(-20.0, 100.0)));
    assert!(u.contains_cell(cell.children()[2]));
    assert!(u.intersects_cell(cell.parent(3)));
    assert!(u.area_m2() > 0.0);
}

#[test]
fn cell_union_set_algebra_round_trips() {
    let a = CellUnion::from(vec![S2CellId::from_geo_point_at_level(gp(10.0, 10.0), 6)]);
    let b = CellUnion::from(vec![
        S2CellId::from_geo_point_at_level(gp(10.0, 10.0), 8),
        S2CellId::from_geo_point_at_level(gp(40.0, -100.0), 7),
    ]);
    let union = a.union(&b);
    let inter = a.intersection(&b);
    let diff = a.difference(&b);
    assert!(union.contains_cell(a.cell_ids()[0]));
    for &c in inter.cell_ids() {
        assert!(a.contains_cell(c) && b.contains_cell(c));
    }
    assert!(diff.intersection(&b).is_empty());
    assert_eq!(diff.union(&inter), a, "A = (A∖B) ∪ (A∩B) when B ⊆ ranges");
}

#[test]
fn adaptive_covering_of_rect_covers_and_respects_budget() {
    let rect = Rect::from_degrees(10.0, 170.0, 20.0, -170.0);
    let coverer = RegionCoverer {
        min_level: 0,
        max_level: 14,
        max_cells: 20,
    };
    let cover = coverer.adaptive_covering(&rect);
    assert!(!cover.is_empty());
    assert!(cover.len() <= 20);
    for c in cover.iter() {
        assert!(rect.intersects_cell(c));
    }
    for &(lat, lon) in &[(15.0, 179.5), (15.0, -179.5), (10.5, 171.0), (19.5, -171.0)] {
        assert!(
            cover.contains_point(pt(lat, lon)),
            "({lat},{lon}) escaped the antimeridian cover"
        );
    }
    assert!(!cover.contains_point(pt(15.0, 0.0)) || rect.contains_lat_lng(gp(15.0, 0.0)));

    let interior = coverer.interior_covering(&rect);
    for c in interior.iter() {
        assert!(rect.contains_cell(c));
        assert!(cover.contains_cell(c));
    }
}

/// Unit vector for a WGS84 coordinate (the test-local version of the s2
/// conversion, so edge samples follow the true great circle).
fn xyz(p: GeoPoint) -> [f64; 3] {
    let (lat, lon) = (p.lat.to_radians(), p.lon.to_radians());
    [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()]
}

/// Spherical interpolation between two coordinates at fraction `t`.
fn slerp(a: GeoPoint, b: GeoPoint, t: f64) -> GeoPoint {
    let (va, vb) = (xyz(a), xyz(b));
    let dot: f64 = va.iter().zip(&vb).map(|(x, y)| x * y).sum();
    let theta = dot.clamp(-1.0, 1.0).acos();
    let (wa, wb) = (((1.0 - t) * theta).sin(), (t * theta).sin());
    let v: Vec<f64> = va.iter().zip(&vb).map(|(x, y)| wa * x + wb * y).collect();
    let lat = v[2].atan2(v[0].hypot(v[1]));
    let lon = v[1].atan2(v[0]);
    gp(lat.to_degrees(), lon.to_degrees())
}

#[test]
fn polyline_cover_is_boundary_complete() {
    let waypoints = [gp(35.0, 139.0), gp(35.5, 139.8), gp(36.0, 139.4)];
    let line = Polyline::from_geo_points(&waypoints);
    let level = 10;
    let cover = CellUnion::from(RegionCoverer::at_level(level).covering(&line));
    assert!(!cover.is_empty());
    for w in waypoints.windows(2) {
        for k in 0..=200 {
            let t = f64::from(k) / 200.0;
            let p = S2Point::from_geo_point(slerp(w[0], w[1], t));
            assert!(cover.contains_point(p), "sampled line point not covered");
        }
    }
}

#[test]
fn r2_line_geometries_cover_through_the_public_entry_points() {
    use hausdorff_planar::space::r2::{Geometry, LineString};
    let line = Geometry::LineString(LineString::from(vec![(10.0, 50.0), (10.5, 50.4)]));
    let fixed = cover_r2_geometry_at_level(&line, 11).unwrap();
    assert!(!fixed.is_empty());
    let start = S2CellId::from_geo_point_at_level(gp(50.0, 10.0), 11);
    assert!(fixed.contains(&start));

    let coverer = RegionCoverer {
        min_level: 0,
        max_level: 12,
        max_cells: 16,
    };
    let adaptive = cover_r2_geometry(&line, &coverer).unwrap();
    assert!(!adaptive.is_empty());
    assert!(adaptive.len() <= 16);
    assert!(adaptive.contains_point(pt(50.2, 10.25)));
}
