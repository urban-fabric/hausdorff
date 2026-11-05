use std::f64::consts::PI;

use super::*;
use hausdorff_base::earth::GeoPoint;

fn pt(lon: f64, lat: f64) -> S2Point {
    S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
}

fn leaf(lon: f64, lat: f64) -> S2CellId {
    S2CellId::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
}

/// A small CCW lon/lat box `[lon0,lon1] × [lat0,lat1]` as a spherical loop
/// (interior to the left — the small side).
fn square(lon0: f64, lat0: f64, lon1: f64, lat1: f64) -> Loop {
    let l = Loop::new(vec![
        pt(lon0, lat0),
        pt(lon1, lat0),
        pt(lon1, lat1),
        pt(lon0, lat1),
    ]);
    if l.area() > 2.0 * PI {
        let mut v = l.vertices().to_vec();
        v.reverse();
        Loop::new(v)
    } else {
        l
    }
}

#[test]
fn cell_loop_is_oriented_and_contains_its_center() {
    for &(lon, lat) in &[(0.0, 0.0), (90.0, 10.0), (-120.0, -45.0), (179.0, 80.0)] {
        let cell = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(lat, lon), 6);
        let cl = cell_loop(cell);
        assert!(cl.contains_point(cell.center_point()), "({lon},{lat})");
    }
}

#[test]
fn covering_contains_interior_excludes_exterior() {
    let region = square(-1.0, -1.0, 1.0, 1.0);
    assert!(
        region.contains_point(pt(0.0, 0.0)),
        "origin inside the loop"
    );
    assert!(!region.contains_point(pt(20.0, 20.0)), "far point outside");

    let cover = RegionCoverer::at_level(8).covering(&region);
    assert!(!cover.is_empty());
    for &c in &cover {
        assert_eq!(c.level(), 8);
        assert!(region.intersects_cell(c));
    }
    assert!(
        cover.iter().any(|c| c.contains(leaf(0.0, 0.0))),
        "interior point not covered"
    );
    assert!(
        !cover.iter().any(|c| c.contains(leaf(45.0, 45.0))),
        "exterior point wrongly covered"
    );
}

#[test]
fn multi_level_cover_bounds_levels_and_covers_interior() {
    let region = square(-2.0, -2.0, 2.0, 2.0);
    let cover = RegionCoverer {
        min_level: 2,
        max_level: 12,
        max_cells: 64,
    }
    .covering(&region);
    assert!(!cover.is_empty());
    for &c in &cover {
        assert!((2..=12).contains(&c.level()));
        assert!(region.intersects_cell(c));
    }
    assert!(cover.iter().any(|c| c.contains(leaf(0.0, 0.0))));
}

#[test]
fn polygon_hole_is_excluded_from_the_cover() {
    let poly = Polygon::new(
        square(-4.0, -4.0, 4.0, 4.0),
        vec![square(-1.0, -1.0, 1.0, 1.0)],
    );
    let cover = RegionCoverer::at_level(10).covering(&poly);
    assert!(
        cover.iter().any(|c| c.contains(leaf(0.0, 3.0))),
        "ring point should be covered"
    );
    assert!(
        !cover.iter().any(|c| c.contains(leaf(0.0, 0.0))),
        "deep-hole point should not be covered"
    );
}

#[test]
fn at_level_clamps_to_the_leaf() {
    assert_eq!(RegionCoverer::at_level(7).max_level, 7);
    let clamped = RegionCoverer::at_level(40);
    assert_eq!(clamped.min_level, MAX_LEVEL);
    assert_eq!(clamped.max_level, MAX_LEVEL);
}

/// Sample points on a lon/lat grid over `[lon0,lon1] × [lat0,lat1]`,
/// keeping those the region contains.
fn interior_samples<R: Region>(
    region: &R,
    lon0: f64,
    lat0: f64,
    lon1: f64,
    lat1: f64,
) -> Vec<S2Point> {
    let mut samples = Vec::new();
    let steps = 20;
    for i in 0..=steps {
        for j in 0..=steps {
            let lon = lon0 + (lon1 - lon0) * f64::from(i) / f64::from(steps);
            let lat = lat0 + (lat1 - lat0) * f64::from(j) / f64::from(steps);
            let p = pt(lon, lat);
            if region.contains_point(p) {
                samples.push(p);
            }
        }
    }
    samples
}

#[test]
fn adaptive_covering_of_cap_respects_budget_and_covers_samples() {
    let cap = crate::space::s2::Cap::from_center_angle(
        pt(-100.0, 40.0),
        hausdorff_base::space::s1::Angle::from_radians(0.1),
    );
    let coverer = RegionCoverer {
        min_level: 0,
        max_level: 30,
        max_cells: 16,
    };
    let cover = coverer.adaptive_covering(&cap);
    assert!(!cover.is_empty());
    assert!(cover.len() <= 16, "budget exceeded: {}", cover.len());
    for c in cover.iter() {
        assert!(c.level() <= 30);
        assert!(cap.intersects_cell(c), "cover cell {c} misses the cap");
    }
    let samples = interior_samples(&cap, -110.0, 30.0, -90.0, 50.0);
    assert!(samples.len() > 50, "degenerate sample set");
    for p in samples {
        assert!(cover.contains_point(p), "interior point not covered");
    }
}

#[test]
fn adaptive_covering_of_polygon_is_multi_level_and_complete() {
    let poly = Polygon::new(
        square(-4.0, -4.0, 4.0, 4.0),
        vec![square(-1.0, -1.0, 1.0, 1.0)],
    );
    let coverer = RegionCoverer {
        min_level: 2,
        max_level: 14,
        max_cells: 64,
    };
    let cover = coverer.adaptive_covering(&poly);
    assert!(!cover.is_empty());
    assert!(cover.len() <= 64);
    let mut levels = std::collections::BTreeSet::new();
    for c in cover.iter() {
        assert!(c.level() <= 14);
        assert!(poly.intersects_cell(c));
        levels.insert(c.level());
    }
    assert!(levels.len() > 1, "expected a multi-level cover: {levels:?}");
    for p in interior_samples(&poly, -5.0, -5.0, 5.0, 5.0) {
        assert!(cover.contains_point(p), "interior point not covered");
    }
}

#[test]
fn interior_covering_cells_are_contained_in_the_region() {
    let cap = crate::space::s2::Cap::from_center_angle(
        pt(10.0, 10.0),
        hausdorff_base::space::s1::Angle::from_radians(0.15),
    );
    let coverer = RegionCoverer {
        min_level: 0,
        max_level: 12,
        max_cells: 32,
    };
    let interior = coverer.interior_covering(&cap);
    assert!(!interior.is_empty(), "fat cap should admit interior cells");
    assert!(interior.len() <= 32);
    for c in interior.iter() {
        assert!(
            cap.contains_cell(c),
            "interior cell {c} pokes outside the cap"
        );
    }
    let cover = coverer.adaptive_covering(&cap);
    for c in interior.iter() {
        assert!(
            cover.contains_cell(c),
            "interior covering must be inside the covering"
        );
    }
}

#[test]
fn adaptive_covering_with_tight_budget_stays_coarse() {
    let region = square(-30.0, -30.0, 30.0, 30.0);
    let small = RegionCoverer {
        min_level: 0,
        max_level: 30,
        max_cells: 6,
    }
    .adaptive_covering(&region);
    assert!(!small.is_empty());
    assert!(small.len() <= 6, "budget exceeded: {}", small.len());
    let large = RegionCoverer {
        min_level: 0,
        max_level: 30,
        max_cells: 60,
    }
    .adaptive_covering(&region);
    assert!(
        large.len() > small.len(),
        "bigger budget should refine more"
    );
    for p in interior_samples(&region, -29.0, -29.0, 29.0, 29.0) {
        assert!(small.contains_point(p));
        assert!(large.contains_point(p));
    }
}

#[test]
fn boundary_coincident_neighbour_is_not_dropped() {
    let c = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(12.0, 34.0), 10);
    let region = cell_loop(c);
    for n in c.edge_neighbors() {
        assert!(
            region.intersects_cell(n),
            "edge-neighbour sharing a boundary with the region was dropped"
        );
    }
}
