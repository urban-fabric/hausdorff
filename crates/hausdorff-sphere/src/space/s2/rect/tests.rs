use super::*;
use crate::space::s2::RegionCoverer;

fn gp(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

#[test]
fn from_degrees_contains_and_excludes() {
    let r = Rect::from_degrees(40.0, -5.0, 45.0, 5.0);
    assert!(r.is_valid() && !r.is_empty() && !r.is_full());
    assert!(r.contains_lat_lng(gp(42.0, 0.0)));
    assert!(r.contains_lat_lng(gp(40.0, -5.0)), "boundary is closed");
    assert!(!r.contains_lat_lng(gp(39.9, 0.0)));
    assert!(!r.contains_lat_lng(gp(42.0, 5.1)));
    assert!((r.lo().lat - 40.0).abs() < 1e-12 && (r.hi().lon - 5.0).abs() < 1e-12);
}

#[test]
fn empty_and_full() {
    assert!(Rect::empty().is_empty());
    assert!(Rect::empty().is_valid());
    assert!(!Rect::empty().contains_lat_lng(gp(0.0, 0.0)));
    assert!(Rect::full().is_full());
    for &(lat, lon) in &[(0.0, 0.0), (90.0, 0.0), (-90.0, 0.0), (15.0, 180.0)] {
        assert!(Rect::full().contains_lat_lng(gp(lat, lon)));
    }
    assert!(Rect::full().contains_rect(&Rect::from_degrees(1.0, 2.0, 3.0, 4.0)));
}

#[test]
fn antimeridian_wrap_contains_and_set_ops() {
    let r = Rect::from_degrees(10.0, 170.0, 20.0, -170.0);
    assert!(r.is_valid());
    assert!(r.contains_lat_lng(gp(15.0, 175.0)));
    assert!(r.contains_lat_lng(gp(15.0, -175.0)));
    assert!(r.contains_lat_lng(gp(15.0, 180.0)));
    assert!(r.contains_lat_lng(gp(15.0, -180.0)));
    assert!(!r.contains_lat_lng(gp(15.0, 0.0)));
    assert!(!r.contains_lat_lng(gp(15.0, -160.0)));
    assert!(!r.contains_lat_lng(gp(25.0, 175.0)));

    let other = Rect::from_degrees(12.0, -175.0, 18.0, -150.0);
    let inter = r.intersection(&other);
    assert!(inter.contains_lat_lng(gp(15.0, -172.0)));
    assert!(!inter.contains_lat_lng(gp(15.0, -160.0)));
    assert!(!inter.contains_lat_lng(gp(15.0, 175.0)));

    let union = r.union(&other);
    assert!(union.contains_lat_lng(gp(15.0, 175.0)));
    assert!(union.contains_lat_lng(gp(15.0, -160.0)));
    assert!(!union.contains_lat_lng(gp(15.0, 0.0)));
    assert!(union.contains_rect(&r) && union.contains_rect(&other));

    let disjoint = Rect::from_degrees(12.0, -10.0, 18.0, 10.0);
    assert!(r.intersection(&disjoint).is_empty());
}

#[test]
fn union_intersection_of_plain_rects() {
    let a = Rect::from_degrees(0.0, 0.0, 10.0, 10.0);
    let b = Rect::from_degrees(5.0, 5.0, 15.0, 15.0);
    let u = a.union(&b);
    let i = a.intersection(&b);
    assert!(u.contains_rect(&a) && u.contains_rect(&b) && u.contains_rect(&i));
    assert!(a.contains_rect(&i) && b.contains_rect(&i));
    assert!(i.contains_lat_lng(gp(7.0, 7.0)));
    assert!(!i.contains_lat_lng(gp(3.0, 3.0)));
    assert!(u.contains_lat_lng(gp(3.0, 3.0)) && u.contains_lat_lng(gp(13.0, 13.0)));
    assert_eq!(a.union(&Rect::empty()), a);
    assert!(a.intersection(&Rect::empty()).is_empty());
}

#[test]
fn cell_bound_contains_cell_samples() {
    for &(lat, lon, level) in &[
        (0.0, 0.0, 3),
        (37.7749, -122.4194, 10),
        (-33.8688, 151.2093, 6),
        (84.0, 5.0, 4),
        (-84.0, 5.0, 4),
        (0.1, 179.9, 8),
        (45.0, 45.0, 1),
    ] {
        let cell = S2CellId::from_geo_point_at_level(gp(lat, lon), level);
        let bound = Rect::cell_bound(cell);
        let mut samples: Vec<S2Point> = vec![cell.center_point()];
        let verts = cell.vertices();
        samples.extend_from_slice(&verts);
        for k in 0..4 {
            let mid = verts[k]
                .vector()
                .add(verts[(k + 1) % 4].vector())
                .normalize();
            samples.push(S2Point::from_vector(mid));
        }
        for s in samples {
            assert!(
                bound.contains_lat_lng(s.to_geo_point()),
                "cell {cell} sample escaped its bound"
            );
        }
    }
    for face in 0..6u8 {
        let cell = S2CellId::from_face(face);
        let bound = Rect::cell_bound(cell);
        assert!(bound.contains_lat_lng(cell.center()));
        for v in cell.vertices() {
            assert!(bound.contains_lat_lng(v.to_geo_point()));
        }
    }
}

#[test]
fn region_impl_covers_rect_and_respects_strict_contains() {
    let r = Rect::from_degrees(40.0, -5.0, 45.0, 5.0);
    let cover = RegionCoverer::at_level(7).covering(&r);
    assert!(!cover.is_empty());
    for &c in &cover {
        assert!(r.intersects_cell(c));
    }
    for &(lat, lon) in &[(40.5, -4.5), (42.5, 0.0), (44.9, 4.9)] {
        let leaf = S2CellId::from_geo_point(gp(lat, lon));
        assert!(
            cover.iter().any(|c| c.contains(leaf)),
            "interior point ({lat},{lon}) not covered"
        );
    }
    let outside = S2CellId::from_geo_point(gp(50.0, 20.0));
    assert!(!cover.iter().any(|c| c.contains(outside)));

    let deep_inside = S2CellId::from_geo_point_at_level(gp(42.5, 0.0), 10);
    assert!(r.contains_cell(deep_inside));
    assert!(r.intersects_cell(deep_inside));
    let far = S2CellId::from_geo_point_at_level(gp(-40.0, 120.0), 10);
    assert!(!r.intersects_cell(far));
    assert!(!r.contains_cell(far));
}

#[test]
fn antimeridian_rect_adaptive_cover() {
    let r = Rect::from_degrees(10.0, 170.0, 20.0, -170.0);
    let cover = RegionCoverer {
        min_level: 0,
        max_level: 12,
        max_cells: 32,
    }
    .adaptive_covering(&r);
    assert!(!cover.is_empty());
    assert!(cover.len() <= 32);
    for &(lat, lon) in &[(15.0, 178.0), (15.0, -178.0), (11.0, 171.0), (19.0, -171.0)] {
        assert!(
            cover.contains_point(S2Point::from_geo_point(gp(lat, lon))),
            "wrapped point ({lat},{lon}) not covered"
        );
    }
}

#[test]
fn bounding_cap_contains_rect_samples() {
    for r in [
        Rect::from_degrees(40.0, -5.0, 45.0, 5.0),
        Rect::from_degrees(10.0, 170.0, 20.0, -170.0),
        Rect::from_degrees(70.0, -180.0, 90.0, 180.0),
        Rect::from_degrees(-90.0, -180.0, -70.0, 180.0),
    ] {
        let cap = r.bounding_cap();
        for p in r.corner_points() {
            assert!(cap.contains_point(p), "corner escaped the bounding cap");
        }
        let mid =
            GeoPoint::from_lat_lon(0.5 * (r.lo().lat + r.hi().lat), r.lng.center().to_degrees());
        assert!(cap.contains_point(S2Point::from_geo_point(mid)));
    }
    assert!(Rect::empty().bounding_cap().is_empty());
}

#[test]
fn from_geo_points_bounds_all_vertices() {
    let pts = [gp(40.0, -5.0), gp(45.0, 5.0), gp(42.0, 0.0), gp(41.0, -3.0)];
    let bound = Rect::from_geo_points(&pts);
    assert!(!bound.is_empty() && bound.is_valid());
    for &p in &pts {
        assert!(bound.contains_lat_lng(p), "vertex {p:?} escaped its bound");
    }
    // A tight, non-wrapping box: lat [40,45], lng [-5,5].
    assert!(!bound.is_lng_inverted());
    assert!((bound.lo().lat - 40.0).abs() < 1e-9 && (bound.hi().lat - 45.0).abs() < 1e-9);
    assert!((bound.lo().lon - -5.0).abs() < 1e-9 && (bound.hi().lon - 5.0).abs() < 1e-9);
    assert_eq!(Rect::from_geo_points(&[]), Rect::empty());
}

#[test]
fn from_geo_points_bounds_antimeridian_ring_tightly() {
    // A ring straddling ±180° (a 20°-wide box centred on the antimeridian). The
    // vertex longitudes are 175, -175 — a min/max box would span 350° of the
    // globe; the wrap-aware bound spans only 10° the short way.
    let ring = [
        gp(10.0, 175.0),
        gp(10.0, -175.0),
        gp(20.0, -175.0),
        gp(20.0, 175.0),
    ];
    let bound = Rect::from_geo_points(&ring);
    assert!(bound.is_valid());
    assert!(bound.is_lng_inverted(), "antimeridian ring must wrap");
    for &p in &ring {
        assert!(bound.contains_lat_lng(p));
    }
    // Points inside the short arc are covered; the far side of the globe is not.
    assert!(bound.contains_lat_lng(gp(15.0, 180.0)));
    assert!(bound.contains_lat_lng(gp(15.0, -180.0)));
    assert!(!bound.contains_lat_lng(gp(15.0, 0.0)));
    assert!(!bound.contains_lat_lng(gp(15.0, 90.0)));
    // Its longitude corners are the ring's east/west edges, inverted (lo > hi).
    assert!((bound.lo().lon - 175.0).abs() < 1e-9 && (bound.hi().lon - -175.0).abs() < 1e-9);
}

#[test]
fn add_point_grows_the_short_way() {
    let base = Rect::from_point(gp(0.0, 179.0));
    assert!(!base.is_lng_inverted());
    let grown = base.add_point(gp(0.0, -179.0));
    // 179 → -179 is 2° across the antimeridian, not 358° around.
    assert!(grown.is_lng_inverted());
    assert!(grown.contains_lat_lng(gp(0.0, 180.0)));
    assert!(!grown.contains_lat_lng(gp(0.0, 0.0)));
}
