//! Tests for the SoA batch geometry types. The batch kernels must produce
//! exactly the same results as the per-element scalar APIs — whether built with
//! the `parallel` feature or not (run the suite both ways to cover the seam).

use crate::earth::{BBox, GeoPoint, LengthMeters, LonLatBatch, haversine_m};
use hausdorff_index::index::h3::H3Cell;
use hausdorff_planar::space::r2::{Coord, PointBatch};

fn sample_points() -> Vec<GeoPoint> {
    vec![
        GeoPoint::from_lat_lon(52.52, 13.405),
        GeoPoint::from_lat_lon(48.8566, 2.3522),
        GeoPoint::from_lat_lon(-33.8688, 151.2093),
        GeoPoint::from_lat_lon(0.0, 0.0),
        GeoPoint::from_lat_lon(40.7128, -74.006),
    ]
}

#[test]
fn lonlat_batch_construction_and_columns() {
    let pts = sample_points();
    let b = LonLatBatch::from_points(&pts);
    assert_eq!(b.len(), pts.len());
    assert!(!b.is_empty());
    for (i, p) in pts.iter().enumerate() {
        assert_eq!(b.point(i), *p);
        assert_eq!(b.lons()[i], p.lon);
        assert_eq!(b.lats()[i], p.lat);
    }
    assert!(LonLatBatch::from_columns(vec![1.0, 2.0], vec![3.0, 4.0]).is_some());
    assert!(LonLatBatch::from_columns(vec![1.0], vec![3.0, 4.0]).is_none());
    assert!(LonLatBatch::default().is_empty());
}

#[test]
fn lonlat_batch_haversine_matches_scalar() {
    let pts = sample_points();
    let b = LonLatBatch::from_points(&pts);
    let target = GeoPoint::from_lat_lon(51.5074, -0.1278);
    let fan = b.haversine_to(target);
    let raw = b.haversine_to_m(target);
    assert_eq!(fan.len(), pts.len());
    assert_eq!(raw.len(), pts.len());
    for (i, p) in pts.iter().enumerate() {
        let expected = haversine_m(*p, target);
        assert_eq!(fan[i].value(), expected);
        assert_eq!(raw[i], expected);
    }
}

#[test]
fn lonlat_batch_masks_match_scalar_predicates() {
    let pts = sample_points();
    let b = LonLatBatch::from_points(&pts);
    let bbox = BBox::new(-40.0, -10.0, 55.0, 20.0);
    let target = GeoPoint::from_lat_lon(51.5074, -0.1278);

    assert_eq!(
        b.bbox_contains_mask(bbox),
        pts.iter().map(|&p| bbox.contains(p)).collect::<Vec<_>>()
    );

    let radius_m = 1_000_000.0;
    assert_eq!(
        b.within_radius_m(target, radius_m),
        pts.iter()
            .map(|&p| p.is_finite() && target.is_finite() && haversine_m(p, target) <= radius_m)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        b.within_radius(target, LengthMeters::new(radius_m)),
        b.within_radius_latlon_m(target.lat, target.lon, radius_m)
    );

    let bad = LonLatBatch::from_columns(vec![0.0, 10.0, f64::NAN], vec![0.0, 91.0, 10.0])
        .expect("same length");
    assert_eq!(bad.valid_lonlat_mask(), vec![true, false, false]);
}

#[test]
fn lonlat_batch_h3_matches_scalar() {
    let pts = sample_points();
    let b = LonLatBatch::from_points(&pts);
    for res in [0u8, 5, 9, 15] {
        let cells = b.h3_cells(res);
        for (i, p) in pts.iter().enumerate() {
            assert_eq!(
                cells[i],
                H3Cell::for_point(p.lat, p.lon, res),
                "res {res} point {i}"
            );
        }
    }
    let bad = LonLatBatch::from_columns(vec![0.0, 200.0], vec![91.0, 0.0]).unwrap();
    let cells = bad.h3_cells(7);
    assert!(cells[0].is_none() && cells[1].is_none());
}

#[test]
fn lonlat_batch_bounding_box_matches_scalar() {
    let pts = sample_points();
    let b = LonLatBatch::from_points(&pts);
    assert_eq!(b.bounding_box(), Some(BBox::from_coords(&pts)));
    assert_eq!(LonLatBatch::default().bounding_box(), None);
}

#[test]
fn point_batch_matches_scalar() {
    let coords = [
        Coord::new(1.0, 2.0),
        Coord::new(-3.0, 5.0),
        Coord::new(4.0, -1.0),
        Coord::new(0.0, 0.0),
    ];
    let b = PointBatch::from_coords(&coords);
    assert_eq!(b.len(), coords.len());
    for (i, c) in coords.iter().enumerate() {
        assert_eq!(b.coord(i), *c);
        assert_eq!(b.xs()[i], c.x);
        assert_eq!(b.ys()[i], c.y);
    }
    let r = b.bounding_rect().expect("non-empty");
    assert_eq!(r.min(), Coord::new(-3.0, -1.0));
    assert_eq!(r.max(), Coord::new(4.0, 5.0));
    assert_eq!(PointBatch::default().bounding_rect(), None);
    assert!(PointBatch::from_columns(vec![1.0], vec![]).is_none());
}

#[test]
fn point_batch_distances_to_matches_scalar() {
    let coords: Vec<Coord> = (0..13)
        .map(|i| Coord::new(f64::from(i) * 1.5 - 4.0, f64::from(i % 5) * -2.25 + 1.0))
        .collect();
    let b = PointBatch::from_coords(&coords);
    let target = Coord::new(0.5, -0.5);

    let got = b.distances_to(target);
    assert_eq!(got.len(), coords.len());
    for (i, c) in coords.iter().enumerate() {
        let (dx, dy) = (c.x - target.x, c.y - target.y);
        let want = (dx * dx + dy * dy).sqrt();
        assert_eq!(got[i], want, "distance mismatch at {i}");
    }
    assert!(PointBatch::default().distances_to(target).is_empty());
}
