//! Public-API tests for the self-contained WGS84 earth core.

mod point {
    use crate::earth::*;

    #[test]
    fn round_trip_within_quantization_tolerance() {
        for &(lat, lon) in &[
            (40.0, -74.0),
            (40.123_456_7, -74.987_654_3),
            (-89.999_999_9, 179.999_999_9),
            (0.0, 0.0),
        ] {
            let p = GeoPoint::from_lat_lon(lat, lon);
            let back = FixedGeoPoint::from(p).to_geo();
            assert!((back.lat - lat).abs() <= 1e-7, "lat {lat} -> {}", back.lat);
            assert!((back.lon - lon).abs() <= 1e-7, "lon {lon} -> {}", back.lon);
        }
    }
}

mod bbox {
    use crate::earth::*;

    #[test]
    fn bbox_contains_and_parses() {
        let bbox: BBox = "40.0,-74.0,41.0,-73.0".parse().unwrap();
        assert!(bbox.contains(GeoPoint::from_lat_lon(40.5, -73.5)));
        assert!(!bbox.contains(GeoPoint::from_lat_lon(0.0, 0.0)));
        assert!(!bbox.contains(GeoPoint::from_lat_lon(f64::NAN, -73.5)));
        assert_eq!(
            "40.0,-74.0".parse::<BBox>(),
            Err(BBoxParseError::WrongArity(2))
        );
    }
}

mod bearing {
    use crate::earth::*;

    #[test]
    fn bearing_normalizes_into_unit_turn() {
        assert!((Bearing::from_degrees(450.0).degrees() - 90.0).abs() < 1e-9);
        assert!((Bearing::from_degrees(-90.0).degrees() - 270.0).abs() < 1e-9);
    }

    #[test]
    fn initial_bearing_points_east_and_north() {
        let origin = GeoPoint::from_lat_lon(0.0, 0.0);
        let east = Bearing::initial(origin, GeoPoint::from_lat_lon(0.0, 1.0));
        let north = Bearing::initial(origin, GeoPoint::from_lat_lon(1.0, 0.0));
        assert!(
            (east.degrees() - 90.0).abs() < 1e-6,
            "east: {}",
            east.degrees()
        );
        assert!(north.degrees() < 1e-6 || (360.0 - north.degrees()) < 1e-6);
    }
}

mod geodesic {
    use crate::earth::*;

    #[test]
    fn destination_is_inverse_of_haversine() {
        let berlin = GeoPoint::from_lat_lon(52.52, 13.405);
        for deg in [0.0, 37.0, 90.0, 180.0, 271.5, 359.0] {
            let dest = destination_m(berlin, Bearing::from_degrees(deg), 1_500.0);
            assert!(
                (haversine_m(berlin, dest) - 1_500.0).abs() < 0.5,
                "deg {deg}"
            );
        }
    }
}

mod radial_grid {
    use crate::earth::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn counts_and_dimensions() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(52.52, 13.405),
            15,
            LengthMeters::new(100.0),
            32,
        );
        assert_eq!(g.ring_count(), 15);
        assert_eq!(g.points_per_ring(), 32);
        assert_eq!(g.len(), 480);
        assert!(!g.is_empty());
        assert_eq!(g.to_points().len(), 480);
        assert_eq!(g.points().count(), 480);
        assert_eq!(g.rings().count(), 15);
    }

    #[test]
    fn radius_and_bearing_schedule() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            3,
            LengthMeters::new(100.0),
            4,
        );
        assert!(approx(g.radius(0).value(), 100.0, 1e-9));
        assert!(approx(g.radius(2).value(), 300.0, 1e-9));
        assert!(approx(g.bearing(0).degrees(), 0.0, 1e-9));
        assert!(approx(g.bearing(1).degrees(), 90.0, 1e-9));
        assert!(approx(g.bearing(2).degrees(), 180.0, 1e-9));
        assert!(approx(g.bearing(3).degrees(), 270.0, 1e-9));
    }

    #[test]
    fn zero_points_per_ring_is_empty_without_divide_by_zero() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(52.52, 13.405),
            3,
            LengthMeters::new(100.0),
            0,
        );
        assert!(g.is_empty());
        assert_eq!(g.len(), 0);
        assert_eq!(g.points().count(), 0);
        assert!(g.bearing(0).degrees().is_finite());
    }

    #[test]
    fn ring_major_ordering() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            2,
            LengthMeters::new(100.0),
            4,
        );
        let pts: Vec<RingPoint> = g.points().collect();
        assert_eq!((pts[0].ring, pts[0].spoke), (0, 0));
        assert_eq!((pts[3].ring, pts[3].spoke), (0, 3));
        assert_eq!((pts[4].ring, pts[4].spoke), (1, 0));
        let ring1 = g.ring(1);
        assert_eq!(ring1.len(), 4);
        assert_eq!(ring1[0], pts[4].point);
    }

    #[test]
    fn raw_points_keep_copy_radius_metadata() {
        fn require_copy<T: Copy>(value: T) -> T {
            value
        }

        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            1,
            LengthMeters::new(100.0),
            4,
        );
        let raw: RawRingPoint = require_copy(g.raw_ring_point(0, 0));
        assert!(approx(raw.radius, 100.0, 1e-9));
        assert!(approx(raw.with_unit_radius().radius.value(), 100.0, 1e-9));
    }

    #[test]
    fn cardinal_directions_point_the_right_way() {
        let center = GeoPoint::from_lat_lon(0.0, 0.0);
        let g = RadialGrid::new(center, 1, LengthMeters::new(1_000.0), 4);
        let north = g.point(0, 0);
        let east = g.point(0, 1);
        let south = g.point(0, 2);
        let west = g.point(0, 3);
        assert!(north.lat > center.lat && approx(north.lon, 0.0, 1e-6));
        assert!(east.lon > center.lon && approx(east.lat, 0.0, 1e-6));
        assert!(south.lat < center.lat && approx(south.lon, 0.0, 1e-6));
        assert!(west.lon < center.lon && approx(west.lat, 0.0, 1e-6));
    }

    #[test]
    fn points_sit_at_the_ring_radius() {
        let center = GeoPoint::from_lat_lon(52.52, 13.405);
        let g = RadialGrid::new(center, 15, LengthMeters::new(100.0), 32);
        for rp in g.points() {
            let d = haversine_m(center, rp.point);
            assert!(
                approx(d, rp.radius_m(), 0.5),
                "ring {} spoke {}: {d} m vs nominal {} m",
                rp.ring,
                rp.spoke,
                rp.radius_m(),
            );
        }
    }

    #[test]
    fn linear_with_initial_constant_offsets_the_inner_ring() {
        let g = RadialGrid::linear(
            GeoPoint::from_lat_lon(0.0, 0.0),
            3,
            LengthMeters::new(50.0),
            LengthMeters::new(100.0),
            8,
        );
        assert!(approx(g.first_radius().value(), 50.0, 1e-9));
        assert!(approx(g.radius(0).value(), 50.0, 1e-9));
        assert!(approx(g.radius(1).value(), 150.0, 1e-9));
        assert!(approx(g.radius(2).value(), 250.0, 1e-9));
    }

    #[test]
    fn with_first_radius_matches_linear_constructor() {
        let center = GeoPoint::from_lat_lon(52.52, 13.405);
        let a = RadialGrid::new(center, 4, LengthMeters::new(100.0), 8)
            .with_first_radius(LengthMeters::new(50.0));
        let b = RadialGrid::linear(
            center,
            4,
            LengthMeters::new(50.0),
            LengthMeters::new(100.0),
            8,
        );
        assert_eq!(a, b);
    }

    #[test]
    fn new_keeps_first_ring_at_step() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            3,
            LengthMeters::new(100.0),
            8,
        );
        assert!(approx(g.radius(0).value(), 100.0, 1e-9));
        assert!(approx(g.radius(2).value(), 300.0, 1e-9));
    }

    #[test]
    fn bearing_offset_rotates_the_first_spoke() {
        let center = GeoPoint::from_lat_lon(0.0, 0.0);
        let g = RadialGrid::new(center, 1, LengthMeters::new(1_000.0), 4).with_bearing_offset(90.0);
        assert!(approx(g.bearing(0).degrees(), 90.0, 1e-9));
        let p = g.point(0, 0);
        assert!(p.lon > center.lon && approx(p.lat, 0.0, 1e-6));
    }

    #[test]
    fn odd_points_per_ring_is_supported() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            2,
            LengthMeters::new(100.0),
            5,
        );
        assert_eq!(g.len(), 10);
        assert!(approx(g.bearing(1).degrees(), 72.0, 1e-9));
    }

    #[test]
    fn empty_grid_yields_no_points() {
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            0,
            LengthMeters::new(100.0),
            8,
        );
        assert!(g.is_empty());
        assert_eq!(g.points().count(), 0);
        assert_eq!(g.to_points().len(), 0);
    }
}
