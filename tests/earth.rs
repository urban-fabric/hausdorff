//! Public-API tests for `hausdorff::earth` (spherical WGS84), relocated out of `src/`.

mod point {
    use hausdorff::earth::*;

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
            // Round-to-nearest at 1e-7 ⇒ error never exceeds half a step.
            assert!((back.lat - lat).abs() <= 1e-7, "lat {lat} -> {}", back.lat);
            assert!((back.lon - lon).abs() <= 1e-7, "lon {lon} -> {}", back.lon);
        }
    }
}

mod bbox {
    use hausdorff::earth::*;

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
    use hausdorff::earth::*;

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
    use hausdorff::earth::*;

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
    use euclid::LengthMeters;
    use hausdorff::earth::*;

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
        assert!(approx(g.radius(0).as_meters(), 100.0, 1e-9));
        assert!(approx(g.radius(2).as_meters(), 300.0, 1e-9));
        assert!(approx(g.bearing(0).degrees(), 0.0, 1e-9));
        assert!(approx(g.bearing(1).degrees(), 90.0, 1e-9));
        assert!(approx(g.bearing(2).degrees(), 180.0, 1e-9));
        assert!(approx(g.bearing(3).degrees(), 270.0, 1e-9));
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
        // `ring(i)` matches the corresponding flat slice.
        let ring1 = g.ring(1);
        assert_eq!(ring1.len(), 4);
        assert_eq!(ring1[0], pts[4].point);
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
                approx(d, rp.radius.as_meters(), 0.5),
                "ring {} spoke {}: {d} m vs nominal {} m",
                rp.ring,
                rp.spoke,
                rp.radius.as_meters(),
            );
        }
    }

    #[test]
    fn linear_with_initial_constant_offsets_the_inner_ring() {
        // first ring at 50 m, then +100 m per ring -> 50, 150, 250 m.
        let g = RadialGrid::linear(
            GeoPoint::from_lat_lon(0.0, 0.0),
            3,
            LengthMeters::new(50.0),
            LengthMeters::new(100.0),
            8,
        );
        assert!(approx(g.first_radius().as_meters(), 50.0, 1e-9));
        assert!(approx(g.radius(0).as_meters(), 50.0, 1e-9));
        assert!(approx(g.radius(1).as_meters(), 150.0, 1e-9));
        assert!(approx(g.radius(2).as_meters(), 250.0, 1e-9));
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
        // The simple `new` schedule: ring r at (r + 1) * step.
        let g = RadialGrid::new(
            GeoPoint::from_lat_lon(0.0, 0.0),
            3,
            LengthMeters::new(100.0),
            8,
        );
        assert!(approx(g.radius(0).as_meters(), 100.0, 1e-9));
        assert!(approx(g.radius(2).as_meters(), 300.0, 1e-9));
    }

    #[test]
    fn bearing_offset_rotates_the_first_spoke() {
        let center = GeoPoint::from_lat_lon(0.0, 0.0);
        let g = RadialGrid::new(center, 1, LengthMeters::new(1_000.0), 4).with_bearing_offset(90.0);
        assert!(approx(g.bearing(0).degrees(), 90.0, 1e-9));
        let p = g.point(0, 0);
        assert!(p.lon > center.lon && approx(p.lat, 0.0, 1e-6)); // first spoke now points East
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

mod similarity {
    use hausdorff::earth::*;

    fn pt(lat: f64, lon: f64) -> GeoPoint {
        GeoPoint::from_lat_lon(lat, lon)
    }

    #[test]
    fn hausdorff_is_symmetric_and_zero_on_equal_sets() {
        let a = [pt(0.0, 0.0), pt(0.0, 1.0), pt(1.0, 1.0)];
        assert_eq!(hausdorff_distance_m(&a, &a), 0.0);

        let b = [pt(0.0, 0.0), pt(0.0, 1.0)];
        assert!((hausdorff_distance_m(&a, &b) - hausdorff_distance_m(&b, &a)).abs() < 1e-9);
    }

    #[test]
    fn hausdorff_equals_known_great_circle_gap() {
        // b is a's first point plus one extra; the mismatch is the distance from
        // a's lone far point (0,1) to the nearest point of b (0,0).
        let a = [pt(0.0, 0.0), pt(0.0, 1.0)];
        let b = [pt(0.0, 0.0)];
        let expected = haversine_m(pt(0.0, 1.0), pt(0.0, 0.0));
        assert!((hausdorff_distance_m(&a, &b) - expected).abs() < 1e-6);
    }

    #[test]
    fn empty_sets_are_zero() {
        let a = [pt(1.0, 2.0)];
        assert_eq!(hausdorff_distance_m(&a, &[]), 0.0);
        assert_eq!(frechet_distance_m(&[], &a), 0.0);
    }

    #[test]
    fn frechet_zero_on_identical_path() {
        let path = [pt(0.0, 0.0), pt(0.0, 1.0), pt(1.0, 1.0)];
        assert!(frechet_distance_m(&path, &path) < 1e-6);
    }

    #[test]
    fn frechet_of_parallel_offset_lines_is_the_offset() {
        // Two horizontal lines one degree of latitude apart: the leash never
        // needs to exceed that constant north-south gap.
        let lower = [pt(0.0, 0.0), pt(0.0, 1.0), pt(0.0, 2.0)];
        let upper = [pt(1.0, 0.0), pt(1.0, 1.0), pt(1.0, 2.0)];
        let gap = haversine_m(pt(0.0, 0.0), pt(1.0, 0.0));
        assert!((frechet_distance_m(&lower, &upper) - gap).abs() < 1.0);
    }

    #[test]
    fn frechet_dominates_endpoint_distance() {
        // The leash is at least the distance between the two start points and the
        // two end points.
        let p = [pt(0.0, 0.0), pt(2.0, 3.0)];
        let q = [pt(1.0, 1.0), pt(3.0, 5.0)];
        let d = frechet_distance_m(&p, &q);
        assert!(d >= haversine_m(p[0], q[0]) - 1e-6);
        assert!(d >= haversine_m(p[1], q[1]) - 1e-6);
    }

    #[test]
    fn typed_wrappers_match_raw() {
        let a = [pt(0.0, 0.0), pt(0.0, 1.0)];
        let b = [pt(0.5, 0.5)];
        assert_eq!(
            hausdorff_distance(&a, &b).as_meters(),
            hausdorff_distance_m(&a, &b)
        );
        assert_eq!(
            frechet_distance(&a, &b).as_meters(),
            frechet_distance_m(&a, &b)
        );
    }
}
