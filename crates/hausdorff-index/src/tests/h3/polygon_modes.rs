use super::*;

fn cell_set(cells: &[H3Cell]) -> std::collections::BTreeSet<u64> {
    cells.iter().map(|cell| cell.as_u64()).collect()
}

/// `(exterior, holes)` rings in `(latitude, longitude)` degrees.
type PolygonRings = (Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>);

#[test]
fn polygon_to_cells_modes_match_golden() {
    let nyc_exterior = vec![
        (40.70, -74.05),
        (40.70, -73.90),
        (40.85, -73.90),
        (40.85, -74.05),
    ];
    let nyc_hole = vec![
        (40.74, -74.01),
        (40.74, -73.97),
        (40.78, -73.97),
        (40.78, -74.01),
    ];
    let berlin_triangle = vec![(52.40, 13.20), (52.60, 13.50), (52.30, 13.60)];

    let cases: Vec<PolygonRings> = vec![
        (nyc_exterior.clone(), Vec::new()),
        (nyc_exterior, vec![nyc_hole]),
        (berlin_triangle, Vec::new()),
    ];

    let mut f = Fnv::new();
    for (case, (exterior, holes)) in cases.into_iter().enumerate() {
        let polygon = H3Polygon::new(exterior, holes).expect("polygon");
        for res in [5u8, 6, 7, 8] {
            let centroid = polygon
                .to_cells_with_mode(res, ContainmentMode::Centroid)
                .expect("centroid cells");
            let full = polygon
                .to_cells_with_mode(res, ContainmentMode::Full)
                .expect("full cells");
            let overlapping = polygon
                .to_cells_with_mode(res, ContainmentMode::Overlapping)
                .expect("overlapping cells");

            // The default mode is Centroid, and the inherent method agrees with
            // the associated function for Full.
            assert_eq!(polygon.to_cells(res).expect("default cells"), centroid);
            assert_eq!(
                H3Cell::polygon_to_cells_with_mode(&polygon, res, ContainmentMode::Full),
                Ok(full.clone())
            );

            let centroid = cell_set(&centroid);
            let full = cell_set(&full);
            let overlapping = cell_set(&overlapping);
            // Metamorphic nesting: full subset centroid subset overlapping.
            assert!(full.is_subset(&centroid), "case {case} res {res}");
            assert!(centroid.is_subset(&overlapping), "case {case} res {res}");

            // Pin each mode's coverage (sorted) to the golden. Order matches the
            // capture: Centroid, Full, Overlapping.
            for set in [&centroid, &full, &overlapping] {
                for &x in set {
                    f.u64(x);
                }
            }
        }
    }
    assert_eq!(f.finish(), golden::POLYGON_MODES, "polygon modes golden");
}

#[test]
fn polygon_modes_reject_invalid_resolution() {
    let polygon = H3Polygon::new(
        vec![(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        Vec::<Vec<(f64, f64)>>::new(),
    )
    .expect("polygon");
    for mode in [
        ContainmentMode::Centroid,
        ContainmentMode::Full,
        ContainmentMode::Overlapping,
    ] {
        assert_eq!(
            polygon.to_cells_with_mode(16, mode),
            Err(H3PolygonError::InvalidResolution(16))
        );
    }
}

#[test]
fn transmeridian_polygon_to_cells_matches_golden() {
    let exterior = vec![
        (-19.0, 178.0),
        (-19.0, -178.0),
        (-15.0, -178.0),
        (-15.0, 178.0),
    ];
    let polygon = H3Polygon::new(exterior, Vec::<Vec<(f64, f64)>>::new()).expect("polygon");
    assert!(polygon.is_transmeridian());

    assert!(polygon.contains_lat_lng(-17.0, 179.5));
    assert!(polygon.contains_lat_lng(-17.0, -179.5));
    assert!(!polygon.contains_lat_lng(-17.0, 170.0));
    assert!(!polygon.contains_lat_lng(-17.0, -170.0));
    assert!(!polygon.contains_lat_lng(-10.0, 179.5));

    let mut f = Fnv::new();
    for res in [2u8, 3, 4, 5] {
        let ours = polygon.to_cells(res).expect("cells");
        // The cover spans both sides of the antimeridian (metamorphic).
        assert!(
            ours.iter().any(|cell| cell.center_lat_lng().1 > 170.0)
                && ours.iter().any(|cell| cell.center_lat_lng().1 < -170.0),
            "cover should span both sides of the date line at res {res}"
        );
        let mut set: Vec<u64> = ours.iter().map(|c| c.as_u64()).collect();
        set.sort_unstable();
        for x in set {
            f.u64(x);
        }
    }
    assert_eq!(f.finish(), golden::TRANSMERIDIAN, "transmeridian golden");
}

#[test]
fn transmeridian_polygon_modes_are_nested() {
    let polygon = H3Polygon::new(
        vec![
            (-19.0, 178.0),
            (-19.0, -178.0),
            (-15.0, -178.0),
            (-15.0, 178.0),
        ],
        Vec::<Vec<(f64, f64)>>::new(),
    )
    .expect("polygon");

    for res in [3u8, 4] {
        let centroid = cell_set(&polygon.to_cells(res).expect("centroid cells"));
        let full = cell_set(
            &polygon
                .to_cells_with_mode(res, ContainmentMode::Full)
                .expect("full cells"),
        );
        let overlapping = cell_set(
            &polygon
                .to_cells_with_mode(res, ContainmentMode::Overlapping)
                .expect("overlapping cells"),
        );

        assert!(!full.is_empty(), "res {res}");
        assert!(full.is_subset(&centroid), "res {res}");
        assert!(centroid.is_subset(&overlapping), "res {res}");

        for raw in &overlapping {
            let cell = H3Cell::try_new(*raw).expect("cell");
            let inside = full.contains(raw)
                || cell
                    .boundary_lat_lng()
                    .iter()
                    .any(|&(lat, lng)| polygon.contains_lat_lng(lat, lng));
            let near_boundary = cell.grid_disk(1).iter().any(|n| {
                let (lat, lng) = n.center_lat_lng();
                polygon.contains_lat_lng(lat, lng)
            });
            assert!(
                inside || near_boundary,
                "overlapping cell {cell} is far from the polygon"
            );
        }
    }
}
