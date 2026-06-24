use super::*;

#[test]
fn polygon_contains_points_and_holes() {
    let polygon = H3Polygon::new(
        vec![
            (40.70, -74.05),
            (40.70, -73.90),
            (40.85, -73.90),
            (40.85, -74.05),
            (40.70, -74.05),
        ],
        vec![vec![
            (40.74, -74.01),
            (40.74, -73.97),
            (40.78, -73.97),
            (40.78, -74.01),
        ]],
    )
    .expect("polygon");

    assert!(polygon.contains_lat_lng(40.82, -73.98));
    assert!(!polygon.contains_lat_lng(40.76, -73.99));
    assert!(!polygon.contains_lat_lng(41.0, -73.98));
}

#[test]
fn polygon_to_cells_uses_centroid_containment() {
    let polygon = H3Polygon::new(
        vec![
            (40.70, -74.05),
            (40.70, -73.90),
            (40.85, -73.90),
            (40.85, -74.05),
        ],
        Vec::<Vec<(f64, f64)>>::new(),
    )
    .expect("polygon");

    let cells = polygon.to_cells(7).expect("cells");
    assert!(!cells.is_empty());
    assert!(cells.windows(2).all(|pair| pair[0] < pair[1]));
    for cell in &cells {
        let (lat, lng) = cell.center_lat_lng();
        assert!(
            polygon.contains_lat_lng(lat, lng),
            "cell center outside polygon: {cell} ({lat}, {lng})"
        );
    }

    assert_eq!(
        polygon.to_cells(16),
        Err(H3PolygonError::InvalidResolution(16))
    );
}

#[test]
fn cells_to_boundary_rings_dissolves_internal_edges() {
    let cell = H3Cell::for_point(40.758, -73.985, 8).expect("cell");
    let rings = H3Cell::cells_to_boundary_rings([cell]).expect("single ring");
    assert_eq!(rings.len(), 1);
    assert_eq!(rings[0].points().len(), cell.boundary_lat_lng().len());
    assert!((rings[0].area_rads2() - cell.area_rads2()).abs() < 1e-15);

    let neighbor = cell.grid_ring(1).into_iter().next().expect("neighbor");
    let rings = H3Cell::cells_to_boundary_rings([cell, neighbor]).expect("merged ring");
    assert_eq!(rings.len(), 1);
    assert!(rings[0].points().len() < cell.boundary_lat_lng().len() * 2);
    let area = cell.area_rads2() + neighbor.area_rads2();
    assert!((rings[0].area_rads2() - area).abs() < 1e-12);

    assert_eq!(
        H3Cell::cells_to_boundary_rings([cell, cell]),
        Err(H3PolygonError::DuplicateInput)
    );
    assert_eq!(
        H3Cell::cells_to_boundary_rings([cell, cell.center_child(9).expect("child")]),
        Err(H3PolygonError::HeterogeneousResolution)
    );
}

#[test]
fn polygon_rejects_invalid_rings() {
    assert_eq!(
        H3Polygon::new(vec![(0.0, 0.0), (1.0, 1.0)], Vec::<Vec<(f64, f64)>>::new()),
        Err(H3PolygonError::RingTooShort)
    );
    assert_eq!(
        H3Polygon::new(
            vec![(0.0, 200.0), (0.0, 201.0), (1.0, 200.5)],
            Vec::<Vec<(f64, f64)>>::new()
        ),
        Err(H3PolygonError::InvalidCoordinate)
    );
}

#[test]
fn rejects_non_cell_indices() {
    assert!(H3Cell::try_new(0).is_none());
    let valid = H3Cell::for_point(40.7, -74.0, 9).expect("valid cell");
    let raw = valid.as_u64();
    let edge_mode = (raw & !(0xf << 59)) | (2 << 59);
    assert!(H3Cell::try_new(edge_mode).is_none());
    assert!(H3Cell::try_new(raw).is_some());
}
