use super::*;

#[test]
fn line_cover_endpoints_are_vertex_cells_and_path_is_connected() {
    let resolution = 9;
    let vertices = vec![(40.70, -74.00), (40.70, -73.95), (40.72, -73.90)];
    let line = H3Line::new(vertices.clone()).expect("line");
    let cells = line.to_cells(resolution).expect("cells");

    let first = H3Cell::for_point(vertices[0].0, vertices[0].1, resolution).expect("first cell");
    let &(last_lat, last_lng) = vertices.last().expect("last vertex");
    let last = H3Cell::for_point(last_lat, last_lng, resolution).expect("last cell");
    assert_eq!(cells.first().copied(), Some(first));
    assert_eq!(cells.last().copied(), Some(last));

    for &(lat, lng) in &vertices {
        let cell = H3Cell::for_point(lat, lng, resolution).expect("vertex cell");
        assert!(
            cells.contains(&cell),
            "vertex cell {cell} missing from cover"
        );
    }

    let mut unique = cells.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(unique.len(), cells.len(), "cover has duplicate cells");
    for pair in cells.windows(2) {
        assert!(
            pair[0].is_neighbor_with(pair[1]),
            "cover step {} -> {} is not a grid neighbour",
            pair[0],
            pair[1]
        );
    }

    assert_eq!(H3Cell::line_to_cells(&line, resolution), Ok(cells));
}

#[test]
fn line_cover_within_one_cell_is_that_cell() {
    let line = H3Line::new(vec![(40.700, -74.000), (40.701, -74.001)]).expect("line");
    let cell = H3Cell::for_point(40.700, -74.000, 1).expect("cell");
    assert_eq!(line.to_cells(1), Ok(vec![cell]));

    let dot = H3Line::new(vec![(51.5, -0.12)]).expect("dot");
    let dot_cell = H3Cell::for_point(51.5, -0.12, 9).expect("dot cell");
    assert_eq!(dot.to_cells(9), Ok(vec![dot_cell]));
}

#[test]
fn line_cover_matches_concatenated_segment_grid_paths() {
    let resolution = 7;
    let vertices = vec![(40.00, -74.00), (40.00, -73.80), (40.15, -73.65)];
    let line = H3Line::new(vertices.clone()).expect("line");
    let cells = line.to_cells(resolution).expect("cells");

    let vertex_cells: Vec<H3Cell> = vertices
        .iter()
        .map(|&(lat, lng)| H3Cell::for_point(lat, lng, resolution).expect("vertex cell"))
        .collect();
    let mut seen = std::collections::BTreeSet::new();
    let mut expected = Vec::new();
    if seen.insert(vertex_cells[0]) {
        expected.push(vertex_cells[0]);
    }
    for pair in vertex_cells.windows(2) {
        if pair[0] == pair[1] {
            continue;
        }
        for cell in pair[0].grid_path_cells(pair[1]).expect("grid path") {
            if seen.insert(cell) {
                expected.push(cell);
            }
        }
    }
    assert_eq!(cells, expected);
}

#[test]
fn line_cover_validates_input() {
    assert_eq!(
        H3Line::new(Vec::<(f64, f64)>::new()),
        Err(H3LineError::Empty)
    );
    assert_eq!(
        H3Line::new(vec![(40.0, -74.0), (f64::NAN, 0.0)]),
        Err(H3LineError::InvalidCoordinate)
    );
    assert_eq!(
        H3Line::new(vec![(91.0, 0.0)]),
        Err(H3LineError::InvalidCoordinate)
    );
    let line = H3Line::new(vec![(40.0, -74.0), (41.0, -73.0)]).expect("line");
    assert_eq!(line.to_cells(16), Err(H3LineError::InvalidResolution(16)));
}
