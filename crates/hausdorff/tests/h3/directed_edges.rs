use super::*;

#[test]
fn directed_edges_match_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");

        let ours_edges = ours.edges();
        assert_eq!(ours_edges.len(), if ours.is_pentagon() { 5 } else { 6 });

        for edge in &ours_edges {
            let edge = *edge;
            let destination = edge.destination();
            let (origin, destination_from_cells) = edge.cells();

            // Self-consistency: an edge is valid, parses back from its token,
            // its origin is this cell, cells() agrees with destination(), and
            // it is recoverable from (origin, destination).
            assert!(H3DirectedEdge::is_valid(edge.as_u64()));
            assert_eq!(H3DirectedEdge::from_str(&edge.to_string()).ok(), Some(edge));
            assert_eq!(origin, ours);
            assert_eq!(destination_from_cells, destination);
            assert_eq!(edge.origin(), ours);
            assert_eq!(ours.edge(destination), Some(edge));
            assert_eq!(H3DirectedEdge::new(ours, destination), Some(edge));

            f.u64(edge.as_u64());
            f.u64(edge.reverse().as_u64());
            f.u64(destination.as_u64());
            f.u64(edge.origin().as_u64());
        }

        for neighbor in ours.grid_ring(1) {
            f.opt_u64(ours.edge(neighbor).map(H3DirectedEdge::as_u64));
        }
    }
    assert_eq!(f.finish(), golden::EDGES, "directed-edge golden");
}

#[test]
fn directed_edge_geometry_matches_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");
        for edge in ours.edges() {
            for (lat, lon) in edge.boundary_lat_lng() {
                f.f64q(lat);
                f.f64q(lon);
            }
            f.f64q(edge.length_rads());
            f.f64q(edge.length_km());
            f.f64q(edge.length_m());
            // Self-consistency between the typed length and the m accessor.
            assert_eq!(edge.length().value(), edge.length_m());
        }
    }
    assert_eq!(f.finish(), golden::EDGE_GEOM, "edge geometry golden");
}

#[test]
fn directed_edges_reject_invalid_inputs() {
    let cell = H3Cell::for_point(40.758, -73.985, 9).expect("cell");
    let far = H3Cell::for_point(34.0522, -118.2437, 9).expect("far");
    let pentagon = H3Cell::pentagons(1).pop().expect("pentagon");

    assert!(H3DirectedEdge::try_new(0).is_none());
    assert!(!H3DirectedEdge::is_valid(cell.as_u64()));
    assert!(H3DirectedEdge::new(cell, cell).is_none());
    assert!(H3DirectedEdge::new(cell, far).is_none());
    assert!(H3DirectedEdge::from_origin_and_direction(cell, 0).is_none());
    assert!(H3DirectedEdge::from_origin_and_direction(cell, 7).is_none());
    assert!(H3DirectedEdge::from_origin_and_direction(pentagon, 1).is_none());
}
