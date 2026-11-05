use super::*;

#[test]
fn vertices_match_golden() {
    let mut f = Fnv::new();
    for raw in small_sample() {
        let ours = H3Cell::try_new(raw).expect("ours");

        for v in ours.vertexes() {
            f.u64(v.as_u64());
        }
        assert_eq!(
            ours.vertex_count(),
            if ours.is_pentagon() { 5 } else { 6 },
            "vertex count {raw:#x}"
        );

        for vertex in 0..ours.vertex_count() {
            let ours_vertex = ours.vertex(vertex).expect("our vertex");
            let owner = ours_vertex.owner();
            let (lat, lon) = ours_vertex.lat_lng();

            // Self-consistency: the vertex is valid and round-trips through both
            // its raw bits and its textual token.
            assert!(H3Vertex::is_valid(ours_vertex.as_u64()));
            assert_eq!(H3Vertex::try_new(ours_vertex.as_u64()), Some(ours_vertex));
            assert_eq!(
                H3Vertex::from_str(&ours_vertex.to_string()).ok(),
                Some(ours_vertex)
            );

            f.u64(ours_vertex.as_u64());
            f.u64(owner.as_u64());
            f.byte(ours_vertex.vertex());
            f.f64q(lat);
            f.f64q(lon);
        }

        assert!(ours.vertex(ours.vertex_count()).is_none());
    }
    assert_eq!(f.finish(), golden::VERTICES, "vertices golden");
}

#[test]
fn vertices_reject_invalid_inputs() {
    let cell = H3Cell::for_point(40.758, -73.985, 9).expect("cell");
    let pentagon = H3Cell::pentagons(1).pop().expect("pentagon");

    assert!(H3Vertex::try_new(0).is_none());
    assert!(!H3Vertex::is_valid(cell.as_u64()));
    assert!(H3Vertex::new(cell, 6).is_none());
    assert!(H3Vertex::new(pentagon, 5).is_none());
}
