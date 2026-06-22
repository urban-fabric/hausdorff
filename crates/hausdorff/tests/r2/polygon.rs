use hausdorff::space::r2::*;

#[test]
fn polygon_closes_rings() {
    let ext = LineString::new(vec![
        Coord::new(0.0, 0.0),
        Coord::new(1.0, 0.0),
        Coord::new(1.0, 1.0),
    ]);
    let p = Polygon::new(ext, vec![]);
    assert!(p.exterior().is_closed());
    assert_eq!(p.exterior().coords_count(), 4);
    assert!(p.interiors().is_empty());
}

#[test]
fn exterior_mut_recloses() {
    let mut p = Polygon::new(
        LineString::new(vec![
            Coord::new(0.0, 0.0),
            Coord::new(2.0, 0.0),
            Coord::new(1.0, 2.0),
        ]),
        vec![],
    );
    p.exterior_mut(|ring| {
        ring.0.reverse();
    });
    assert!(p.exterior().is_closed());
}
