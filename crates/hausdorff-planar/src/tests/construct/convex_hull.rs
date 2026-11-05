use crate::space::r2::*;

#[test]
fn hull_of_square_with_interior_point() {
    let mp = MultiPoint::from(vec![
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (2.0, 2.0),
    ]);
    let hull = mp.convex_hull();
    assert_eq!(hull.unsigned_area(), 16.0);
    assert!(hull.exterior().is_closed());
    assert_eq!(hull.exterior().coords_count(), 5);
}

#[test]
fn hull_is_counter_clockwise() {
    let mp = MultiPoint::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
    let hull = mp.convex_hull();
    assert!(hull.signed_area() > 0.0);
}
