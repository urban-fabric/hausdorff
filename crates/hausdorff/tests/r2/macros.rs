use hausdorff::space::r2::*;
use hausdorff::{coord, line_string, point, polygon};

#[test]
fn coord_and_point_macros() {
    assert_eq!(coord! { x: 1.0, y: 2.0 }, Coord::new(1.0, 2.0));
    assert_eq!(point! { x: 1.0, y: 2.0 }, Point::new(1.0, 2.0));
}

#[test]
fn line_string_macro_tuples() {
    let ls = line_string![(x: 0.0, y: 0.0), (x: 1.0, y: 1.0)];
    assert_eq!(ls.coords_count(), 2);
    assert_eq!(ls.0[1], Coord::new(1.0, 1.0));
}

#[test]
fn polygon_macro_exterior_only() {
    let p = polygon![(x: 0.0, y: 0.0), (x: 1.0, y: 0.0), (x: 0.0, y: 1.0)];
    assert!(p.exterior().is_closed());
    assert_eq!(p.exterior().coords_count(), 4);
    assert!(p.interiors().is_empty());
}

#[test]
fn polygon_macro_with_hole() {
    let p = polygon!(
        exterior: [(x: 0.0, y: 0.0), (x: 4.0, y: 0.0), (x: 4.0, y: 4.0), (x: 0.0, y: 4.0)],
        interiors: [[(x: 1.0, y: 1.0), (x: 2.0, y: 1.0), (x: 2.0, y: 2.0), (x: 1.0, y: 2.0)]]
    );
    assert_eq!(p.interiors().len(), 1);
    assert!(p.interiors()[0].is_closed());
}
