use crate::space::r2::*;

fn ls(coords: &[(f64, f64)]) -> LineString {
    LineString::new(coords.iter().map(|&(x, y)| Coord::new(x, y)).collect())
}

#[test]
fn winding_detection() {
    let ccw = ls(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]);
    assert_eq!(ccw.winding_order(), Some(WindingOrder::CounterClockwise));
    assert!(ccw.is_ccw());

    let mut rev = ccw.0.clone();
    rev.reverse();
    assert_eq!(
        LineString::new(rev).winding_order(),
        Some(WindingOrder::Clockwise)
    );
}

#[test]
fn orient_default_makes_exterior_ccw() {
    let cw = Polygon::new(
        ls(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]),
        vec![],
    );
    assert!(cw.exterior().is_cw());
    let oriented = cw.orient(Direction::Default);
    assert!(oriented.exterior().is_ccw());
}
