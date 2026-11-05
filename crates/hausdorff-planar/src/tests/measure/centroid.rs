use crate::space::r2::*;

fn poly(coords: &[(f64, f64)]) -> Polygon {
    Polygon::new(LineString::from(coords.to_vec()), vec![])
}

#[test]
fn square_centroid_is_center() {
    let p = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
    let c = p.centroid().unwrap();
    assert!((c.x() - 1.0).abs() < 1e-12);
    assert!((c.y() - 1.0).abs() < 1e-12);
}

#[test]
fn centroid_is_winding_independent() {
    let ccw = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
    let cw = poly(&[(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0)]);
    assert_eq!(ccw.centroid(), cw.centroid());
}

#[test]
fn multipoint_centroid_is_mean() {
    let mp = MultiPoint::from(vec![(0.0, 0.0), (2.0, 0.0), (1.0, 3.0)]);
    let c = mp.centroid().unwrap();
    assert!((c.x() - 1.0).abs() < 1e-12);
    assert!((c.y() - 1.0).abs() < 1e-12);
}
