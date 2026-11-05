use crate::space::r2::*;

#[test]
fn closing_a_ring() {
    let mut ls = LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]);
    assert!(!ls.is_closed());
    ls.close();
    assert!(ls.is_closed());
    assert_eq!(ls.coords_count(), 4);
    ls.close();
    assert_eq!(ls.coords_count(), 4);
}

#[test]
fn lines_iterates_segments() {
    let ls = LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]);
    let segs: Vec<_> = ls.lines().collect();
    assert_eq!(segs.len(), 2);
    assert_eq!(segs[0], Line::new((0.0, 0.0), (1.0, 0.0)));
}
