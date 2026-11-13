//! `r1::Interval` — the real line, validated through the public `space::r1`
//! surface against golang/geo's `r1.Interval` semantics (hand-checked cases).

use crate::space::r1::Interval;

const EMPTY: Interval = Interval::EMPTY;

#[test]
fn empty_is_empty_but_a_point_is_not() {
    assert!(EMPTY.is_empty());
    assert!(Interval::new(1.0, 0.0).is_empty());
    assert!(EMPTY.length() < 0.0);
    // The zero-length interval [0,0] is a *point*, not empty.
    assert!(!Interval::from_point(0.0).is_empty());
    assert!(!Interval::new(0.0, 1.0).is_empty());
}

#[test]
fn accessors() {
    let p = Interval::from_point(3.0);
    assert_eq!((p.lo(), p.hi()), (3.0, 3.0));
    assert_eq!(p.center(), 3.0);
    assert_eq!(p.length(), 0.0);
    let u = Interval::new(0.0, 1.0);
    assert_eq!(u.center(), 0.5);
    assert_eq!(u.length(), 1.0);
}

#[test]
fn contains_includes_endpoints_interior_excludes_them() {
    let u = Interval::new(0.0, 1.0);
    assert!(u.contains(0.0) && u.contains(1.0) && u.contains(0.5));
    assert!(!u.contains(-0.1) && !u.contains(1.1));
    assert!(u.interior_contains(0.5));
    assert!(!u.interior_contains(0.0) && !u.interior_contains(1.0));
}

#[test]
fn contains_interval_and_empty_is_contained_by_all() {
    let u = Interval::new(0.0, 1.0);
    assert!(u.contains_interval(Interval::new(0.2, 0.8)));
    assert!(u.contains_interval(u)); // contains itself
    assert!(!u.contains_interval(Interval::new(-0.1, 0.5)));
    assert!(u.contains_interval(EMPTY));

    assert!(u.interior_contains_interval(Interval::new(0.2, 0.8)));
    assert!(!u.interior_contains_interval(Interval::new(0.0, 0.8))); // touches lo
    assert!(u.interior_contains_interval(EMPTY));
}

#[test]
fn intersects_and_interior_intersects() {
    let u = Interval::new(0.0, 1.0);
    assert!(u.intersects(Interval::new(0.5, 1.5)));
    assert!(u.intersects(Interval::new(-0.5, 0.0))); // boundary touch at 0
    assert!(!u.intersects(Interval::new(1.5, 2.0)));
    assert!(!u.intersects(EMPTY));

    assert!(u.interior_intersects(Interval::new(0.5, 1.5)));
    assert!(!u.interior_intersects(Interval::new(-0.5, 0.0))); // only the boundary
    assert!(!u.interior_intersects(EMPTY));
}

#[test]
fn intersection_and_union() {
    let a = Interval::new(0.0, 1.0);
    let b = Interval::new(0.5, 2.0);
    assert_eq!(a.intersection(b), Interval::new(0.5, 1.0));
    assert!(a.intersection(Interval::new(2.0, 3.0)).is_empty()); // disjoint → empty

    assert_eq!(a.union(b), Interval::new(0.0, 2.0));
    assert_eq!(a.union(EMPTY), a);
    assert_eq!(EMPTY.union(a), a);
    assert!(EMPTY.union(EMPTY).is_empty());
}

#[test]
fn add_point_grows_to_include() {
    assert_eq!(EMPTY.add_point(2.0), Interval::from_point(2.0));
    let u = Interval::new(0.0, 1.0);
    assert_eq!(u.add_point(0.5), u); // already inside
    assert_eq!(u.add_point(2.0), Interval::new(0.0, 2.0));
    assert_eq!(u.add_point(-1.0), Interval::new(-1.0, 1.0));
}

#[test]
fn expanded_and_clamp_point() {
    let u = Interval::new(0.0, 1.0);
    assert_eq!(u.expanded(0.5), Interval::new(-0.5, 1.5));
    assert_eq!(u.expanded(-0.25), Interval::new(0.25, 0.75));
    assert!(EMPTY.expanded(1.0).is_empty()); // empty stays empty

    assert_eq!(u.clamp_point(0.5), 0.5);
    assert_eq!(u.clamp_point(-3.0), 0.0);
    assert_eq!(u.clamp_point(7.0), 1.0);
    // A NaN coordinate stays NaN instead of being clamped to a plausible-looking `hi`.
    assert!(u.clamp_point(f64::NAN).is_nan());
    // The empty interval clamps to its `lo` rather than panicking.
    assert_eq!(EMPTY.clamp_point(0.5), EMPTY.lo);
}

#[test]
fn approx_equal() {
    let u = Interval::new(0.0, 1.0);
    assert!(u.approx_equal(Interval::new(0.0, 1.0 + 1e-16)));
    assert!(!u.approx_equal(Interval::new(0.0, 1.01)));
    assert!(EMPTY.approx_equal(EMPTY));
    assert!(EMPTY.approx_equal(Interval::from_point(1.0))); // length 0 ≤ 2ε
}
