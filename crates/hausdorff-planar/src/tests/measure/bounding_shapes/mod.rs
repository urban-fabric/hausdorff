//! Bounding-shape descriptors, checked against closed forms and enclosure
//! invariants.
//!
//! These have exact answers on simple input, which is unusual and worth
//! exploiting: the minimum bounding circle of a square of side `s` is its
//! circumcircle, radius `s/√2`; the oriented envelope of a rotated rectangle is
//! that rectangle, whatever the rotation. So most of this suite compares against
//! arithmetic rather than against itself.
//!
//! Where a closed form is not available the oracle is **enclosure plus
//! minimality**, which together pin the answer without computing it:
//!
//! - every input coordinate must be inside the shape — otherwise it is not a
//!   bounding shape at all; and
//! - shrinking it must break that — otherwise it is not minimal.
//!
//! The second half is the one worth having. Enclosure alone is satisfied by any
//! shape large enough, so a test that only checked containment would pass for an
//! implementation that returned something absurdly big.
//!
//! The rotation test on the oriented envelope is the sharpest case here, because
//! rotation invariance is exactly the property the axis-aligned bounding box
//! lacks and the reason this descriptor exists. A rotated square's bounding box
//! grows by up to √2; its oriented envelope must not change size at all.
//!
//! Split by descriptor, because the single file had grown past the 400-line house
//! rule — partly from the three exhaustive-enumeration suites added when the
//! inscribed-circle defect was found. The fixtures and the deterministic grid
//! generator the enumeration suites share live here.

use crate::space::r2::{
    Area, BoundingRect, BoundingShapes, Circle, Contains, Coord, CoordsIter, Geometry, LineString,
    MultiPoint, Point, Polygon, oriented_envelope_area,
};

fn square(side: f64) -> Polygon {
    Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (side, 0.0),
            (side, side),
            (0.0, side),
            (0.0, 0.0),
        ]),
        vec![],
    )
}

/// Rotate a ring about the origin, for the invariance checks.
fn rotate(poly: &Polygon, radians: f64) -> Polygon {
    let (s, c) = radians.sin_cos();
    let ring: Vec<Coord> = poly
        .exterior()
        .0
        .iter()
        .map(|p| Coord::new(p.x * c - p.y * s, p.x * s + p.y * c))
        .collect();
    Polygon::new(LineString::new(ring), vec![])
}

fn close(got: f64, want: f64, tol: f64, what: &str) {
    assert!(
        (got - want).abs() <= tol * want.abs().max(1.0),
        "{what}: got {got}, want {want}"
    );
}

/// A deterministic LCG, so a failure is reproducible and bisectable.
///
/// A real RNG would make this suite a lottery: a counterexample found on CI and
/// not locally is worse than no coverage, because it cannot be investigated.
fn lcg(state: &mut u64) -> usize {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    (*state >> 33) as usize
}

/// Integer coordinates on a small grid, which is the adversarial choice: it makes
/// collinear and cocircular triples common rather than rare, and those are the
/// degeneracies that break circle predicates.
fn grid_points(state: &mut u64, max_n: usize) -> Vec<Coord> {
    let n = 1 + lcg(state) % max_n;
    let span = [2_usize, 3, 5, 9][lcg(state) % 4];
    (0..n)
        .map(|_| {
            let (x, y) = (lcg(state) % span, lcg(state) % span);
            Coord::new(x as f64, y as f64)
        })
        .collect()
}

mod circle;
mod degenerate;
mod envelope;
mod inscribed;
