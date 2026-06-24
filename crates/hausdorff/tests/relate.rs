//! Conformance for the native DE-9IM relate. With geo retired as the oracle,
//! correctness is pinned by (1) an independent cross-check of the `intersects`
//! bit against the dedicated Intersects kernel and (2) the full DE-9IM
//! logical-consistency web — symmetry, the within/contains and coveredby/covers
//! duals across orderings, and the predicate implications — enforced across a
//! battery of fixed pairs and ~3,000 randomized stress pairs.

use hausdorff::space::r2::Intersects as _;
use hausdorff::space::r2::relate::{geometry_dimension, relate};
use hausdorff::{Coord, Geometry, Line, LineString, MultiPolygon, Point, Polygon};

fn poly(coords: &[(f64, f64)]) -> Geometry {
    Geometry::Polygon(Polygon::new(LineString::from(coords.to_vec()), vec![]))
}
fn sq(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    poly(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
}
fn pt(x: f64, y: f64) -> Geometry {
    Geometry::Point(Point::new(x, y))
}
fn line(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    Geometry::Line(Line::new(Coord::new(x0, y0), Coord::new(x1, y1)))
}
fn ls(coords: &[(f64, f64)]) -> Geometry {
    Geometry::LineString(LineString::from(coords.to_vec()))
}

/// Enforce the DE-9IM consistency web for one ordered pair (no oracle): every
/// predicate must be mutually consistent, symmetric/dual under reordering, and
/// the `intersects` bit must agree with the independent Intersects kernel.
fn check(a: &Geometry, b: &Geometry) {
    let m = relate(a, b);
    let mr = relate(b, a);
    let (da, db) = (geometry_dimension(a), geometry_dimension(b));
    let ctx = |name: &str| format!("{name}\n  A={a:?}\n  B={b:?}\n  IM={m:?}");

    // (1) Independent anchor: the relate intersects-bit must equal the dedicated
    // Intersects kernel (a separate implementation), and ⟺ ¬disjoint.
    assert_eq!(
        m.is_intersects(),
        a.intersects(b),
        "{}",
        ctx("intersects vs kernel")
    );
    assert_eq!(
        m.is_intersects(),
        !m.is_disjoint(),
        "{}",
        ctx("intersects == !disjoint")
    );

    // (2) Symmetry of symmetric predicates under reordering.
    assert_eq!(
        mr.is_intersects(),
        m.is_intersects(),
        "{}",
        ctx("intersects symmetric")
    );
    assert_eq!(
        mr.is_disjoint(),
        m.is_disjoint(),
        "{}",
        ctx("disjoint symmetric")
    );
    assert_eq!(
        mr.is_equal_topo(),
        m.is_equal_topo(),
        "{}",
        ctx("equals symmetric")
    );
    assert_eq!(
        mr.is_touches(),
        m.is_touches(),
        "{}",
        ctx("touches symmetric")
    );
    assert_eq!(
        mr.is_overlaps(db, da),
        m.is_overlaps(da, db),
        "{}",
        ctx("overlaps symmetric")
    );
    assert_eq!(
        mr.is_crosses(db, da),
        m.is_crosses(da, db),
        "{}",
        ctx("crosses symmetric")
    );

    // (3) within/contains and coveredby/covers are duals across the ordering.
    assert_eq!(
        m.is_within(),
        mr.is_contains(),
        "{}",
        ctx("within(a,b) == contains(b,a)")
    );
    assert_eq!(
        m.is_coveredby(),
        mr.is_covers(),
        "{}",
        ctx("coveredby(a,b) == covers(b,a)")
    );

    // (4) The predicate implication web.
    if m.is_contains() {
        assert!(
            m.is_intersects() && m.is_covers(),
            "{}",
            ctx("contains ⟹ intersects ∧ covers")
        );
    }
    if m.is_within() {
        assert!(
            m.is_intersects() && m.is_coveredby(),
            "{}",
            ctx("within ⟹ intersects ∧ coveredby")
        );
    }
    if m.is_covers() || m.is_coveredby() || m.is_crosses(da, db) {
        assert!(
            m.is_intersects(),
            "{}",
            ctx("covers/coveredby/crosses ⟹ intersects")
        );
    }
    if m.is_equal_topo() {
        assert!(
            m.is_contains() && m.is_within(),
            "{}",
            ctx("equals ⟹ contains ∧ within")
        );
        assert!(
            !m.is_touches() && !m.is_crosses(da, db) && !m.is_overlaps(da, db),
            "{}",
            ctx("equals ⟹ ¬touches,crosses,overlaps")
        );
    }
    if m.is_touches() {
        assert!(m.is_intersects(), "{}", ctx("touches ⟹ intersects"));
        assert!(
            !m.is_contains() && !m.is_within() && !m.is_overlaps(da, db) && !m.is_crosses(da, db),
            "{}",
            ctx("touches ⟹ ¬contains,within,overlaps,crosses")
        );
    }
    if m.is_overlaps(da, db) {
        assert!(m.is_intersects(), "{}", ctx("overlaps ⟹ intersects"));
        assert!(
            !m.is_contains() && !m.is_within() && !m.is_equal_topo(),
            "{}",
            ctx("overlaps ⟹ ¬contains,within,equals")
        );
    }
    if m.is_disjoint() {
        assert!(
            !m.is_within() && !m.is_contains() && !m.is_covers() && !m.is_coveredby(),
            "{}",
            ctx("disjoint ⟹ no containment")
        );
        assert!(
            !m.is_touches()
                && !m.is_crosses(da, db)
                && !m.is_overlaps(da, db)
                && !m.is_equal_topo(),
            "{}",
            ctx("disjoint ⟹ no touch/cross/overlap/equal")
        );
    }
}

/// Check both orderings.
fn check2(a: &Geometry, b: &Geometry) {
    check(a, b);
    check(b, a);
}

#[path = "relate/fixed.rs"]
mod fixed;

#[path = "relate/stress.rs"]
mod stress;
