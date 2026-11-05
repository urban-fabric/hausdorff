//! Spherical boolean ops against closed forms, set algebra, and the planar engine.
//!
//! Four independent kinds of oracle, because no single one is sufficient for a
//! boolean engine:
//!
//! 1. **Closed-form areas.** A lune of angular width `w` has area `2w`
//!    steradians. Geometry built from one has an area known before any code runs,
//!    so the result is checked absolutely rather than only for self-consistency.
//! 2. **Set algebra in steradians.** Inclusion–exclusion, idempotence, De Morgan,
//!    and the monotone bounds tie all four operations to each other and to `area`.
//!    An implementation satisfies all of them at once only by being right, and no
//!    reference implementation is needed.
//! 3. **Agreement with the planar engine away from the seams.** Over a small patch
//!    the sphere is nearly flat, so the two must agree to within the curvature of
//!    the patch. That is a genuinely independent oracle — a different algorithm
//!    (Martinez–Rueda sweep) in a different crate.
//! 4. **Rotation invariance across the antimeridian.** The only oracle that shows
//!    the module was worth writing. Longitude is a rotation and the sphere is
//!    rotationally symmetric, so shifting a shape in longitude must not change any
//!    area. In the plane it does.
//!
//! Areas compare with a relative tolerance. `Loop::area` is a Gauss–Bonnet
//! turning-angle sum whose error grows with vertex count, and the boundary is
//! reassembled with more vertices than the inputs had. Structural errors — a
//! dropped component, an inverted ring, a hole on the wrong shell — are off by a
//! whole piece of geometry and fail by orders of magnitude.

use super::{Operation, SphericalBooleanOps, boolean_op, symmetric_difference, total_area};
use crate::space::s2::{Loop, Polygon, S2Point};
use hausdorff_base::earth::GeoPoint;

/// Relative tolerance for an area identity.
const TOL: f64 = 1e-6;

fn p(lat: f64, lon: f64) -> S2Point {
    S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
}

/// A lat/lon "rectangle" as a spherical polygon, counter-clockwise.
///
/// Not a rectangle on the sphere — its north and south edges are great-circle
/// arcs, not parallels — so its area is *not* the graticule formula. That is why
/// the absolute checks use lunes, whose area is exact for great-circle edges.
fn quad(lat0: f64, lon0: f64, lat1: f64, lon1: f64) -> Polygon {
    Polygon::from_shell(Loop::new(vec![
        p(lat0, lon0),
        p(lat0, lon1),
        p(lat1, lon1),
        p(lat1, lon0),
    ]))
}

/// A lune: the wedge between two meridians, pole to pole. Area is exactly
/// `2·width` steradians however it is subdivided.
fn lune(lon0: f64, lon1: f64) -> Polygon {
    let mid = (lon0 + lon1) / 2.0;
    // Counter-clockwise, matching `quad`: east along the south edge, north up the
    // eastern meridian, west along the north edge, south down the western one. The
    // mirror ordering winds clockwise, and `Loop::area` would then report the
    // complement (4π − 2w) rather than the lune.
    Polygon::from_shell(Loop::new(vec![
        p(-89.999_999, lon0),
        p(-89.999_999, mid),
        p(-89.999_999, lon1),
        p(0.0, lon1),
        p(89.999_999, lon1),
        p(89.999_999, mid),
        p(89.999_999, lon0),
        p(0.0, lon0),
    ]))
}

/// A polygonal cap around the north pole: `n` vertices on the parallel at `lat`,
/// in increasing longitude so the pole lies to the left and is enclosed.
///
/// Spelled out as `n` vertices rather than as `quad(lat, -180, lat1, 180)`,
/// which looks like a cap and is not: −180° and +180° are the same meridian, so
/// that quad's four corners are two distinct points and the ring degenerates to a
/// there-and-back sliver. Its area came back as the whole sphere, which then made
/// a union exceed 4π — a broken fixture reading as an engine defect.
fn polar_cap(lat: f64, n: u32) -> Polygon {
    Polygon::from_shell(Loop::new(
        (0..n)
            .map(|k| p(lat, f64::from(k) * 360.0 / f64::from(n)))
            .collect(),
    ))
}

fn close(got: f64, want: f64, what: &str) {
    let bound = TOL * want.abs().max(got.abs()).max(1e-3);
    assert!(
        (got - want).abs() <= bound,
        "{what}: got {got}, want {want} (difference {} exceeds {bound})",
        (got - want).abs()
    );
}

/// Every set-algebra relation, on one pair.
fn assert_set_algebra(name: &str, a: &[Polygon], b: &[Polygon]) {
    let (area_a, area_b) = (total_area(a), total_area(b));
    let u = total_area(&a.union(b));
    let i = total_area(&a.intersection(b));
    let dab = total_area(&a.difference(b));
    let dba = total_area(&b.difference(a));
    let s = total_area(&a.symmetric_difference(b));

    for (label, v) in [
        ("union", u),
        ("intersection", i),
        ("a\\b", dab),
        ("b\\a", dba),
        ("symdiff", s),
    ] {
        assert!(
            v.is_finite() && v >= -TOL,
            "{name}: {label} area must be finite and non-negative, got {v}"
        );
        assert!(
            v <= 4.0 * std::f64::consts::PI + TOL,
            "{name}: {label} area {v} exceeds the whole sphere"
        );
    }

    close(
        u + i,
        area_a + area_b,
        &format!("{name}: |a u b| + |a n b| = |a| + |b|"),
    );
    close(dab + i, area_a, &format!("{name}: |a\\b| + |a n b| = |a|"));
    close(dba + i, area_b, &format!("{name}: |b\\a| + |a n b| = |b|"));
    close(s, u - i, &format!("{name}: |a sym b| = |a u b| - |a n b|"));

    assert!(
        u >= area_a - TOL * area_a.max(1.0) && u >= area_b - TOL * area_b.max(1.0),
        "{name}: union {u} is smaller than an operand ({area_a}, {area_b})"
    );
    assert!(
        i <= area_a + TOL * area_a.max(1.0) && i <= area_b + TOL * area_b.max(1.0),
        "{name}: intersection {i} exceeds an operand ({area_a}, {area_b})"
    );

    // Commutativity — the arrangement treats the operands asymmetrically (one
    // supplies reversed arcs in a difference), so this is not free.
    close(
        total_area(&b.union(a)),
        u,
        &format!("{name}: union commutes"),
    );
    close(
        total_area(&b.intersection(a)),
        i,
        &format!("{name}: intersection commutes"),
    );

    close(
        total_area(&a.union(a)),
        area_a,
        &format!("{name}: a u a = a"),
    );
    close(
        total_area(&a.intersection(a)),
        area_a,
        &format!("{name}: a n a = a"),
    );
    assert!(
        total_area(&a.difference(a)) <= TOL.max(area_a * TOL),
        "{name}: a \\ a must be empty, got {}",
        total_area(&a.difference(a))
    );
}

#[test]
fn overlapping_lunes_have_closed_form_areas() {
    // A lune of width w has area exactly 2w steradians, and the intersection of
    // two overlapping lunes is itself a lune — so every quantity is known in
    // advance, with no reference implementation involved.
    let a = [lune(0.0, 30.0)];
    let b = [lune(20.0, 50.0)];
    close(
        total_area(&a),
        2.0 * 30.0_f64.to_radians(),
        "a lune of width 30 deg has area 2w",
    );
    close(
        total_area(&b),
        2.0 * 30.0_f64.to_radians(),
        "the other lune",
    );
    close(
        total_area(&a.intersection(&b)),
        2.0 * 10.0_f64.to_radians(),
        "the intersection is a lune of width 10 deg",
    );
    close(
        total_area(&a.union(&b)),
        2.0 * 50.0_f64.to_radians(),
        "the union is a lune of width 50 deg",
    );
    close(
        total_area(&a.difference(&b)),
        2.0 * 20.0_f64.to_radians(),
        "a minus b is a lune of width 20 deg",
    );
    close(
        total_area(&symmetric_difference(&a, &b)),
        2.0 * 40.0_f64.to_radians(),
        "the symmetric difference is two lunes of 20 deg",
    );
}

#[test]
fn set_algebra_holds_across_configurations() {
    let cases: Vec<(&str, Vec<Polygon>, Vec<Polygon>)> = vec![
        (
            "overlapping quads at a mid latitude",
            vec![quad(10.0, 10.0, 30.0, 30.0)],
            vec![quad(20.0, 20.0, 40.0, 40.0)],
        ),
        (
            "one quad inside another",
            vec![quad(0.0, 0.0, 40.0, 40.0)],
            vec![quad(10.0, 10.0, 20.0, 20.0)],
        ),
        (
            "identical quads",
            vec![quad(5.0, 5.0, 25.0, 25.0)],
            vec![quad(5.0, 5.0, 25.0, 25.0)],
        ),
        (
            "disjoint quads",
            vec![quad(0.0, 0.0, 10.0, 10.0)],
            vec![quad(50.0, 50.0, 60.0, 60.0)],
        ),
        (
            "edge-sharing quads",
            vec![quad(0.0, 0.0, 20.0, 20.0)],
            vec![quad(0.0, 20.0, 20.0, 40.0)],
        ),
        (
            "overlapping lunes",
            vec![lune(0.0, 40.0)],
            vec![lune(20.0, 60.0)],
        ),
        (
            "quads straddling the antimeridian",
            vec![quad(10.0, 170.0, 30.0, -170.0)],
            vec![quad(20.0, 175.0, 40.0, -160.0)],
        ),
        (
            "two components against one",
            vec![quad(0.0, 0.0, 10.0, 10.0), quad(0.0, 30.0, 10.0, 40.0)],
            vec![quad(5.0, 5.0, 15.0, 35.0)],
        ),
        (
            "southern hemisphere quads",
            vec![quad(-30.0, -30.0, -10.0, -10.0)],
            vec![quad(-20.0, -20.0, 0.0, 0.0)],
        ),
    ];
    for (name, a, b) in cases {
        assert_set_algebra(name, &a, &b);
    }
}

#[test]
fn containment_collapses_union_and_intersection() {
    let outer = [quad(0.0, 0.0, 40.0, 40.0)];
    for inner in [
        vec![quad(10.0, 10.0, 20.0, 20.0)],
        vec![quad(0.0, 0.0, 10.0, 10.0)],
        vec![quad(0.0, 10.0, 10.0, 20.0)],
        vec![quad(0.0, 0.0, 40.0, 40.0)],
    ] {
        let (ai, ao) = (total_area(&inner), total_area(&outer));
        close(
            total_area(&outer.union(&inner)),
            ao,
            "a u b = b when a is inside b",
        );
        close(
            total_area(&outer.intersection(&inner)),
            ai,
            "a n b = a when a is inside b",
        );
        assert!(
            total_area(&inner.difference(&outer)) <= TOL,
            "a \\ b is empty when a is inside b, got {}",
            total_area(&inner.difference(&outer))
        );
        close(
            total_area(&outer.difference(&inner)),
            ao - ai,
            "b \\ a is the annulus",
        );
    }
}

#[test]
fn difference_producing_a_hole_is_a_polygon_with_a_hole() {
    // The structural check an area identity cannot make: subtracting an interior
    // region must give one shell with one hole, not two shells or a dropped hole.
    let outer = [quad(0.0, 0.0, 40.0, 40.0)];
    let inner = [quad(10.0, 10.0, 20.0, 20.0)];
    let holed = outer.difference(&inner);
    assert_eq!(holed.len(), 1, "an annulus is one polygon, got {holed:?}");
    assert_eq!(
        holed[0].holes().len(),
        1,
        "it must carry exactly one hole, got {}",
        holed[0].holes().len()
    );
    close(
        total_area(&holed),
        total_area(&outer) - total_area(&inner),
        "the holed area is the difference",
    );
    assert!(
        !holed[0].contains_point(p(15.0, 15.0)),
        "the hole interior must be excluded"
    );
    assert!(
        holed[0].contains_point(p(5.0, 5.0)),
        "the annulus must be included"
    );
}

#[test]
fn polygons_with_holes_are_handled_as_operands() {
    // The case the spherical predicate used to refuse outright.
    let annulus = Polygon::new(
        Loop::new(vec![p(0.0, 0.0), p(0.0, 40.0), p(40.0, 40.0), p(40.0, 0.0)]),
        vec![Loop::new(vec![
            p(10.0, 10.0),
            p(10.0, 30.0),
            p(30.0, 30.0),
            p(30.0, 10.0),
        ])],
    );
    let a = [annulus];
    let b = [quad(15.0, 15.0, 25.0, 25.0)];
    // `b` sits inside `a`'s hole, so they do not overlap at all.
    assert!(
        total_area(&a.intersection(&b)) <= TOL,
        "a polygon inside the hole intersects nothing, got {}",
        total_area(&a.intersection(&b))
    );
    close(
        total_area(&a.difference(&b)),
        total_area(&a),
        "subtracting something in the hole changes nothing",
    );
    assert_set_algebra("annulus vs quad in its hole", &a, &b);

    let c = [quad(5.0, 5.0, 15.0, 15.0)];
    assert_set_algebra("annulus vs overlapping quad", &a, &c);
    assert!(
        total_area(&a.intersection(&c)) > 0.0,
        "a quad crossing the ring must intersect it"
    );
}

#[test]
fn empty_operands_behave_as_identity_and_annihilator() {
    let empty: Vec<Polygon> = Vec::new();
    let a = [quad(0.0, 0.0, 20.0, 20.0)];
    close(
        total_area(&a.union(&empty)),
        total_area(&a),
        "a u empty = a",
    );
    close(
        total_area(&empty.union(&a)),
        total_area(&a),
        "empty u a = a",
    );
    assert!(total_area(&a.intersection(&empty)) <= TOL);
    assert!(total_area(&empty.intersection(&a)) <= TOL);
    close(
        total_area(&a.difference(&empty)),
        total_area(&a),
        "a \\ empty = a",
    );
    assert!(total_area(&empty.difference(&a)) <= TOL);
    for v in [
        total_area(&empty.union(&empty)),
        total_area(&empty.intersection(&empty)),
        total_area(&empty.difference(&empty)),
        total_area(&symmetric_difference(&empty, &empty)),
    ] {
        assert!(v.abs() <= TOL, "every op on two empties is empty, got {v}");
    }
}

#[test]
fn de_morgan_and_associativity_hold() {
    // Relations among all three operations at once, so a compensating error in two
    // of them cannot satisfy them.
    let a = [quad(0.0, 0.0, 20.0, 20.0)];
    let b = [quad(10.0, 10.0, 30.0, 30.0)];
    let c = [quad(15.0, 0.0, 35.0, 25.0)];
    close(
        total_area(&boolean_op(&a.union(&b), &c, Operation::Difference)),
        total_area(&boolean_op(
            &a.difference(&c),
            &b.difference(&c),
            Operation::Union,
        )),
        "(a u b) \\ c = (a \\ c) u (b \\ c)",
    );
    // Associativity re-runs the arrangement on its own output.
    close(
        total_area(&boolean_op(&a.union(&b), &c, Operation::Union)),
        total_area(&boolean_op(&a, &b.union(&c), Operation::Union)),
        "(a u b) u c = a u (b u c)",
    );
}

#[test]
fn the_antimeridian_is_not_a_boundary() {
    // The reason this module exists. Longitude is a rotation and the sphere is
    // rotationally symmetric, so shifting a shape in longitude cannot change any
    // area. In the plane, crossing +-180 does.
    let reference = total_area(&[quad(10.0, 0.0, 30.0, 20.0)]);
    for shift in [0.0, 45.0, 90.0, 170.0, 175.0, -170.0, 180.0] {
        close(
            total_area(&[quad(10.0, shift, 30.0, shift + 20.0)]),
            reference,
            &format!("a quad shifted {shift} deg keeps its area"),
        );
    }

    // And every operation must be rotation-invariant across the seam: the same
    // overlapping pair at the prime meridian and straddling +-180.
    let (a0, b0) = (
        [quad(10.0, 0.0, 30.0, 20.0)],
        [quad(20.0, 10.0, 40.0, 30.0)],
    );
    let (a1, b1) = (
        [quad(10.0, 170.0, 30.0, -170.0)],
        [quad(20.0, 180.0, 40.0, -160.0)],
    );
    for op in [
        Operation::Union,
        Operation::Intersection,
        Operation::Difference,
    ] {
        close(
            total_area(&boolean_op(&a1, &b1, op)),
            total_area(&boolean_op(&a0, &b0, op)),
            &format!("{op:?} must be rotation-invariant across the seam"),
        );
    }
}

#[test]
fn results_agree_with_the_planar_engine_on_a_small_patch() {
    // An independent oracle: a Martinez-Rueda sweep, in another crate. Over a
    // patch this small the sphere is nearly flat, so the two must agree to within
    // the curvature of the patch — which is what the tolerance represents, not
    // slack in either engine.
    use hausdorff_planar::space::r2::{Area, BooleanOps, LineString, MultiPolygon, Polygon as R2};

    let to_r2 = |lat0: f64, lon0: f64, lat1: f64, lon1: f64| {
        MultiPolygon::new(vec![R2::new(
            LineString::from(vec![(lon0, lat0), (lon1, lat0), (lon1, lat1), (lon0, lat1)]),
            vec![],
        )])
    };
    let (a2, b2) = (to_r2(1.0, 1.0, 2.0, 2.0), to_r2(1.5, 1.5, 2.5, 2.5));
    let (asph, bsph) = ([quad(1.0, 1.0, 2.0, 2.0)], [quad(1.5, 1.5, 2.5, 2.5)]);
    // Square degrees to steradians: one square degree at latitude phi covers
    // (pi/180)^2 * cos(phi).
    let scale = (std::f64::consts::PI / 180.0).powi(2) * 1.75_f64.to_radians().cos();
    for (name, planar, spherical) in [
        (
            "union",
            a2.union(&b2).unsigned_area(),
            total_area(&asph.union(&bsph)),
        ),
        (
            "intersection",
            a2.intersection(&b2).unsigned_area(),
            total_area(&asph.intersection(&bsph)),
        ),
        (
            "difference",
            a2.difference(&b2).unsigned_area(),
            total_area(&asph.difference(&bsph)),
        ),
    ] {
        let expected = planar * scale;
        let relative = (spherical - expected).abs() / expected;
        assert!(
            relative < 2e-3,
            "{name}: spherical {spherical} vs planar-scaled {expected} \
             (relative {relative}) — the engines should agree on a flat patch"
        );
    }
}

// ---- adversarial inputs ----
//
// These found two real defects, both in how orientation was *measured*:
//
//  - input winding was inferred from the sign of `polygon_signed_area_vos`, which
//    is ambiguous above a hemisphere (the van Oosterom sum returns the spherical
//    excess, so `A > 2π` reports `A − 4π` and reads as clockwise). A 200°-wide
//    lune was silently inverted into its own complement, and its union with
//    anything came back smaller than itself.
//  - the same ambiguity classified a single traced output ring as a hole.
//
// Both are now decided structurally rather than by measurement: a shell's
// direction is read off the `Polygon` type, and a lone traced ring is a shell by
// construction because every kept arc carries the result on its left.

/// Every set-algebra relation that must hold for any pair, as a reusable battery.
#[track_caller]
fn assert_relations(name: &str, a: &[Polygon], b: &[Polygon]) {
    let (aa, ab) = (total_area(a), total_area(b));
    let u = total_area(&a.union(b));
    let i = total_area(&a.intersection(b));
    let dab = total_area(&a.difference(b));
    let tol = 1e-6 * (aa + ab).max(1e-9);

    assert!(
        (dab + i - aa).abs() <= tol,
        "{name}: |a\\b| + |a n b| should be |a|: {dab} + {i} vs {aa}"
    );
    assert!(
        u >= aa - tol && u >= ab - tol,
        "{name}: union {u} is smaller than an operand ({aa}, {ab}) — the shape of \
         the inverted-winding defect"
    );
    assert!(
        i <= aa + tol && i <= ab + tol,
        "{name}: intersection {i} exceeds an operand ({aa}, {ab})"
    );
    for v in [u, i, dab] {
        assert!(v.is_finite() && v >= -tol, "{name}: {v} is not a sane area");
        assert!(
            v <= 4.0 * std::f64::consts::PI + tol,
            "{name}: {v} exceeds the whole sphere"
        );
    }
}

#[test]
fn an_operand_larger_than_a_hemisphere_is_not_inverted() {
    // The defect this suite found. A lune of width w has area exactly 2w, so a
    // 200° lune is 6.98 sr against a hemisphere's 6.28 — and its area is the
    // closed form, independent of the engine.
    for span in [200.0_f64, 250.0, 300.0] {
        let wide = [lune(-span / 2.0, span / 2.0)];
        close(
            total_area(&wide),
            2.0 * span.to_radians(),
            &format!("a {span}-deg lune keeps its own area"),
        );
        // Union with something strictly inside it must be the lune itself. Before
        // the fix this returned 4π − area: the complement.
        let inner = [quad(-10.0, -10.0, 10.0, 10.0)];
        close(
            total_area(&wide.union(&inner)),
            total_area(&wide),
            &format!("a {span}-deg lune unioned with a quad inside it is itself"),
        );
        assert!(
            total_area(&wide.union(&inner)) >= total_area(&wide) - 1e-9,
            "a union can never be smaller than an operand"
        );
    }
    // Difference is deliberately *not* asserted above. Subtracting an interior
    // quad gives a lune with a hole — two rings — and a result that is both
    // multi-ring and larger than a hemisphere is the one case `nest` does not
    // claim: with more than one ring it decides shell-versus-hole from
    // `polygon_signed_area_vos`, a fan sum whose apex is the ring's first vertex,
    // and above a hemisphere that apex determines the sign. All three spans here
    // return the complement. See the "Known limit on the output side" note on
    // `nest` and the scope section in this module's docs.
    //
    // What *is* asserted is where that limit starts. Up to a hemisphere the same
    // operation is exact, so the hole-producing path itself is sound — the defect
    // is the measurement above 2π, not the assembly.
    let under = [lune(-85.0, 85.0)];
    let inner = [quad(-10.0, -10.0, 10.0, 10.0)];
    assert!(
        total_area(&under) < 2.0 * std::f64::consts::PI,
        "the 170-deg lune must be under a hemisphere for this to be the supported case"
    );
    close(
        total_area(&under.difference(&inner)),
        total_area(&under) - total_area(&inner),
        "under a hemisphere, a lune minus an interior quad is exact",
    );
}

#[test]
fn set_algebra_survives_adversarial_configurations() {
    // Everything here is under a hemisphere, which is where the engine is fully
    // supported — see the crate docs for the over-hemisphere limit on intersection.
    let holed = Polygon::new(
        Loop::new(vec![p(0.0, 0.0), p(0.0, 30.0), p(30.0, 30.0), p(30.0, 0.0)]),
        vec![Loop::new(vec![
            p(0.0, 10.0),
            p(0.0, 20.0),
            p(15.0, 20.0),
            p(15.0, 10.0),
        ])],
    );
    let many_a: Vec<Polygon> = (0..6)
        .map(|k| quad(0.0, f64::from(k) * 20.0, 10.0, f64::from(k) * 20.0 + 10.0))
        .collect();
    let many_b: Vec<Polygon> = (0..6)
        .map(|k| {
            quad(
                5.0,
                f64::from(k) * 20.0 + 5.0,
                15.0,
                f64::from(k) * 20.0 + 15.0,
            )
        })
        .collect();

    let cases: Vec<(&str, Vec<Polygon>, Vec<Polygon>)> = vec![
        (
            "polar caps overlapping",
            vec![quad(70.0, -180.0, 89.0, 0.0)],
            vec![quad(70.0, -90.0, 89.0, 90.0)],
        ),
        (
            "tiny polygons, 1e-9 deg",
            vec![quad(0.0, 0.0, 1e-9, 1e-9)],
            vec![quad(5e-10, 5e-10, 1.5e-9, 1.5e-9)],
        ),
        (
            "tiny polygons at the snap radius",
            vec![quad(0.0, 0.0, 1e-11, 1e-11)],
            vec![quad(5e-12, 5e-12, 1.5e-11, 1.5e-11)],
        ),
        (
            "thin sliver across the antimeridian",
            vec![quad(0.0, 170.0, 0.001, -170.0)],
            vec![quad(-0.0005, 175.0, 0.0015, -175.0)],
        ),
        (
            "corner touch at the antimeridian",
            vec![quad(0.0, 170.0, 10.0, 180.0)],
            vec![quad(10.0, 180.0, 20.0, -170.0)],
        ),
        ("six components against six", many_a, many_b),
        (
            "hole flush with its shell",
            vec![holed],
            vec![quad(5.0, 5.0, 25.0, 25.0)],
        ),
        (
            "quad over the north pole",
            vec![quad(80.0, -90.0, 85.0, 90.0)],
            vec![quad(82.0, 0.0, 88.0, 180.0)],
        ),
        (
            "a 90-deg lune against a polar cap",
            vec![lune(0.0, 90.0)],
            vec![polar_cap(80.0, 8)],
        ),
    ];
    for (name, a, b) in cases {
        assert_relations(name, &a, &b);
    }
}

#[test]
fn a_lone_traced_ring_is_taken_as_a_shell_not_measured() {
    // The second defect: with one output ring there is nothing to nest, and every
    // kept arc carries the result on its left — so the ring *is* the shell. Reading
    // that off the construction rather than from a signed area is what makes an
    // over-hemisphere result come out right.
    let wide = [lune(-150.0, 150.0)];
    let inner = [quad(-5.0, -5.0, 5.0, 5.0)];
    let u = wide.union(&inner);
    assert_eq!(
        u.len(),
        1,
        "a lune containing a quad unions to one component"
    );
    assert!(
        u[0].holes().is_empty(),
        "and it has no hole — the quad is interior, not subtracted"
    );
    close(total_area(&u), total_area(&wide), "the union is the lune");
}

#[test]
fn a_vertex_lying_on_the_other_operands_edge_is_noded() {
    // A T-junction, and the floating-point reason it is not symmetric in
    // longitude. `polar_cap(.., 8)` puts vertices at 0° and 90°, exactly on the
    // meridian edges of `lune(0, 90)`. At longitude 0 the meridian is exactly the
    // `y = 0` plane (`sin 0° = 0` exactly), so the vertex is an exact degeneracy;
    // at 90° it is off by about 1e-17 radians (`cos 90° = 6.1e-17`), which is
    // five orders of magnitude inside the snap radius but is still a definite
    // non-crossing to an exact predicate. Noding the second one only when the
    // arithmetic happened to be exact left a vertex the walk could enter and not
    // leave, and `trace` dropped every ring: all three operations returned
    // nothing, silently losing both operands.
    //
    // The oracle is a closed form, not a reference run. A regular octagon about
    // the pole is unchanged by a 45° rotation, and the lune is exactly two of
    // those 45° sectors wide, so the overlap is exactly a quarter of the cap's
    // area at *every* rotation — which is also why the rotated fixtures below
    // must agree with the aligned one to the last bit of tolerance.
    let lune90 = [lune(0.0, 90.0)];
    let cap = [polar_cap(80.0, 8)];
    let (area_lune, area_cap) = (total_area(&lune90), total_area(&cap));

    close(area_lune, 2.0 * 90.0_f64.to_radians(), "the lune is 2w");
    assert!(
        area_cap < 2.0 * std::f64::consts::PI * (1.0 - 80.0_f64.to_radians().sin()),
        "an inscribed octagon must be smaller than the true cap, got {area_cap}"
    );

    let quarter = area_cap / 4.0;
    close(
        total_area(&lune90.intersection(&cap)),
        quarter,
        "a 90-deg wedge takes exactly a quarter of an 8-fold symmetric cap",
    );
    close(
        total_area(&lune90.difference(&cap)),
        area_lune - quarter,
        "the lune keeps the rest",
    );
    close(
        total_area(&lune90.union(&cap)),
        area_lune + area_cap - quarter,
        "and the union is inclusion-exclusion",
    );
    assert_relations(
        "lune against a cap with vertices on its meridians",
        &lune90,
        &cap,
    );

    // Rotating the cap moves its vertices off the lune's meridians, which is the
    // configuration that always worked. The aligned case must now give the same
    // answer — if noding depends on which longitude the degeneracy lands at, this
    // is where it shows.
    for offset in [1.0_f64, 5.0, 22.5, 45.0] {
        let rotated = [Polygon::from_shell(Loop::new(
            (0..8u32)
                .map(|k| p(80.0, offset + f64::from(k) * 45.0))
                .collect(),
        ))];
        for op in [
            Operation::Union,
            Operation::Intersection,
            Operation::Difference,
        ] {
            close(
                total_area(&boolean_op(&lune90, &cap, op)),
                total_area(&boolean_op(&lune90, &rotated, op)),
                &format!("{op:?} must not care that the cap is rotated by {offset} deg"),
            );
        }
    }
}

#[test]
fn repeated_calls_on_one_input_give_bit_identical_results() {
    // `classify` groups arcs in a map that it then *iterates*, so that map's order
    // is the order of `kept`, which is the order `trace` seeds its walks in, which
    // fixes the first vertex of every traced ring — and `nest` takes a ring's
    // shell-versus-hole sign from a fan sum apexed at that vertex. Through a
    // `HashMap` the entire chain hung off the hash seed, and `RandomState` draws a
    // fresh key per map instance, so two identical calls in one process really did
    // disagree: a 200°-wide lune minus an interior quad returned 6.858880 on some
    // calls and 5.707491 — its own complement — on others.
    //
    // Determinism is worth a test of its own even where the value is wrong. A
    // stable wrong answer is a documented limit; an unstable one is a different
    // and worse thing, and it hides behind any single-run assertion.
    let cases: Vec<(&str, Vec<Polygon>, Vec<Polygon>)> = vec![
        (
            "over-hemisphere lune minus an interior quad",
            vec![lune(-100.0, 100.0)],
            vec![quad(-10.0, -10.0, 10.0, 10.0)],
        ),
        (
            "overlapping quads",
            vec![quad(10.0, 10.0, 30.0, 30.0)],
            vec![quad(20.0, 20.0, 40.0, 40.0)],
        ),
        (
            "six components against six",
            (0..6)
                .map(|k| quad(0.0, f64::from(k) * 20.0, 10.0, f64::from(k) * 20.0 + 10.0))
                .collect(),
            (0..6)
                .map(|k| {
                    quad(
                        5.0,
                        f64::from(k) * 20.0 + 5.0,
                        15.0,
                        f64::from(k) * 20.0 + 15.0,
                    )
                })
                .collect(),
        ),
        (
            "lune against a polar cap",
            vec![lune(0.0, 90.0)],
            vec![polar_cap(80.0, 8)],
        ),
    ];

    for (name, a, b) in cases {
        for op in [
            Operation::Union,
            Operation::Intersection,
            Operation::Difference,
        ] {
            // A fingerprint of the structure as well as the area: a reordering that
            // left the total intact but split a shell differently would still be a
            // change in the result.
            let fingerprint = || {
                let out = boolean_op(&a, &b, op);
                (
                    out.len(),
                    out.iter().map(|q| q.holes().len()).collect::<Vec<_>>(),
                    total_area(&out).to_bits(),
                )
            };
            let first = fingerprint();
            for trial in 1..16 {
                assert_eq!(
                    fingerprint(),
                    first,
                    "{name}: {op:?} must not depend on hash-map iteration order \
                     (trial {trial} disagreed with the first call)"
                );
            }
        }
    }
}
