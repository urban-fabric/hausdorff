//! Shared test support: the builders, comparison helpers and deterministic RNG
//! that more than one suite needs.
//!
//! This exists so a suite can be filed by *subject* without dragging its
//! scaffolding along. Before it, the four suites now split across `measure`,
//! `construct` and `transform` were held in one directory called `r2_query_ops`
//! purely because they shared these helpers — the directory was a consequence of
//! the scaffolding, not a statement about what was being tested.
//!
//! What belongs here is narrow: something used by two or more suites, with one
//! obvious definition. An oracle specific to one suite stays with that suite —
//! the reference closest-point implementation lives next to the closest-point
//! tests, because nothing else can use it.
//!
//! ## Deliberately not consolidated yet
//!
//! Two sets of duplicates are left alone, because merging them is a change to what
//! the tests *do* rather than to where they live, and a half-migration would leave
//! two competing ways to do one thing — which is worse than ten local ones.
//!
//! **Polygon builders.** Ten files define their own `fn poly()`, under four
//! mutually incompatible signatures: two return `Geometry`, two return `Polygon`
//! from flat coordinates, two take rings as slices and three take rings as `Vec`s.
//! Picking a canonical form means rewriting call sites and deciding what the name
//! should mean. `ls` (6 files) and `sq` (5) are the same story.
//!
//! (The LCG that used to be listed here is now [`Lcg`] below.)

use crate::space::r2::{Coord, LineString, Point, Polygon};

/// A polygon's rings as plain coordinate pairs: the exterior, then each hole.
pub(super) type Rings = (Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>);

/// A polygon fixture beside the rings it was built from, so a test can compare a
/// result against the coordinates that produced it.
pub(super) type PolyFixture = (Polygon, Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>);

/// Absolute tolerance for a coordinate comparison in these suites.
pub(super) const EPS: f64 = 1e-9;

/// A `LineString` from plain coordinate pairs.
pub(super) fn line_string(coords: &[(f64, f64)]) -> LineString {
    LineString::new(coords.iter().map(|&(x, y)| Coord::new(x, y)).collect())
}

/// A `Polygon` from an exterior ring and holes, all as plain pairs.
pub(super) fn polygon(exterior: &[(f64, f64)], holes: &[Vec<(f64, f64)>]) -> Polygon {
    Polygon::new(
        line_string(exterior),
        holes.iter().map(|h| line_string(h)).collect(),
    )
}

pub(super) fn closed(ring: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut v = ring.to_vec();
    if let (Some(&first), Some(&last)) = (v.first(), v.last())
        && first != last
    {
        v.push(first);
    }
    v
}

pub(super) fn assert_points_eq(h: Point, r: (f64, f64), ctx: &str) {
    assert!(
        (h.x() - r.0).abs() < EPS && (h.y() - r.1).abs() < EPS,
        "{ctx}: ({}, {}) vs ref ({}, {})",
        h.x(),
        h.y(),
        r.0,
        r.1
    );
}

// ---------------------------------------------------------------------------
// Independent `ClosestPoint` reference oracle.
//
// `RefClosest` mirrors the documented three-way result: the query lies *on* the
// geometry (`Intersection`), a unique nearest spot exists (`SinglePoint`), or
// the geometry is degenerate / empty (`Indeterminate`). Everything below is
// computed from first principles — segment projection and ray-cast
// point-in-ring — with no dependence on hausdorff's own closest-point code.
// ---------------------------------------------------------------------------

pub(super) fn dist2(a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dy) = (a.0 - b.0, a.1 - b.1);
    dx * dx + dy * dy
}

/// True iff `p` lies on the bounded segment `a`–`b` (collinear, within bbox),
/// matching the kernel's exact-zero cross-product test.
pub(super) fn query_points() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0),
        (4.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (5.0, 5.0),
        (-3.0, 1.0),
        (10.0, -10.0),
        (0.5, 3.75),
        (4.0, 4.0),
        (100.0, 0.0),
    ]
}

/// A small deterministic linear congruential generator, so a stress battery is
/// reproducible without a `rand` dependency.
///
/// One definition replacing five. `bool_ops`, `relate`, `min_distance`, `join` and
/// `prepared` each carried the same generator — same multiplier, same increment,
/// same `>> 33` — behind three different wrappers. Two of those wrappers took the
/// raw draw as `u32` and one as `u64`, which looked like a difference and was not:
/// `state >> 33` leaves 31 bits, so both forms hold the same value.
///
/// The *mappings* did genuinely differ, so both survive here under names that say
/// which is which, rather than two incompatible methods both called `coord`.
///
/// `euclid-random` is deliberately not used: it is a contracts crate of key/sampler
/// *traits* and defines no engine, so there is nothing there to call.
pub(crate) struct Lcg(pub(crate) u64);

impl Lcg {
    /// The raw draw: 31 bits, so it is non-negative in any integer type.
    pub(crate) fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) as u32
    }

    /// A coordinate in `[0, n)`.
    pub(crate) fn modulo(&mut self, n: i64) -> f64 {
        (i64::from(self.next_u32()) % n) as f64
    }

    /// A coordinate in `[-n/2, n/2)`, centred on the origin.
    ///
    /// Distinct from [`modulo`](Self::modulo) on purpose: the suites that use this
    /// one want geometry straddling the origin, and swapping the two silently
    /// moves every fixture into the first quadrant.
    pub(crate) fn centred(&mut self, n: i64) -> f64 {
        (i64::from(self.next_u32()) % n) as f64 - (n as f64) / 2.0
    }

    /// A coordinate in `[0, 65536)`, drawn from bits 40..56.
    ///
    /// A third extraction rather than a reuse of [`next_u32`](Self::next_u32):
    /// `query::index`'s own test module was the sixth copy of this generator, and
    /// it takes a different window of the state. Reproducing the window exactly
    /// keeps that suite's point set unchanged while still leaving one LCG step in
    /// the crate.
    pub(crate) fn bits16(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f64::from((self.0 >> 40) as u32 & 0xFFFF)
    }
}
