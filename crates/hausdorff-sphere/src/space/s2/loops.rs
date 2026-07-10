//! `s2::Loop` — a simple closed loop on the **unit sphere**, i.e. a spherical
//! polygon boundary. By convention the interior is the region to the *left* of
//! the directed edges, so a counter-clockwise loop encloses its interior.
//!
//! Everything here is generic unit-sphere geometry: [`area`](Loop::area) is in
//! **steradians** and [`perimeter`](Loop::perimeter) is an [`Angle`], both
//! radius-independent. Convert to metres/m² for a particular planet with
//! [`hausdorff_base::earth::sphere::Sphere`]. Computed from unit vectors — no projection, no
//! pole/antimeridian special cases.

use std::f64::consts::PI;

use hausdorff_base::space::s1::Angle;
use hausdorff_planar::algebra::{Orientation, convex_winding};

use super::edges::{Crossing, crossing_sign};
use super::measures::turn_angle;
use super::{Rect, S2Point, arc};

/// A closed loop of vertices on the unit sphere.
///
/// `convex_sign` caches whether the loop is **geometrically convex** (every
/// vertex turns the same way) and, if so, its turn direction (`+1` CCW, `−1`
/// CW). It is a deterministic function of the vertices, computed once at
/// construction, and gates the proof-backed convex containment fast path (see
/// [`Loop::contains_point_convex`]).
#[derive(Clone, Debug)]
pub struct Loop {
    vertices: Vec<S2Point>,
    convex_sign: Option<f64>,
}

/// Equality is on the vertex ring only; `convex_sign` is a cached function of it.
impl PartialEq for Loop {
    fn eq(&self, other: &Self) -> bool {
        self.vertices == other.vertices
    }
}

impl Loop {
    /// A loop over `vertices`. The closing edge from the last vertex back to the
    /// first is implicit; a trailing vertex equal to the first is dropped.
    #[must_use]
    pub fn new(mut vertices: Vec<S2Point>) -> Self {
        if vertices.len() >= 2 && vertices[0] == *vertices.last().unwrap() {
            vertices.pop();
        }
        // Convexity (and its winding sign) is the generic consistent-turn-sign
        // classification over the spherical `Orientation` predicate `(a×b)·c`.
        let convex_sign = convex_winding(&vertices);
        Self {
            vertices,
            convex_sign,
        }
    }

    /// True iff the loop is **convex** — every vertex turns the same way, so the
    /// region it bounds is geodesically convex (the intersection of the
    /// hemispheres of its edge great circles). Backed by `riemann`
    /// `Hausdorff/SpherePolygon.sphereInPolygon_convex` (the inner-side region is
    /// closed under great-circle chords). By Gauss–Bonnet a convex spherical
    /// polygon has interior area `< 2π`, so it lies within a hemisphere.
    #[must_use]
    pub fn is_convex(&self) -> bool {
        self.convex_sign.is_some()
    }

    /// The cached convex winding sign: `Some(+1)` counter-clockwise, `Some(−1)`
    /// clockwise, `None` when non-convex or degenerate. Used by [`ConvexLoop`] to
    /// validate and normalize a loop to CCW.
    ///
    /// [`ConvexLoop`]: super::ConvexLoop
    #[must_use]
    pub(crate) fn convex_winding(&self) -> Option<f64> {
        self.convex_sign
    }

    /// The loop vertices (no repeated closing vertex).
    #[must_use]
    pub fn vertices(&self) -> &[S2Point] {
        &self.vertices
    }

    /// The number of vertices (= the number of edges).
    #[must_use]
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// The directed edges `(vᵢ, vᵢ₊₁)`, wrapping at the end.
    fn edges(&self) -> impl Iterator<Item = (S2Point, S2Point)> + '_ {
        let n = self.vertices.len();
        (0..n).map(move |i| (self.vertices[i], self.vertices[(i + 1) % n]))
    }

    /// Area of the loop interior in **steradians**, in `[0, 4π]` — the region to
    /// the left of the directed edges. Robust for loops of any size (including
    /// great-circle loops) via the Gauss–Bonnet theorem: `area = 2π − Σ turn`.
    ///
    /// Relative accuracy degrades for sub-square-metre loops (the turning-angle
    /// method has a small absolute error per vertex); for tiny triangles use the
    /// [`super::point_area`] measures directly.
    #[must_use]
    pub fn area(&self) -> f64 {
        let n = self.vertices.len();
        if n < 3 {
            return 0.0;
        }
        let mut turning = 0.0;
        for i in 0..n {
            let prev = self.vertices[(i + n - 1) % n];
            let cur = self.vertices[i];
            let next = self.vertices[(i + 1) % n];
            turning += turn_angle(prev, cur, next).radians();
        }
        (2.0 * PI - turning).clamp(0.0, 4.0 * PI)
    }

    /// True iff `p` lies in the loop interior. Dispatches to the proof-backed
    /// [`contains_point_convex`](Self::contains_point_convex) hemisphere test when
    /// the loop is counter-clockwise convex (the common case — S2 cells, rects,
    /// triangles, cap polygons are all CCW convex), otherwise to the general
    /// spherical winding-number test. Both agree on the interior away from the
    /// boundary; results exactly on an edge/vertex are not guaranteed.
    #[must_use]
    pub fn contains_point(&self, p: S2Point) -> bool {
        if let Some(inside) = self.contains_point_convex(p) {
            return inside;
        }
        self.contains_point_winding(p)
    }

    /// **Convex fast path** — containment by the hemisphere-intersection test.
    /// Returns `Some(inside)` when the loop is counter-clockwise convex and
    /// `None` when it is not (caller should use [`contains_point`](Self::contains_point)).
    ///
    /// A convex spherical polygon is the intersection of the hemispheres bounded
    /// by its edge great circles, so `p` is inside iff it is on the inner (left)
    /// side of *every* directed edge: `orient(vᵢ, vᵢ₊₁, p) = (vᵢ × vᵢ₊₁)·p > 0`.
    /// This is `n` exact dot products with **no transcendentals**, versus the
    /// general winding test's `atan2` per edge. Proof-backed by `riemann`
    /// `Hausdorff/SpherePolygon.sphereInPolygon_hemispheres` (membership ⇔ the
    /// edge-plane hemisphere intersection) — the spherical instance of the generic
    /// [`left_of_all`](hausdorff_planar::algebra::left_of_all) /
    /// `OrientGeom.leftOfAll_convex`.
    ///
    /// Only the counter-clockwise winding is fast-pathed: for a clockwise convex
    /// loop the left-of-edges interior is the *large* complement region, which is
    /// not a hemisphere intersection, so that case defers to the winding test.
    #[must_use]
    pub fn contains_point_convex(&self, p: S2Point) -> Option<bool> {
        match self.convex_sign {
            Some(s) if s > 0.0 => Some(self.edges().all(|(a, b)| S2Point::orient(a, b, p) > 0.0)),
            _ => None,
        }
    }

    /// The general spherical winding-number containment test — the signed solid
    /// angle the loop subtends at `p`: for a counter-clockwise loop an interior
    /// point sees ≈ +2π and an exterior point ≈ 0 or −2π. Robust for loops of any
    /// size or shape; results exactly on an edge/vertex are not guaranteed.
    fn contains_point_winding(&self, p: S2Point) -> bool {
        if self.vertices.len() < 3 {
            return false;
        }
        let pv = p.vector();
        let mut sum = 0.0;
        for (a, b) in self.edges() {
            let c1 = pv.cross(a.vector());
            let c2 = pv.cross(b.vector());
            sum += c1.cross(c2).dot(pv).atan2(c1.dot(c2));
        }
        sum > PI
    }

    /// True iff `p` lies in the loop interior **or on its boundary** — the
    /// boundary-inclusive companion to the open-interior
    /// [`contains_point`](Self::contains_point). Where `contains_point`'s
    /// proof-backed hemisphere test leaves a point exactly on an edge/vertex
    /// "not guaranteed", this resolves it: on-boundary counts as inside. It is
    /// what the OGC predicates that include the boundary (Intersects, Covers,
    /// CoveredBy, Disjoint) test against.
    #[must_use]
    pub fn contains_point_closed(&self, p: S2Point) -> bool {
        self.contains_point(p) || self.edges().any(|(a, b)| point_on_arc(a, b, p))
    }

    /// Perimeter — the sum of the great-circle edge lengths — as an [`Angle`].
    /// Multiply by a planet radius (see [`hausdorff_base::earth::sphere::Sphere`]) for metres.
    #[must_use]
    pub fn perimeter(&self) -> Angle {
        let total: f64 = self
            .edges()
            .map(|(a, b)| a.vector().angle(b.vector()).radians())
            .sum();
        Angle::from_radians(total)
    }

    /// A conservative lat/lng rectangle containing the loop boundary.
    #[must_use]
    pub fn rect_bound(&self) -> Rect {
        if self.vertices.is_empty() {
            return Rect::empty();
        }
        if self.vertices.len() == 1 {
            let g = self.vertices[0].to_geo_point();
            return Rect::from_degrees(g.lat, g.lon, g.lat, g.lon);
        }
        self.edges()
            .map(|(a, b)| arc::rect_bound(a, b))
            .fold(Rect::empty(), |acc, r| acc.union(&r))
    }

    /// True iff this loop and `other` share at least one point (OGC
    /// `Intersects`): their boundaries cross or touch, or one loop contains the
    /// other (the nested case). Boundary-inclusive — coincident edges and
    /// boundary-touching or -coincident vertices count, so a loop intersects
    /// itself. `O(n·m)` in the edge counts.
    #[must_use]
    pub fn intersects(&self, other: &Loop) -> bool {
        if self.vertices.len() < 3 || other.vertices.len() < 3 {
            return false;
        }
        for (a, b) in self.edges() {
            for (c, d) in other.edges() {
                match crossing_sign(a, b, c, d) {
                    Crossing::Cross => return true,
                    // An endpoint lies on the other edge's great circle. Confirm
                    // it lies on the *arc* (segment), not merely the circle,
                    // before counting the boundary touch.
                    Crossing::Maybe => {
                        if point_on_arc(a, b, c)
                            || point_on_arc(a, b, d)
                            || point_on_arc(c, d, a)
                            || point_on_arc(c, d, b)
                        {
                            return true;
                        }
                    }
                    Crossing::DoNotCross => {}
                }
            }
        }
        // Nesting or coincidence: a boundary-coincident vertex counts as shared,
        // so identical loops intersect.
        self.contains_point_closed(other.vertices[0])
            || other.contains_point_closed(self.vertices[0])
    }
}

/// Great-circle-plane tolerance for [`point_on_arc`]: a point within this angle
/// (radians) of an edge's great circle counts as on it. `1e-12 rad` ≈ 6 µm on
/// Earth — far below any real geometry's resolution, but wide enough to absorb
/// the floating-point rounding in the triple product `(a × b) · p`, which is a
/// mathematical zero for a point on the circle (and exactly a vertex) yet lands
/// near `1e-16` in `f64`.
const ON_ARC_EPS: f64 = 1e-12;

/// True iff `p` lies on the great-circle arc `a → b`, endpoints included — the
/// boundary-membership primitive the open-interior `contains_point` omits. `p`
/// must be within [`ON_ARC_EPS`] of the great circle through `a, b` *and* within
/// the minor arc: no farther from either endpoint than they are from each other
/// (`cos` decreases with angle, so a larger dot means a smaller subtended angle).
fn point_on_arc(a: S2Point, b: S2Point, p: S2Point) -> bool {
    let (av, bv, pv) = (a.vector(), b.vector(), p.vector());
    // Normal to the edge plane; `|n| = sin(angle(a, b))`, so `(n · p) / |n|` is
    // `sin` of `p`'s angular distance to the plane. Compare squared to avoid a
    // `sqrt`: `|n·p| ≤ ε·|n|`.
    let n = av.cross(bv);
    let nn = n.dot(n);
    if nn == 0.0 {
        return false; // degenerate edge (a and b antipodal or equal)
    }
    let d = n.dot(pv);
    if d * d > ON_ARC_EPS * ON_ARC_EPS * nn {
        return false;
    }
    let ab = av.dot(bv);
    pv.dot(av) >= ab && pv.dot(bv) >= ab
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::earth::GeoPoint;
    use hausdorff_planar::algebra::Orientation;

    fn sp(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    /// Build the loop, reversing if it captured the large (area > 2π) side, so a
    /// small convex polygon ends up CCW (interior to the left, the small side).
    fn ccw(verts: Vec<S2Point>) -> Loop {
        let l = Loop::new(verts.clone());
        if l.area() > 2.0 * PI {
            let mut v = verts;
            v.reverse();
            Loop::new(v)
        } else {
            l
        }
    }

    #[test]
    fn a_loop_intersects_itself() {
        // The exact failing rect from galois's geometry_invariants proptest.
        let l = ccw(vec![
            sp(-16.898893641227332, -78.04642515443419),
            sp(-16.898893641227332, -76.43908903820845),
            sp(-7.387428552600339, -76.43908903820845),
            sp(-7.387428552600339, -78.04642515443419),
        ]);
        assert!(
            l.contains_point_closed(l.vertices()[0]),
            "a vertex lies on the loop boundary"
        );
        assert!(
            point_on_arc(l.vertices()[0], l.vertices()[1], l.vertices()[0]),
            "an endpoint is on its own edge's arc"
        );
        assert!(
            l.intersects(&l),
            "a non-empty loop shares all its points with itself"
        );
    }

    #[test]
    fn convexity_is_detected() {
        let sq = ccw(vec![
            sp(-5.0, -5.0),
            sp(-5.0, 5.0),
            sp(5.0, 5.0),
            sp(5.0, -5.0),
        ]);
        assert!(sq.is_convex(), "axis-aligned spherical box is convex");
        assert_eq!(sq.convex_sign, Some(1.0), "normalized to CCW (small side)");

        // An L-shape has a reflex vertex → non-convex.
        let l = Loop::new(vec![
            sp(0.0, 0.0),
            sp(6.0, 0.0),
            sp(6.0, 6.0),
            sp(3.0, 6.0),
            sp(3.0, 3.0),
            sp(0.0, 3.0),
        ]);
        assert!(!l.is_convex(), "L-shape is concave");

        // Fewer than three vertices bound nothing.
        assert!(!Loop::new(vec![sp(0.0, 0.0), sp(0.0, 10.0)]).is_convex());
    }

    #[test]
    fn convex_fast_path_only_fires_ccw() {
        let verts = vec![sp(-5.0, -5.0), sp(-5.0, 5.0), sp(5.0, 5.0), sp(5.0, -5.0)];
        let ccw_loop = ccw(verts.clone());
        // CCW convex ⟹ the proof-backed hemisphere test is used.
        assert!(ccw_loop.contains_point_convex(sp(0.0, 0.0)).is_some());

        // The reversed (CW) loop is still geometrically convex, but its
        // left-of-edges interior is the large complement region — not a
        // hemisphere intersection — so the fast path defers to the winding test.
        let mut rev = ccw_loop.vertices().to_vec();
        rev.reverse();
        let cw_loop = Loop::new(rev);
        assert!(cw_loop.is_convex());
        assert_eq!(cw_loop.convex_sign, Some(-1.0));
        assert!(cw_loop.contains_point_convex(sp(0.0, 0.0)).is_none());
    }

    #[test]
    fn convex_fast_path_agrees_with_winding_off_boundary() {
        // A convex spherical quad. Over a lat/lon grid, the proof-backed
        // hemisphere fast path (via `contains_point`) must agree with the general
        // winding-number test everywhere clear of the boundary. "Clear" is
        // decided by the orientation margin: a point with every edge orientation
        // > tol is unambiguously inside, one with an edge < −tol unambiguously
        // outside; the thin band in between is where either test is documented as
        // unreliable, so we skip it.
        let quad = ccw(vec![
            sp(-8.0, -8.0),
            sp(-8.0, 8.0),
            sp(8.0, 8.0),
            sp(8.0, -8.0),
        ]);
        assert!(quad.contains_point_convex(sp(0.0, 0.0)) == Some(true));

        let tol = 1e-3;
        let mut checked_inside = 0;
        let mut checked_outside = 0;
        for i in -250..=250 {
            for j in -250..=250 {
                let p = sp(f64::from(i) * 0.1, f64::from(j) * 0.1);
                let margin = quad
                    .edges()
                    .map(|(a, b)| S2Point::orient(a, b, p))
                    .fold(f64::INFINITY, f64::min);
                if margin > tol {
                    // unambiguously interior
                    assert!(quad.contains_point(p), "interior point reported outside");
                    assert!(quad.contains_point_winding(p), "winding disagrees inside");
                    checked_inside += 1;
                } else if margin < -tol {
                    // unambiguously exterior
                    assert!(!quad.contains_point(p), "exterior point reported inside");
                    assert!(!quad.contains_point_winding(p), "winding disagrees outside");
                    checked_outside += 1;
                }
            }
        }
        assert!(
            checked_inside > 100 && checked_outside > 100,
            "grid should exercise both regions (in={checked_inside}, out={checked_outside})"
        );
    }

    #[test]
    fn rect_bound_contains_loop_boundary_samples() {
        use super::super::arc;

        let l = ccw(vec![
            sp(45.0, 0.0),
            sp(45.0, 90.0),
            sp(10.0, 90.0),
            sp(10.0, 0.0),
        ]);
        let bound = l.rect_bound();
        assert!(bound.hi().lat > 50.0, "bound should include the high arc");
        for (a, b) in l.edges() {
            for i in 0..=100 {
                let p = arc::interpolate(a, b, f64::from(i) / 100.0).to_geo_point();
                assert!(bound.contains_lat_lng(p), "{p:?} escaped {bound:?}");
            }
        }
    }
}
