//! `s2::ConvexLoop` — a spherical [`Loop`] **validated convex once**, after which
//! every convex-only fast algorithm is *total*: no per-call convexity check, no
//! `Option`, no fallback. The parse-don't-validate companion to [`Loop`]'s
//! runtime [`is_convex`](Loop::is_convex) / auto-dispatching
//! [`contains_point`](Loop::contains_point): a general loop that *happens* to be
//! convex still gets the fast path transparently, while a `ConvexLoop` makes the
//! invariant a compile-time guarantee and unlocks algorithms that are only valid
//! (or only fast) for convex regions.
//!
//! A convex spherical polygon is the intersection of the hemispheres bounded by
//! its edge great circles (`riemann` `Hausdorff/SpherePolygon`). It is
//! geodesically convex — closed under great-circle chords
//! (`sphereInPolygon_convex`) — and by Gauss–Bonnet has interior area `< 2π`, so
//! it lies within an open hemisphere. A `ConvexLoop` is always stored
//! counter-clockwise (interior = the convex, "small" side).

use hausdorff_base::space::s1::Angle;
use hausdorff_planar::algebra::{Orientation, contains_convex};

use super::{Loop, S2Point};

/// A spherical loop known to be convex and wound counter-clockwise.
///
/// Construct with [`ConvexLoop::new`] / [`ConvexLoop::from_vertices`], which
/// validate convexity and normalize the winding. Once built, the convexity
/// invariant holds for the lifetime of the value.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvexLoop {
    /// Invariant: `inner.is_convex()` and `inner.convex_winding() == Some(1.0)`
    /// (counter-clockwise), with `≥ 3` vertices.
    inner: Loop,
}

impl ConvexLoop {
    /// Validate `l` as convex, normalizing to counter-clockwise so the interior
    /// is the convex (small) side. Returns `None` if `l` is non-convex or
    /// degenerate (`< 3` vertices, or collinear).
    ///
    /// A clockwise convex loop is reversed: its left-of-edges interior is the
    /// *large* complement region, but the convex region it bounds is the small
    /// side, which the reversed (CCW) ring describes.
    #[must_use]
    pub fn new(l: Loop) -> Option<Self> {
        match l.convex_winding() {
            Some(s) if s > 0.0 => Some(Self { inner: l }),
            Some(_) => {
                let mut v = l.vertices().to_vec();
                v.reverse();
                let rev = Loop::new(v);
                debug_assert_eq!(rev.convex_winding(), Some(1.0));
                Some(Self { inner: rev })
            }
            None => None,
        }
    }

    /// Validate a vertex ring as a convex loop. See [`new`](Self::new).
    #[must_use]
    pub fn from_vertices(vertices: Vec<S2Point>) -> Option<Self> {
        Self::new(Loop::new(vertices))
    }

    /// Borrow the underlying [`Loop`].
    #[must_use]
    pub fn as_loop(&self) -> &Loop {
        &self.inner
    }

    /// Consume into the underlying [`Loop`].
    #[must_use]
    pub fn into_loop(self) -> Loop {
        self.inner
    }

    /// The loop vertices, counter-clockwise (no repeated closing vertex).
    #[must_use]
    pub fn vertices(&self) -> &[S2Point] {
        self.inner.vertices()
    }

    /// The number of vertices (= the number of edges).
    #[must_use]
    pub fn num_vertices(&self) -> usize {
        self.inner.num_vertices()
    }

    /// Interior area in **steradians** (always `< 2π` — a convex polygon fits in a
    /// hemisphere).
    #[must_use]
    pub fn area(&self) -> f64 {
        self.inner.area()
    }

    /// Total boundary length as an [`Angle`].
    #[must_use]
    pub fn perimeter(&self) -> Angle {
        self.inner.perimeter()
    }

    /// True iff `p` lies in the interior — **O(log n)** via triangle-fan binary
    /// search, the convex-only speedup over the general loop's O(n) winding /
    /// hemisphere scan.
    ///
    /// The vertices of a CCW convex polygon are sorted by angle about the apex
    /// `v₀`, so a binary search locates the fan triangle whose angular wedge holds
    /// `p` (`orient(v₀, vₘ, p)` is monotone), then one orientation test against
    /// that triangle's far edge decides membership. It computes the same predicate
    /// the O(n) hemisphere test does (`riemann`
    /// `Hausdorff/SpherePolygon.sphereInPolygon_hemispheres`: inside ⇔
    /// `orient(vᵢ, vᵢ₊₁, p) > 0` ∀ edge) — the differential `contains` tests pin
    /// the two together. Boundary results (exactly on an edge/vertex) are not
    /// guaranteed. Delegates to the generic
    /// [`contains_convex`](hausdorff_planar::algebra::contains_convex) shared with
    /// the planar `r2::ConvexPolygon`.
    #[must_use]
    pub fn contains_point(&self, p: S2Point) -> bool {
        contains_convex(self.vertices(), p)
    }

    /// True iff every point of `other` lies in `self` — exact convex containment.
    ///
    /// A convex region is closed under great-circle chords
    /// (`riemann` `Hausdorff/SpherePolygon.sphereInPolygon_convex`), so it contains
    /// another convex loop iff it contains that loop's vertices: each edge chord,
    /// hence the whole boundary, hence the (convex) interior follows. `O(m log n)`.
    #[must_use]
    pub fn contains(&self, other: &ConvexLoop) -> bool {
        other.vertices().iter().all(|&q| self.contains_point(q))
    }

    /// True iff `self` and `other` are **disjoint** — the spherical
    /// **separating-axis theorem**: two convex spherical polygons share no point
    /// iff some edge great circle of one has the *entire* other loop strictly on
    /// its outer hemisphere (`orient(vᵢ, vᵢ₊₁, x) = (vᵢ × vᵢ₊₁)·x < 0` for every
    /// vertex `x` of the other). The spherical mirror of
    /// `r2::ConvexPolygon::disjoint`, proof-backed the same way: a convex loop is
    /// the intersection of its edge-plane hemispheres (`riemann`
    /// `Hausdorff/SpherePolygon.sphereInPolygon_hemispheres`), so a loop entirely
    /// outside one of them cannot meet `self`. Valid because a convex loop lies
    /// within a hemisphere (Gauss–Bonnet). A boundary touch counts as intersecting.
    #[must_use]
    pub fn disjoint(&self, other: &ConvexLoop) -> bool {
        separating_edge(self.vertices(), other.vertices())
            || separating_edge(other.vertices(), self.vertices())
    }

    /// True iff `self` and `other` overlap — the negation of [`disjoint`](Self::disjoint).
    #[must_use]
    pub fn intersects(&self, other: &ConvexLoop) -> bool {
        !self.disjoint(other)
    }
}

/// Spherical-SAT half: does some directed edge great circle of the CCW loop `a`
/// put *all* of `b` strictly on its outer hemisphere (`(vᵢ × vᵢ₊₁)·x < 0`)? Such
/// an edge plane separates the two convex loops.
fn separating_edge(a: &[S2Point], b: &[S2Point]) -> bool {
    let n = a.len();
    (0..n).any(|i| {
        let (p, q) = (a[i], a[(i + 1) % n]);
        b.iter().all(|&x| S2Point::orient(p, q, x) < 0.0)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::earth::GeoPoint;
    use hausdorff_planar::algebra::Orientation;
    use std::f64::consts::PI;

    fn sp(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    fn square(half: f64) -> Vec<S2Point> {
        vec![
            sp(-half, -half),
            sp(-half, half),
            sp(half, half),
            sp(half, -half),
        ]
    }

    #[test]
    fn validation_and_normalization() {
        // CCW square accepted as-is.
        let c = ConvexLoop::from_vertices(square(5.0)).expect("convex");
        assert_eq!(c.num_vertices(), 4);
        assert!(c.area() < 2.0 * PI);

        // CW input (reversed) is normalized to CCW — same convex region.
        let mut rev = square(5.0);
        rev.reverse();
        let c2 = ConvexLoop::from_vertices(rev).expect("convex after normalize");
        assert!((c.area() - c2.area()).abs() < 1e-9);
        assert!(c2.contains_point(sp(0.0, 0.0)));

        // Concave L-shape rejected.
        assert!(
            ConvexLoop::from_vertices(vec![
                sp(0.0, 0.0),
                sp(6.0, 0.0),
                sp(6.0, 6.0),
                sp(3.0, 6.0),
                sp(3.0, 3.0),
                sp(0.0, 3.0),
            ])
            .is_none()
        );

        // Degenerate (< 3 vertices) rejected.
        assert!(ConvexLoop::from_vertices(vec![sp(0.0, 0.0), sp(0.0, 5.0)]).is_none());
    }

    #[test]
    fn binary_search_agrees_with_general_loop() {
        // A many-vertex convex polygon (regular-ish ring near the equator) so the
        // O(log n) search exercises more than the trivial base case.
        let n = 16;
        let mut verts = Vec::new();
        for k in 0..n {
            let theta = 2.0 * PI * f64::from(k) / f64::from(n);
            // small circle of angular radius ~10° about (0,0)
            verts.push(sp(10.0 * theta.sin(), 10.0 * theta.cos()));
        }
        let cvx = ConvexLoop::from_vertices(verts.clone()).expect("convex ring");
        let general = cvx.as_loop().clone();

        let tol = 1e-3;
        let (mut inside, mut outside) = (0, 0);
        for i in -200..=200 {
            for j in -200..=200 {
                let p = sp(f64::from(i) * 0.1, f64::from(j) * 0.1);
                // Margin from the proven O(n) hemisphere predicate.
                let v = cvx.vertices();
                let m = v.len();
                let margin = (0..m)
                    .map(|e| S2Point::orient(v[e], v[(e + 1) % m], p))
                    .fold(f64::INFINITY, f64::min);
                if margin > tol {
                    assert!(cvx.contains_point(p), "log-n missed an interior point");
                    assert!(general.contains_point(p), "general loop disagrees inside");
                    inside += 1;
                } else if margin < -tol {
                    assert!(!cvx.contains_point(p), "log-n kept an exterior point");
                    assert!(!general.contains_point(p), "general loop disagrees outside");
                    outside += 1;
                }
            }
        }
        assert!(inside > 100 && outside > 100, "in={inside} out={outside}");
    }

    #[test]
    fn convex_containment_of_loops() {
        let outer = ConvexLoop::from_vertices(square(10.0)).unwrap();
        let inner = ConvexLoop::from_vertices(square(3.0)).unwrap();
        assert!(outer.contains(&inner), "small box inside big box");
        assert!(!inner.contains(&outer), "big box not inside small box");

        // A box straddling the boundary is not contained.
        let straddle =
            ConvexLoop::from_vertices(vec![sp(8.0, 8.0), sp(8.0, 14.0), sp(12.0, 14.0), sp(12.0, 8.0)])
                .unwrap();
        assert!(!outer.contains(&straddle));
    }

    #[test]
    fn sat_disjointness() {
        let a = ConvexLoop::from_vertices(square(3.0)).unwrap();

        // Overlapping: a box straddling a's NE corner ([0,6]² ∩ [-3,3]² = [0,3]²).
        let over =
            ConvexLoop::from_vertices(vec![sp(0.0, 0.0), sp(0.0, 6.0), sp(6.0, 6.0), sp(6.0, 0.0)])
                .unwrap();
        assert!(a.intersects(&over));
        assert!(!a.disjoint(&over));

        // Disjoint: a box well to the east (lon 40–46).
        let far = ConvexLoop::from_vertices(vec![
            sp(-3.0, 40.0),
            sp(-3.0, 46.0),
            sp(3.0, 46.0),
            sp(3.0, 40.0),
        ])
        .unwrap();
        assert!(a.disjoint(&far));
        assert!(!a.intersects(&far));
        assert_eq!(a.disjoint(&far), far.disjoint(&a)); // symmetric

        // Containment ⟹ intersecting (not disjoint).
        let inner = ConvexLoop::from_vertices(square(1.0)).unwrap();
        assert!(a.contains(&inner));
        assert!(a.intersects(&inner));
    }
}
