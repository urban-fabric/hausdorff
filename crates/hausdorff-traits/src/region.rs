//! Spatial regions — the bounding-volume abstraction shared by every space's
//! `O(1)` containment primitive, and the pruning predicate spatial indexes ride.
//!
//! Each space carries a constant-size region with constant-time point
//! membership: `Interval` on the line (R1) and `Sector` on the circle (S1) — both
//! in `hausdorff-base`; `Rect` / `Ellipse` / `ConfidenceEllipse` in the plane
//! (`hausdorff-planar`); `Aabb` in ℝ³ (`hausdorff-base`); and `Cap` on the sphere
//! (`hausdorff-sphere`). They differ only in their point type, so [`PointRegion`]
//! factors out membership once and the generic spatial filters
//! ([`filter_contained`], [`count_contained`], [`any_contained`]) are written
//! against the trait rather than per space.
//!
//! [`Prunable`] adds the second indexing operation — a pairwise
//! [`intersects`](Prunable::intersects) test with **no false negatives**, so a
//! `false` is always safe to prune on. It is exact for the interval / metric-ball
//! family (`Interval`, `Sector`, `Rect`, `Aabb`, and `Cap`); regions whose
//! pairwise intersection has no closed form (the ellipses) implement only
//! [`PointRegion`].
//!
//! The laws behind each impl are machine-checked in `riemann` (the impls cite the
//! theorems beside their types):
//! - membership — `Hausdorff/Ellipse.inEllipse`, `SphereCap.inCap_iff_dot`,
//!   `Circular.inArc_iff_dot`, and interval/box order facts;
//! - pruning — `ExactArith/AABB.{overlaps_of_mem, prune_sound}` (the R-tree
//!   soundness), `SphereCap.caps_disjoint_of_far`, `Circular.arcs_disjoint_of_far`
//!   (far-apart centres share no point — exactly `!intersects`).

use euclid::numeric::LazyNatural;

/// A spatial region: a set of points of type [`Point`](Self::Point) with
/// constant-time membership. The cross-space generalization of every
/// `contains_*` primitive (`contains_coord`, `contains_point`, `contains_angle`,
/// …). Each impl's membership test is proof-backed in `riemann` (see the module
/// docs).
pub trait PointRegion {
    /// The point type this region is tested against (`f64` for an interval,
    /// `Angle` for a sector, `Coord` for the planar regions, `S2Point` for a cap).
    type Point;

    /// True iff `p` lies in the region (boundary inclusive).
    fn contains(&self, p: Self::Point) -> bool;

    /// True iff the region contains no points.
    fn is_empty(&self) -> bool;
}

/// A [`PointRegion`] that supports a pairwise **intersection** test — the
/// bounding-volume pruning predicate behind R-tree / cap / sector indexes.
///
/// The test never yields a false negative: if two regions share any point it
/// returns `true`, so a `false` result certifies they are disjoint and the pair
/// may be pruned. For the interval / metric-ball family it is *exact*, backed by
/// the `riemann` disjointness proofs (`AABB.prune_sound`,
/// `caps_disjoint_of_far`, `arcs_disjoint_of_far`).
pub trait Prunable: PointRegion {
    /// True iff `self` and `other` may share a point.
    fn intersects(&self, other: &Self) -> bool;
}

/// The points of `points` that lie in `region` — a generic spatial filter over
/// any [`PointRegion`], written once for every space.
#[must_use]
pub fn filter_contained<R>(region: &R, points: &[R::Point]) -> Vec<R::Point>
where
    R: PointRegion,
    R::Point: Clone,
{
    points
        .iter()
        .filter(|p| region.contains((*p).clone()))
        .cloned()
        .collect()
}

/// How many of `points` lie in `region`.
#[must_use]
pub fn count_contained<R, N>(region: &R, points: &[R::Point]) -> N
where
    R: PointRegion,
    R::Point: Clone,
    N: LazyNatural,
{
    points.iter().fold(N::zero(), |count, p| {
        if region.contains((*p).clone()) {
            count + N::one()
        } else {
            count
        }
    })
}

/// True iff any of `points` lies in `region` (`false` for an empty region or no
/// points).
#[must_use]
pub fn any_contained<R>(region: &R, points: &[R::Point]) -> bool
where
    R: PointRegion,
    R::Point: Clone,
{
    !region.is_empty() && points.iter().any(|p| region.contains(p.clone()))
}

/// True iff `a` and `b` provably share no point, so an index may prune the pair.
/// The certified contrapositive of [`Prunable::intersects`] — the `riemann`
/// `*_disjoint_of_far` / `prune_sound` soundness.
#[must_use]
pub fn disjoint<R: Prunable>(a: &R, b: &R) -> bool {
    !a.intersects(b)
}
