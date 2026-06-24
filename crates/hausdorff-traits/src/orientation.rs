//! Orientation predicate — the signed turn / volume test shared across the
//! geometry spaces, after the `riemann` Lean `OrientGeom` abstraction.
//!
//! The instances live beside their point types: `Coord` / `Point` (R2, planar
//! signed area `(b−a) × (c−a)`) in `hausdorff-planar`, `Vector` (R3, scalar
//! triple product `(a × b) · c`) in `hausdorff-base`, and `S2Point` (S2) in
//! `hausdorff-sphere`. Both signed measures are bundled here as [`Orientation`]
//! so the generic orientation-based algorithms (containment, hulls, in/out
//! tests) below are written once.
//!
//! The laws the predicate obeys — antisymmetry under a swap, invariance under a
//! cyclic rotation, affinity along a segment — are machine-checked in the
//! `riemann` development (`triDet_swap`/`triDet_cyclic`/`triDet_lerp` for the
//! plane, `triple_*` for ℝ³/S²; unified as `Geometry.OrientGeom`). Any space
//! whose `orient` satisfies them inherits the proved containment/convexity
//! results.

use euclid::numeric::{Boolean as BooleanMask, Compare, Element, Real};

/// The Boolean mask produced by comparisons on scalar `T`.
///
/// For eager scalars this is `bool`; for tensor/lazy backends this is the
/// backend's Boolean predicate value.
pub type Boolean<T> = <T as Compare>::Mask;

#[inline]
fn mask_true<T: Compare>() -> Boolean<T> {
    <Boolean<T> as Element>::one()
}

#[inline]
fn mask_false<T: Compare>() -> Boolean<T> {
    <Boolean<T> as Element>::zero()
}

#[inline]
fn mask_and<M: BooleanMask>(a: M, b: M) -> M {
    a.and(&b)
}

#[inline]
fn mask_or<M: BooleanMask>(a: M, b: M) -> M {
    a.or(&b)
}

#[inline]
fn mask_not<M: BooleanMask>(a: M) -> M {
    a.not()
}

/// Mask-valued convex-winding classification.
///
/// `ccw` and `cw` are mutually exclusive where the ring is non-degenerate and
/// has no mixed turns. `convex` is true where either `ccw` or `cw` is true.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvexWindingMasks<M> {
    /// Every nonzero turn is positive and at least one positive turn exists.
    pub ccw: M,
    /// Every nonzero turn is negative and at least one negative turn exists.
    pub cw: M,
    /// The ring has a consistent nonzero turn direction.
    pub convex: M,
}

/// A space carrying a signed orientation predicate on triples.
pub trait Orientation: Clone {
    /// The scalar field the signed orientation is measured in.
    type Scalar: Real;

    /// Signed orientation of `(a, b, c)`: `> 0` a positive turn
    /// (counter-clockwise in R2, right-handed in R3), `< 0` negative,
    /// `0` degenerate (collinear in R2, coplanar in R3).
    fn orient(a: Self, b: Self, c: Self) -> Self::Scalar;

    /// Boolean mask that is true where `a, b, c` make a strictly positive turn.
    #[inline]
    #[must_use]
    fn turns_positive(a: Self, b: Self, c: Self) -> Boolean<Self::Scalar> {
        Self::orient(a, b, c).gt(&Self::Scalar::zero())
    }

    /// Boolean mask that is true where `a, b, c` are degenerate
    /// (collinear / coplanar).
    #[inline]
    #[must_use]
    fn degenerate(a: Self, b: Self, c: Self) -> Boolean<Self::Scalar> {
        Self::orient(a, b, c).eq(&Self::Scalar::zero())
    }
}

/// **Collinearity** of three planar points — `true` iff `a`, `b`, `c` lie on one
/// line, i.e. the orientation determinant `cross(a, b, c)` is zero within `tol`.
/// The exact `tol = 0` test is [`Orientation::degenerate`]; this is the tolerant
/// convenience for inexact (floating-point) coordinates.
///
/// Proof-backed by `riemann` `Hausdorff/LineGeometry`: `collinear_interp`
/// (points interpolated along `a–b` have `cross = 0`) and its converse
/// `interp_locate_of_collinear` — together `cross a b c = 0 ⇔ c` lies on the
/// line `a–b`.
#[must_use]
pub fn is_collinear<O: Orientation>(a: O, b: O, c: O, tol: O::Scalar) -> Boolean<O::Scalar> {
    O::orient(a, b, c).abs().le(&tol)
}

/// **Strict point-in-triangle**, generic over any [`Orientation`] space.
///
/// `p` is strictly inside triangle `a, b, c` iff it lies on the inner side of
/// all three directed edges — the same-sign orientation test. The barycentric
/// soundness of this certificate is proved in `riemann`
/// (`point_in_triangle_edges`; the spherical mirror `sphereInTriangle`).
#[must_use]
pub fn in_triangle<O: Orientation>(a: O, b: O, c: O, p: O) -> Boolean<O::Scalar> {
    let zero = O::Scalar::zero();
    let whole_positive = O::orient(a.clone(), b.clone(), c.clone()).gt(&zero);
    let abp = O::orient(a.clone(), b.clone(), p.clone());
    let bcp = O::orient(b.clone(), c.clone(), p.clone());
    let cap = O::orient(c, a, p);

    let positive_edges = mask_and(
        mask_and(abp.clone().gt(&zero), bcp.clone().gt(&zero)),
        cap.clone().gt(&zero),
    );
    let negative_edges = mask_and(mask_and(abp.lt(&zero), bcp.lt(&zero)), cap.lt(&zero));

    mask_or(
        mask_and(whole_positive.clone(), positive_edges),
        mask_and(mask_not(whole_positive), negative_edges),
    )
}

/// **Strict point-in-convex-polygon**, generic over any [`Orientation`] space —
/// the n-gon generalization of [`in_triangle`].
///
/// `p` is strictly inside the convex polygon `vertices` (given in
/// counter-clockwise order, `≥ 3` vertices) iff it lies strictly to the **left**
/// of every directed edge `vᵢ → vᵢ₊₁`. The "left of every edge" region is convex
/// — proved once generically in `riemann` (`OrientGeom.leftOfAll_convex`,
/// instantiated planar `leftOfAllStrict_convex` and spherical
/// `sphereLeftOfAll_convex`) — so this is sound for planar polygons and spherical
/// convex polygons alike. Fewer than 3 vertices bound no region (`false`).
#[must_use]
pub fn left_of_all<O: Orientation>(vertices: &[O], p: O) -> Boolean<O::Scalar> {
    let n = vertices.len();
    if n < 3 {
        return mask_false::<O::Scalar>();
    }

    let zero = O::Scalar::zero();
    let mut all_left = mask_true::<O::Scalar>();
    for i in 0..n {
        let left = O::orient(
            vertices[i].clone(),
            vertices[(i + 1) % n].clone(),
            p.clone(),
        )
        .gt(&zero);
        all_left = mask_and(all_left, left);
    }
    all_left
}

/// **Mask-valued convexity classification** of a closed vertex ring, generic
/// over any [`Orientation`] space.
///
/// `ccw` is true where every nonzero turn is positive and at least one positive
/// turn exists; `cw` is the mirrored negative-turn case. Collinear (zero) turns
/// are tolerated as long as the nonzero turns all share a sign. A degenerate
/// ring (`< 3` vertices, or every triple collinear) has both masks false.
///
/// The "inner side of every edge" region of a consistently-wound ring is convex —
/// proved once generically in `riemann` (`OrientGeom.leftOfAll_convex`; planar
/// `leftOfAllStrict_convex`, spherical `sphereInPolygon_convex`). So this
/// classifies exactly the rings on which [`left_of_all`] / [`contains_convex`] are
/// valid, and the sign is the winding direction those tests must be read with.
///
/// **Precondition — simplicity.** For a *simple* (non-self-intersecting) ring,
/// consistent turn signs is equivalent to convexity. A self-intersecting ring can
/// have consistent turns yet not be convex — a pentagram turns the same way at
/// every vertex but winds *twice* (`riemann` notes this for the convex hull, and
/// it is why `Hausdorff/WindingNumber`'s per-edge argument is sound only for
/// `n ≤ 4` without a monotone-turning hypothesis). Callers must ensure the ring is
/// simple; OSM/GeoJSON rings are.
#[must_use]
pub fn convex_winding_masks<O: Orientation>(
    vertices: &[O],
) -> ConvexWindingMasks<Boolean<O::Scalar>> {
    let n = vertices.len();
    if n < 3 {
        let false_mask = mask_false::<O::Scalar>();
        return ConvexWindingMasks {
            ccw: false_mask.clone(),
            cw: false_mask.clone(),
            convex: false_mask,
        };
    }

    let zero = O::Scalar::zero();
    let mut pos = mask_false::<O::Scalar>();
    let mut neg = mask_false::<O::Scalar>();
    for i in 0..n {
        let turn = O::orient(
            vertices[i].clone(),
            vertices[(i + 1) % n].clone(),
            vertices[(i + 2) % n].clone(),
        );
        pos = mask_or(pos, turn.clone().gt(&zero));
        neg = mask_or(neg, turn.lt(&zero));
    }

    let mixed_turns = mask_and(pos.clone(), neg.clone());
    let not_mixed = mask_not(mixed_turns);
    let ccw = mask_and(pos.clone(), not_mixed.clone());
    let cw = mask_and(neg, not_mixed);
    let convex = mask_or(ccw.clone(), cw.clone());

    ConvexWindingMasks { ccw, cw, convex }
}

/// Host-scalar convexity classification of a closed vertex ring.
///
/// This returns `Option` and therefore necessarily branches on concrete scalar
/// masks. Use [`convex_winding_masks`] when writing code that must lower to a
/// lazy/tensor backend.
#[must_use]
pub fn convex_winding<O>(vertices: &[O]) -> Option<O::Scalar>
where
    O: Orientation,
    O::Scalar: Compare<Mask = bool>,
{
    let masks = convex_winding_masks(vertices);
    if masks.ccw {
        Some(O::Scalar::one())
    } else if masks.cw {
        Some(-O::Scalar::one())
    } else {
        None
    }
}

/// True iff the closed vertex ring is a non-degenerate **convex** polygon — every
/// vertex turns the same way. Returns the scalar backend's Boolean mask.
#[must_use]
pub fn is_convex_ring<O: Orientation>(vertices: &[O]) -> Boolean<O::Scalar> {
    convex_winding_masks(vertices).convex
}

/// **Mask-valued point-in-convex-polygon**, generic over any [`Orientation`]
/// space.
///
/// `ccw` must be a **counter-clockwise convex** ring with `≥ 3` vertices and no
/// repeated closing vertex (e.g. validated via [`convex_winding_masks`]'s `ccw`
/// mask, or via [`convex_winding`] `== Some(1.0)` for eager host scalars).
/// It computes the same strict predicate as [`left_of_all`] (inside iff left of
/// every directed edge), returning the scalar backend's [`Boolean`] mask. The
/// implementation intentionally avoids scalar-dependent branching, so tensor
/// backends can lower the result to predicate algebra; it is `O(n)`.
///
/// That left-of-every-edge predicate is proved to **coincide with the topological
/// winding number** on small convex polygons: `riemann`
/// `Hausdorff/WindingNumber.{windingNumber_triangle_interior,
/// windingNumber_quadrilateral_interior}` show a point left of every edge of a CCW
/// triangle / quadrilateral has winding number `1`. So containment here is
/// machine-checked (not only differential-tested) on the dominant small-polygon
/// cases — including the OSM rectangle / building-footprint quads. The general-`n`
/// monotone-turning lemma `windingNumber_eq_one_of_total_lt_three_pi` lifts this
/// to all convex polygons (winding `1` whenever total turning `< 3π`, which a
/// convex interior point's `2π` satisfies — ruling out the pentagram's double
/// wind, the reason per-edge `leftOfAll` alone is insufficient for `n ≥ 5`).
#[must_use]
pub fn contains_convex<O: Orientation>(ccw: &[O], p: O) -> Boolean<O::Scalar> {
    left_of_all(ccw, p)
}
