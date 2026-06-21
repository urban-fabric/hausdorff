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

use euclid::numeric::LazyReal;

/// A space carrying a signed orientation predicate on triples.
pub trait Orientation: Clone {
    /// The scalar field the signed orientation is measured in.
    type Scalar: LazyReal;

    /// Signed orientation of `(a, b, c)`: `> 0` a positive turn
    /// (counter-clockwise in R2, right-handed in R3), `< 0` negative,
    /// `0` degenerate (collinear in R2, coplanar in R3).
    fn orient(a: Self, b: Self, c: Self) -> Self::Scalar;

    /// `true` iff `a, b, c` make a strictly positive turn.
    #[inline]
    #[must_use]
    fn turns_positive(a: Self, b: Self, c: Self) -> bool
    where
        Self::Scalar: PartialOrd,
    {
        Self::orient(a, b, c) > Self::Scalar::zero()
    }

    /// `true` iff `a, b, c` are degenerate (collinear / coplanar).
    #[inline]
    #[must_use]
    fn degenerate(a: Self, b: Self, c: Self) -> bool
    where
        Self::Scalar: PartialOrd,
    {
        Self::orient(a, b, c) == Self::Scalar::zero()
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
pub fn is_collinear<O: Orientation>(a: O, b: O, c: O, tol: O::Scalar) -> bool
where
    O::Scalar: PartialOrd,
{
    O::orient(a, b, c).abs() <= tol
}

/// **Strict point-in-triangle**, generic over any [`Orientation`] space.
///
/// `p` is strictly inside triangle `a, b, c` iff it lies on the inner side of
/// all three directed edges — the same-sign orientation test. The barycentric
/// soundness of this certificate is proved in `riemann`
/// (`point_in_triangle_edges`; the spherical mirror `sphereInTriangle`).
#[must_use]
pub fn in_triangle<O: Orientation>(a: O, b: O, c: O, p: O) -> bool
where
    O::Scalar: PartialOrd,
{
    let zero = O::Scalar::zero();
    if O::orient(a.clone(), b.clone(), c.clone()) > zero {
        O::orient(a.clone(), b.clone(), p.clone()) > zero
            && O::orient(b.clone(), c.clone(), p.clone()) > zero
            && O::orient(c, a, p) > zero
    } else {
        O::orient(a.clone(), b.clone(), p.clone()) < zero
            && O::orient(b.clone(), c.clone(), p.clone()) < zero
            && O::orient(c, a, p) < zero
    }
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
pub fn left_of_all<O: Orientation>(vertices: &[O], p: O) -> bool
where
    O::Scalar: PartialOrd,
{
    let n = vertices.len();
    n >= 3
        && (0..n).all(|i| {
            O::orient(
                vertices[i].clone(),
                vertices[(i + 1) % n].clone(),
                p.clone(),
            ) > O::Scalar::zero()
        })
}

/// **Convexity classification** of a closed vertex ring, generic over any
/// [`Orientation`] space: `Some(+1.0)` if every vertex turns left
/// (counter-clockwise), `Some(−1.0)` if every vertex turns right (clockwise),
/// and `None` if the turns disagree (a reflex vertex — non-convex) or the ring is
/// degenerate (`< 3` vertices, or every triple collinear). Collinear (zero) turns
/// are tolerated as long as the nonzero turns all share a sign.
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
pub fn convex_winding<O: Orientation>(vertices: &[O]) -> Option<O::Scalar>
where
    O::Scalar: PartialOrd,
{
    let n = vertices.len();
    if n < 3 {
        return None;
    }
    let zero = O::Scalar::zero();
    let mut pos = false;
    let mut neg = false;
    for i in 0..n {
        let turn = O::orient(
            vertices[i].clone(),
            vertices[(i + 1) % n].clone(),
            vertices[(i + 2) % n].clone(),
        );
        if turn > zero {
            pos = true;
        } else if turn < zero {
            neg = true;
        }
        if pos && neg {
            return None;
        }
    }
    if pos {
        Some(O::Scalar::one())
    } else if neg {
        Some(-O::Scalar::one())
    } else {
        None
    }
}

/// True iff the closed vertex ring is a non-degenerate **convex** polygon — every
/// vertex turns the same way. The convenience wrapper over [`convex_winding`].
#[must_use]
pub fn is_convex_ring<O: Orientation>(vertices: &[O]) -> bool
where
    O::Scalar: PartialOrd,
{
    convex_winding(vertices).is_some()
}

/// **O(log n) point-in-convex-polygon**, generic over any [`Orientation`] space —
/// the logarithmic companion to the O(n) [`left_of_all`].
///
/// `ccw` must be a **counter-clockwise convex** ring with `≥ 3` vertices and no
/// repeated closing vertex (e.g. validated via [`convex_winding`] `== Some(1.0)`).
/// Because its vertices are sorted by angle about the apex `ccw[0]`, a binary
/// search locates the fan triangle whose angular wedge holds `p`
/// (`orient(ccw[0], ccw[mid], p)` is monotone in `mid`), after which one
/// orientation test against that triangle's far edge decides membership. It
/// computes the *same* predicate as `left_of_all` (strictly inside ⇔ left of every
/// edge — `riemann` `sphereInPolygon_hemispheres` / planar `leftOfAllStrict`) in
/// `O(log n)` rather than `O(n)`. Points exactly on the boundary are not
/// guaranteed either way. A non-CCW or non-convex ring gives meaningless results.
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
pub fn contains_convex<O: Orientation>(ccw: &[O], p: O) -> bool
where
    O::Scalar: PartialOrd,
{
    let n = ccw.len();
    if n < 3 {
        return false;
    }
    let zero = O::Scalar::zero();
    // Quick reject outside the apex wedge: p must be left of the first edge
    // ccw[0]→ccw[1] and right of the apex ray ccw[0]→ccw[n−1].
    if O::orient(ccw[0].clone(), ccw[1].clone(), p.clone()) <= zero
        || O::orient(ccw[0].clone(), ccw[n - 1].clone(), p.clone()) >= zero
    {
        return false;
    }
    let mut lo = 1;
    let mut hi = n - 1;
    while hi - lo > 1 {
        let mid = lo.midpoint(hi);
        if O::orient(ccw[0].clone(), ccw[mid].clone(), p.clone()) > zero {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    // Inside iff strictly left of the wedge's far edge ccw[lo]→ccw[hi].
    O::orient(ccw[lo].clone(), ccw[hi].clone(), p) > zero
}
