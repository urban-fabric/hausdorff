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

/// A space carrying a signed orientation predicate on triples.
pub trait Orientation: Copy {
    /// Signed orientation of `(a, b, c)`: `> 0` a positive turn
    /// (counter-clockwise in R2, right-handed in R3), `< 0` negative,
    /// `0` degenerate (collinear in R2, coplanar in R3).
    fn orient(a: Self, b: Self, c: Self) -> f64;

    /// `true` iff `a, b, c` make a strictly positive turn.
    #[inline]
    #[must_use]
    fn turns_positive(a: Self, b: Self, c: Self) -> bool {
        Self::orient(a, b, c) > 0.0
    }

    /// `true` iff `a, b, c` are degenerate (collinear / coplanar).
    #[inline]
    #[must_use]
    fn degenerate(a: Self, b: Self, c: Self) -> bool {
        Self::orient(a, b, c) == 0.0
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
pub fn is_collinear<O: Orientation>(a: O, b: O, c: O, tol: f64) -> bool {
    O::orient(a, b, c).abs() <= tol
}

/// **Strict point-in-triangle**, generic over any [`Orientation`] space.
///
/// `p` is strictly inside triangle `a, b, c` iff it lies on the inner side of
/// all three directed edges — the same-sign orientation test. The barycentric
/// soundness of this certificate is proved in `riemann`
/// (`point_in_triangle_edges`; the spherical mirror `sphereInTriangle`).
#[must_use]
pub fn in_triangle<O: Orientation>(a: O, b: O, c: O, p: O) -> bool {
    if O::orient(a, b, c) > 0.0 {
        O::orient(a, b, p) > 0.0 && O::orient(b, c, p) > 0.0 && O::orient(c, a, p) > 0.0
    } else {
        O::orient(a, b, p) < 0.0 && O::orient(b, c, p) < 0.0 && O::orient(c, a, p) < 0.0
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
pub fn left_of_all<O: Orientation>(vertices: &[O], p: O) -> bool {
    let n = vertices.len();
    n >= 3 && (0..n).all(|i| O::orient(vertices[i], vertices[(i + 1) % n], p) > 0.0)
}
