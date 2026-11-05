use euclid_numeric::Real;

use super::{Boolean, mask_and, mask_not, mask_or};

/// Signed orientation of three planar points.
///
/// A positive value means `(a, b, c)` makes a counter-clockwise turn, a negative
/// value means clockwise, and zero means collinear.
#[inline]
#[must_use]
pub fn orient2d<T: Real>(ax: T, ay: T, bx: T, by: T, cx: T, cy: T) -> T {
    (bx.clone() - ax.clone()) * (cy - ay.clone()) - (by - ay) * (cx - ax)
}

/// Mask that is true where `(a, b, c)` makes a strictly positive turn.
#[inline]
#[must_use]
pub fn turns_positive2d<T: Real>(ax: T, ay: T, bx: T, by: T, cx: T, cy: T) -> Boolean<T> {
    orient2d(ax, ay, bx, by, cx, cy).gt(&T::zero())
}

/// Mask that is true where the absolute orientation is no larger than `tol`.
#[inline]
#[must_use]
pub fn collinear2d<T: Real>(ax: T, ay: T, bx: T, by: T, cx: T, cy: T, tol: T) -> Boolean<T> {
    orient2d(ax, ay, bx, by, cx, cy).abs().le(&tol)
}

/// Mask that is true where `p` lies strictly inside triangle `a, b, c`.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn in_triangle2d<T: Real>(
    ax: T,
    ay: T,
    bx: T,
    by: T,
    cx: T,
    cy: T,
    px: T,
    py: T,
) -> Boolean<T> {
    let zero = T::zero();
    let whole_positive = orient2d(
        ax.clone(),
        ay.clone(),
        bx.clone(),
        by.clone(),
        cx.clone(),
        cy.clone(),
    )
    .gt(&zero);
    let abp = orient2d(
        ax.clone(),
        ay.clone(),
        bx.clone(),
        by.clone(),
        px.clone(),
        py.clone(),
    );
    let bcp = orient2d(
        bx.clone(),
        by.clone(),
        cx.clone(),
        cy.clone(),
        px.clone(),
        py.clone(),
    );
    let cap = orient2d(cx, cy, ax, ay, px, py);

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

/// Mask that is true where `p` lies strictly left of directed segment `a -> b`.
#[inline]
#[must_use]
pub fn left_of_segment2d<T: Real>(ax: T, ay: T, bx: T, by: T, px: T, py: T) -> Boolean<T> {
    turns_positive2d(ax, ay, bx, by, px, py)
}
