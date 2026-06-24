use euclid::numeric::Real;

use super::{Boolean, mask_and};

/// Inclusive circle membership over scalar columns.
#[inline]
#[must_use]
pub fn circle_contains_xy<T: Real>(cx: T, cy: T, radius: T, x: T, y: T) -> Boolean<T> {
    let dx = x - cx;
    let dy = y - cy;
    let dist2 = dx.clone() * dx + dy.clone() * dy;
    let radius2 = radius.clone() * radius;
    dist2.le(&radius2)
}

/// Inclusive ellipse membership over scalar columns.
///
/// The test avoids division:
/// `semi_y^2 * dx^2 + semi_x^2 * dy^2 <= semi_x^2 * semi_y^2`.
/// Degenerate ellipses with non-positive semi-axis lengths produce a false mask.
#[inline]
#[must_use]
pub fn ellipse_contains_xy<T: Real>(cx: T, cy: T, semi_x: T, semi_y: T, x: T, y: T) -> Boolean<T> {
    let zero = T::zero();
    let dx = x - cx;
    let dy = y - cy;
    let sx2 = semi_x.clone() * semi_x.clone();
    let sy2 = semi_y.clone() * semi_y.clone();
    let lhs = sy2.clone() * dx.clone() * dx + sx2.clone() * dy.clone() * dy;
    let rhs = sx2 * sy2;
    let nonempty = mask_and(semi_x.gt(&zero), semi_y.gt(&zero));

    mask_and(nonempty, lhs.le(&rhs))
}
