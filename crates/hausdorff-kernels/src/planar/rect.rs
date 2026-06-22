use euclid::numeric::Real;

use super::{Boolean, mask_and};

/// Inclusive axis-aligned rectangle membership over scalar columns.
#[inline]
#[must_use]
pub fn rect_contains_xy<T: Real>(min_x: T, min_y: T, max_x: T, max_y: T, x: T, y: T) -> Boolean<T> {
    mask_and(
        mask_and(min_x.le(&x), x.le(&max_x)),
        mask_and(min_y.le(&y), y.le(&max_y)),
    )
}

/// Inclusive overlap test for two axis-aligned rectangles.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn rects_overlap_xy<T: Real>(
    a_min_x: T,
    a_min_y: T,
    a_max_x: T,
    a_max_y: T,
    b_min_x: T,
    b_min_y: T,
    b_max_x: T,
    b_max_y: T,
) -> Boolean<T> {
    mask_and(
        mask_and(a_min_x.le(&b_max_x), b_min_x.le(&a_max_x)),
        mask_and(a_min_y.le(&b_max_y), b_min_y.le(&a_max_y)),
    )
}

/// Squared gap between two axis-aligned rectangles, returning zero where they
/// overlap or touch.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn rect_distance_squared_xy<T: Real>(
    a_min_x: T,
    a_min_y: T,
    a_max_x: T,
    a_max_y: T,
    b_min_x: T,
    b_min_y: T,
    b_max_x: T,
    b_max_y: T,
) -> T {
    let zero = T::zero();
    let dx = (b_min_x - a_max_x.clone())
        .max(a_min_x - b_max_x)
        .max(zero.clone());
    let dy = (b_min_y - a_max_y.clone()).max(a_min_y - b_max_y).max(zero);

    dx.clone() * dx + dy.clone() * dy
}
