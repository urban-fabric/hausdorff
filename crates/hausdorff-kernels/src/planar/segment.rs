use euclid_numeric::Real;

use super::orientation::orient2d;
use super::rect::{rect_contains_xy, rects_overlap_xy};
use super::{Boolean, mask_and};

/// Squared Euclidean distance between two points.
#[inline]
#[must_use]
pub fn distance_squared2d<T: Real>(ax: T, ay: T, bx: T, by: T) -> T {
    let dx = bx - ax;
    let dy = by - ay;
    dx.clone() * dx + dy.clone() * dy
}

/// Mask that is true where `p` lies on the closed segment `a -> b`.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn point_on_segment2d<T: Real>(ax: T, ay: T, bx: T, by: T, px: T, py: T) -> Boolean<T> {
    let min_x = ax.clone().min(bx.clone());
    let min_y = ay.clone().min(by.clone());
    let max_x = ax.clone().max(bx.clone());
    let max_y = ay.clone().max(by.clone());
    let collinear = orient2d(ax, ay, bx, by, px.clone(), py.clone()).eq(&T::zero());

    mask_and(
        collinear,
        rect_contains_xy(min_x, min_y, max_x, max_y, px, py),
    )
}

/// Mask that is true where two closed line segments intersect or touch.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn segments_intersect2d<T: Real>(
    a1x: T,
    a1y: T,
    a2x: T,
    a2y: T,
    b1x: T,
    b1y: T,
    b2x: T,
    b2y: T,
) -> Boolean<T> {
    let a_min_x = a1x.clone().min(a2x.clone());
    let a_min_y = a1y.clone().min(a2y.clone());
    let a_max_x = a1x.clone().max(a2x.clone());
    let a_max_y = a1y.clone().max(a2y.clone());
    let b_min_x = b1x.clone().min(b2x.clone());
    let b_min_y = b1y.clone().min(b2y.clone());
    let b_max_x = b1x.clone().max(b2x.clone());
    let b_max_y = b1y.clone().max(b2y.clone());
    let bbox = rects_overlap_xy(
        a_min_x, a_min_y, a_max_x, a_max_y, b_min_x, b_min_y, b_max_x, b_max_y,
    );

    let zero = T::zero();
    let d1 = orient2d(
        b1x.clone(),
        b1y.clone(),
        b2x.clone(),
        b2y.clone(),
        a1x.clone(),
        a1y.clone(),
    );
    let d2 = orient2d(
        b1x.clone(),
        b1y.clone(),
        b2x.clone(),
        b2y.clone(),
        a2x.clone(),
        a2y.clone(),
    );
    let d3 = orient2d(
        a1x.clone(),
        a1y.clone(),
        a2x.clone(),
        a2y.clone(),
        b1x.clone(),
        b1y.clone(),
    );
    let d4 = orient2d(a1x, a1y, a2x, a2y, b2x, b2y);
    let straddles_b = (d1 * d2).le(&zero);
    let straddles_a = (d3 * d4).le(&zero);

    mask_and(bbox, mask_and(straddles_b, straddles_a))
}

/// Squared distance from point `p` to closed segment `a -> b`.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn point_segment_distance_squared2d<T: Real>(px: T, py: T, ax: T, ay: T, bx: T, by: T) -> T {
    let abx = bx - ax.clone();
    let aby = by - ay.clone();
    let apx = px.clone() - ax.clone();
    let apy = py.clone() - ay.clone();
    let len2 = abx.clone() * abx.clone() + aby.clone() * aby.clone();
    let is_point = len2.eq(&T::zero());
    let safe_len2 = T::select(&is_point, T::one(), len2);
    let t = ((apx * abx.clone() + apy * aby.clone()) / safe_len2).clamp(T::zero(), T::one());
    let proj_x = ax + abx * t.clone();
    let proj_y = ay + aby * t;

    distance_squared2d(px, py, proj_x, proj_y)
}

/// Squared distance between two closed line segments, returning zero where they
/// intersect or touch.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn segment_distance_squared2d<T: Real>(
    a1x: T,
    a1y: T,
    a2x: T,
    a2y: T,
    b1x: T,
    b1y: T,
    b2x: T,
    b2y: T,
) -> T {
    let intersects = segments_intersect2d(
        a1x.clone(),
        a1y.clone(),
        a2x.clone(),
        a2y.clone(),
        b1x.clone(),
        b1y.clone(),
        b2x.clone(),
        b2y.clone(),
    );
    let gap = point_segment_distance_squared2d(
        a1x.clone(),
        a1y.clone(),
        b1x.clone(),
        b1y.clone(),
        b2x.clone(),
        b2y.clone(),
    )
    .min(point_segment_distance_squared2d(
        a2x.clone(),
        a2y.clone(),
        b1x.clone(),
        b1y.clone(),
        b2x.clone(),
        b2y.clone(),
    ))
    .min(point_segment_distance_squared2d(
        b1x.clone(),
        b1y.clone(),
        a1x.clone(),
        a1y.clone(),
        a2x.clone(),
        a2y.clone(),
    ))
    .min(point_segment_distance_squared2d(
        b2x, b2y, a1x, a1y, a2x, a2y,
    ));

    T::select(&intersects, T::zero(), gap)
}
