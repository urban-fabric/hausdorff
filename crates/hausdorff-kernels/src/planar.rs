//! Planar formula kernels over scalar columns.
//!
//! These functions are the expression-friendly peers of Hausdorff's concrete
//! planar records (`Coord`, `Rect`, `Ellipse`, ...). They take coordinates as
//! scalar values and return the scalar backend's comparison mask, so dataframe
//! and tensor adapters can lower them without constructing eager `f64` geometry
//! values or branching on host `bool`.

use euclid_numeric::{Boolean as BooleanMask, Compare};

/// The Boolean mask produced by comparisons on scalar `T`.
///
/// For eager scalars this is `bool`; for dataframe/tensor/lazy backends this is
/// the backend's Boolean predicate value.
pub type Boolean<T> = <T as Compare>::Mask;

#[inline]
pub(crate) fn mask_and<M: BooleanMask>(a: M, b: M) -> M {
    a.and(&b)
}

#[inline]
pub(crate) fn mask_or<M: BooleanMask>(a: M, b: M) -> M {
    a.or(&b)
}

#[inline]
pub(crate) fn mask_not<M: BooleanMask>(a: M) -> M {
    a.not()
}

mod curves;
pub use curves::{circle_contains_xy, ellipse_contains_xy};

mod orientation;
pub use orientation::{collinear2d, in_triangle2d, left_of_segment2d, orient2d, turns_positive2d};

mod rect;
pub use rect::{rect_contains_xy, rect_distance_squared_xy, rects_overlap_xy};

mod segment;
pub use segment::{
    distance_squared2d, point_on_segment2d, point_segment_distance_squared2d,
    segment_distance_squared2d, segments_intersect2d,
};

#[cfg(test)]
mod tests;
