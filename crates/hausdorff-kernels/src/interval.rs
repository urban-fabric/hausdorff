//! One-dimensional interval kernels over scalar columns.
//!
//! These are the backend-neutral peers of `r1::Interval` predicates. They take
//! interval endpoints and points as scalar values and return comparison masks,
//! so dataframe/tensor adapters can keep filters as predicate expressions.

use euclid::numeric::{Boolean as BooleanMask, Compare, Real};

/// The Boolean mask produced by comparisons on scalar `T`.
pub type Boolean<T> = <T as Compare>::Mask;

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

/// Mask that is true where `[lo, hi]` is empty.
#[inline]
#[must_use]
pub fn interval_is_empty<T: Real>(lo: T, hi: T) -> Boolean<T> {
    lo.gt(&hi)
}

/// Mask that is true where `point` lies in closed interval `[lo, hi]`.
#[inline]
#[must_use]
pub fn interval_contains<T: Real>(lo: T, hi: T, point: T) -> Boolean<T> {
    mask_and(lo.le(&point), point.le(&hi))
}

/// Mask that is true where `point` lies in open interval `(lo, hi)`.
#[inline]
#[must_use]
pub fn interval_interior_contains<T: Real>(lo: T, hi: T, point: T) -> Boolean<T> {
    mask_and(lo.lt(&point), point.lt(&hi))
}

/// Mask that is true where closed intervals `[a_lo, a_hi]` and `[b_lo, b_hi]`
/// share at least one point.
#[inline]
#[must_use]
pub fn intervals_intersect<T: Real>(a_lo: T, a_hi: T, b_lo: T, b_hi: T) -> Boolean<T> {
    let a_nonempty = mask_not(interval_is_empty(a_lo.clone(), a_hi.clone()));
    let b_nonempty = mask_not(interval_is_empty(b_lo.clone(), b_hi.clone()));
    let bounds_overlap = mask_and(a_lo.le(&b_hi), b_lo.le(&a_hi));

    mask_and(mask_and(a_nonempty, b_nonempty), bounds_overlap)
}

/// Mask that is true where interval `[outer_lo, outer_hi]` wholly contains
/// `[inner_lo, inner_hi]`. Empty inner intervals are contained in everything.
#[inline]
#[must_use]
pub fn interval_contains_interval<T: Real>(
    outer_lo: T,
    outer_hi: T,
    inner_lo: T,
    inner_hi: T,
) -> Boolean<T> {
    mask_or(
        interval_is_empty(inner_lo.clone(), inner_hi.clone()),
        mask_and(outer_lo.le(&inner_lo), inner_hi.le(&outer_hi)),
    )
}

/// The closed interval intersection bounds.
///
/// A disjoint pair returns empty bounds (`lo > hi`).
#[inline]
#[must_use]
pub fn interval_intersection<T: Real>(a_lo: T, a_hi: T, b_lo: T, b_hi: T) -> (T, T) {
    (a_lo.max(b_lo), a_hi.min(b_hi))
}

/// The smallest interval containing both points.
#[inline]
#[must_use]
pub fn interval_from_point_pair<T: Real>(a: T, b: T) -> (T, T) {
    (a.clone().min(b.clone()), a.max(b))
}

/// `point` clamped into closed interval `[lo, hi]`.
#[inline]
#[must_use]
pub fn interval_clamp_point<T: Real>(lo: T, hi: T, point: T) -> T {
    lo.max(hi.min(point))
}

/// Non-negative gap between two closed intervals, zero where they intersect.
#[inline]
#[must_use]
pub fn interval_gap<T: Real>(a_lo: T, a_hi: T, b_lo: T, b_hi: T) -> T {
    (b_lo - a_hi.clone()).max(a_lo - b_hi).max(T::zero())
}

/// Select between two intervals using a scalar comparison mask.
#[inline]
#[must_use]
pub fn interval_select<T: Real>(
    mask: &Boolean<T>,
    true_bounds: (T, T),
    false_bounds: (T, T),
) -> (T, T) {
    (
        T::select(mask, true_bounds.0, false_bounds.0),
        T::select(mask, true_bounds.1, false_bounds.1),
    )
}

/// The smallest interval containing both intervals, treating empty intervals as
/// identity elements.
#[inline]
#[must_use]
pub fn interval_union<T: Real>(a_lo: T, a_hi: T, b_lo: T, b_hi: T) -> (T, T) {
    let a_empty = interval_is_empty(a_lo.clone(), a_hi.clone());
    let b_empty = interval_is_empty(b_lo.clone(), b_hi.clone());
    let hull = (
        a_lo.clone().min(b_lo.clone()),
        a_hi.clone().max(b_hi.clone()),
    );
    let when_a_nonempty = interval_select(&b_empty, (a_lo, a_hi), hull);
    interval_select(&a_empty, (b_lo, b_hi), when_a_nonempty)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_masks_match_closed_interval_semantics() {
        assert!(!interval_is_empty::<f64>(0.0, 1.0));
        assert!(interval_is_empty::<f64>(1.0, 0.0));
        assert!(interval_contains::<f64>(0.0, 2.0, 2.0));
        assert!(!interval_interior_contains::<f64>(0.0, 2.0, 2.0));

        assert!(intervals_intersect::<f64>(0.0, 2.0, 2.0, 4.0));
        assert!(!intervals_intersect::<f64>(0.0, 2.0, 3.0, 4.0));
        assert!(!intervals_intersect::<f64>(5.0, 4.0, 0.0, 10.0));

        assert!(interval_contains_interval::<f64>(0.0, 4.0, 1.0, 3.0));
        assert!(interval_contains_interval::<f64>(0.0, 4.0, 3.0, 2.0));
        assert!(!interval_contains_interval::<f64>(0.0, 4.0, -1.0, 3.0));
    }

    #[test]
    fn interval_value_kernels_match_known_values() {
        assert_eq!(interval_intersection(0.0, 2.0, 1.0, 3.0), (1.0, 2.0));
        assert_eq!(interval_from_point_pair(3.0, -1.0), (-1.0, 3.0));
        assert_eq!(interval_clamp_point(0.0, 2.0, 5.0), 2.0);
        assert_eq!(interval_gap(0.0, 2.0, 5.0, 8.0), 3.0);
        assert_eq!(interval_gap(0.0, 2.0, 2.0, 8.0), 0.0);
        assert_eq!(interval_union(0.0, 2.0, 1.0, 3.0), (0.0, 3.0));
        assert_eq!(interval_union(1.0, 0.0, 1.0, 3.0), (1.0, 3.0));
        assert_eq!(interval_union(0.0, 2.0, 1.0, 0.0), (0.0, 2.0));
    }

    #[test]
    fn kernels_run_for_f32() {
        assert!(interval_contains(0.0_f32, 2.0, 1.0));
        assert_eq!(interval_gap(0.0_f32, 2.0, 5.0, 8.0), 3.0);
        assert_eq!(interval_union(1.0_f32, 0.0, 1.0, 3.0), (1.0, 3.0));
    }
}
