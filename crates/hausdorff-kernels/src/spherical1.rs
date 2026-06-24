//! One-dimensional spherical angle kernels over scalar columns.
//!
//! These are the backend-neutral peers of `s1::Angle` and `s1::Sector`
//! predicates. Angles are radians unless the function name says otherwise.

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

/// Normalize an angle into `(-pi, pi]`.
#[inline]
#[must_use]
pub fn normalize_radians<T: Real>(radians: T) -> T {
    let pi = T::pi();
    let tau = T::tau();
    let r = radians.rem_euclid(tau.clone());
    T::select(&r.gt(&pi), r.clone() - tau, r)
}

/// Short angular distance between two angles in radians, in `[0, pi]`.
#[inline]
#[must_use]
pub fn angular_distance_radians<T: Real>(a: T, b: T) -> T {
    (a - b).cos().clamp(-T::one(), T::one()).acos()
}

/// Mask that is true where a sector radius denotes the empty sector.
#[inline]
#[must_use]
pub fn sector_is_empty<T: Real>(radius: T) -> Boolean<T> {
    radius.lt(&T::zero())
}

/// Mask that is true where a sector radius denotes the whole circle.
#[inline]
#[must_use]
pub fn sector_is_full<T: Real>(radius: T) -> Boolean<T> {
    radius.ge(&T::pi())
}

/// Angular sector width, clamped to a full turn and zero for empty sectors.
#[inline]
#[must_use]
pub fn sector_width_radians<T: Real>(radius: T) -> T {
    let empty = sector_is_empty(radius.clone());
    let width = (T::two() * radius).min(T::tau());
    T::select(&empty, T::zero(), width)
}

/// Mask that is true where `angle` lies in the sector centered at `center` with
/// angular `radius`.
#[inline]
#[must_use]
pub fn sector_contains_radians<T: Real>(center: T, radius: T, angle: T) -> Boolean<T> {
    let empty = sector_is_empty(radius.clone());
    let full = sector_is_full(radius.clone());
    let inside = (angle - center).cos().ge(&radius.cos());

    mask_or(full, mask_and(mask_not(empty), inside))
}

/// Mask that is true where two angular sectors share at least one angle.
#[inline]
#[must_use]
pub fn sectors_intersect_radians<T: Real>(
    a_center: T,
    a_radius: T,
    b_center: T,
    b_radius: T,
) -> Boolean<T> {
    let a_empty = sector_is_empty(a_radius.clone());
    let b_empty = sector_is_empty(b_radius.clone());
    let nonempty = mask_and(mask_not(a_empty), mask_not(b_empty));
    let full = mask_or(
        sector_is_full(a_radius.clone()),
        sector_is_full(b_radius.clone()),
    );
    let close = angular_distance_radians(a_center, b_center).le(&(a_radius + b_radius));

    mask_and(nonempty, mask_or(full, close))
}

/// Mask that is true where outer sector `a` wholly contains inner sector `b`.
#[inline]
#[must_use]
pub fn sector_contains_sector_radians<T: Real>(
    a_center: T,
    a_radius: T,
    b_center: T,
    b_radius: T,
) -> Boolean<T> {
    let a_empty = sector_is_empty(a_radius.clone());
    let a_full = sector_is_full(a_radius.clone());
    let b_empty = sector_is_empty(b_radius.clone());
    let contains = a_radius.ge(&(angular_distance_radians(a_center, b_center) + b_radius));

    mask_or(
        a_full,
        mask_or(b_empty, mask_and(mask_not(a_empty), contains)),
    )
}

/// Complementary sector `(center + pi, pi - radius)`.
#[inline]
#[must_use]
pub fn sector_complement_radians<T: Real>(center: T, radius: T) -> (T, T) {
    (center + T::pi(), T::pi() - radius)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn angle_normalization_matches_canonical_range() {
        assert!((normalize_radians(3.0 * std::f64::consts::PI) - std::f64::consts::PI).abs() < EPS);
        assert!((normalize_radians(-std::f64::consts::PI) - std::f64::consts::PI).abs() < EPS);
        assert!(
            (normalize_radians(1.5 * std::f64::consts::PI) + std::f64::consts::FRAC_PI_2).abs()
                < EPS
        );
    }

    #[test]
    fn angular_distance_is_short_arc() {
        let d = angular_distance_radians(0.0, 1.5 * std::f64::consts::PI);
        assert!((d - std::f64::consts::FRAC_PI_2).abs() < EPS);
    }

    #[test]
    fn sector_masks_match_s1_semantics() {
        let half = std::f64::consts::FRAC_PI_2;
        assert!(sector_is_empty::<f64>(-1.0));
        assert!(sector_is_full::<f64>(std::f64::consts::PI));
        assert_eq!(sector_width_radians(-1.0), 0.0);
        assert_eq!(sector_width_radians(4.0), std::f64::consts::TAU);

        assert!(sector_contains_radians::<f64>(0.0, half, half));
        assert!(!sector_contains_radians::<f64>(
            0.0,
            half,
            std::f64::consts::PI
        ));
        assert!(!sector_contains_radians::<f64>(0.0, -1.0, 0.0));
        assert!(sector_contains_radians::<f64>(
            0.0,
            std::f64::consts::PI,
            2.0
        ));

        assert!(sectors_intersect_radians::<f64>(0.0, 1.0, 1.5, 0.6));
        assert!(!sectors_intersect_radians::<f64>(0.0, 0.2, 1.0, 0.2));
        assert!(sector_contains_sector_radians::<f64>(0.0, 2.0, 0.5, 0.5));
        assert!(!sector_contains_sector_radians::<f64>(0.0, 0.5, 1.0, 0.5));
    }

    #[test]
    fn kernels_run_for_f32() {
        let d = angular_distance_radians(0.0_f32, std::f32::consts::PI);
        assert!((d - std::f32::consts::PI).abs() < 1e-6);
        assert!(sector_contains_radians(
            0.0_f32,
            std::f32::consts::FRAC_PI_2,
            0.25
        ));
        assert!(sectors_intersect_radians(0.0_f32, 1.0, 1.5, 0.6));
    }
}
