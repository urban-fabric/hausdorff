//! Latitude/longitude bounding-box kernels over scalar columns.

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

/// Mask that is true where `min <= max` fails on either axis.
#[inline]
#[must_use]
pub fn bbox_is_empty<T: Real>(min_lat: T, min_lon: T, max_lat: T, max_lon: T) -> Boolean<T> {
    mask_or(
        mask_not(min_lat.le(&max_lat)),
        mask_not(min_lon.le(&max_lon)),
    )
}

/// Mask that is true where `lat, lon` is finite and lies inside the inclusive
/// latitude/longitude bounding box.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn bbox_contains_latlon<T: Real>(
    min_lat: T,
    min_lon: T,
    max_lat: T,
    max_lon: T,
    lat: T,
    lon: T,
) -> Boolean<T> {
    let finite = mask_and(lat.is_finite(), lon.is_finite());
    let in_lat = mask_and(min_lat.le(&lat), lat.le(&max_lat));
    let in_lon = mask_and(min_lon.le(&lon), lon.le(&max_lon));

    mask_and(finite, mask_and(in_lat, in_lon))
}

/// Mask that is true where two inclusive latitude/longitude bounding boxes
/// share at least one point.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn bboxes_intersect<T: Real>(
    a_min_lat: T,
    a_min_lon: T,
    a_max_lat: T,
    a_max_lon: T,
    b_min_lat: T,
    b_min_lon: T,
    b_max_lat: T,
    b_max_lon: T,
) -> Boolean<T> {
    let a_nonempty = mask_not(bbox_is_empty(
        a_min_lat.clone(),
        a_min_lon.clone(),
        a_max_lat.clone(),
        a_max_lon.clone(),
    ));
    let b_nonempty = mask_not(bbox_is_empty(
        b_min_lat.clone(),
        b_min_lon.clone(),
        b_max_lat.clone(),
        b_max_lon.clone(),
    ));
    let lat_overlap = mask_and(a_min_lat.le(&b_max_lat), b_min_lat.le(&a_max_lat));
    let lon_overlap = mask_and(a_min_lon.le(&b_max_lon), b_min_lon.le(&a_max_lon));

    mask_and(
        mask_and(a_nonempty, b_nonempty),
        mask_and(lat_overlap, lon_overlap),
    )
}

/// Mask that is true where outer bounding box `a` wholly contains inner
/// bounding box `b`. Empty inner boxes are contained in everything.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn bbox_contains_bbox<T: Real>(
    a_min_lat: T,
    a_min_lon: T,
    a_max_lat: T,
    a_max_lon: T,
    b_min_lat: T,
    b_min_lon: T,
    b_max_lat: T,
    b_max_lon: T,
) -> Boolean<T> {
    let b_empty = bbox_is_empty(
        b_min_lat.clone(),
        b_min_lon.clone(),
        b_max_lat.clone(),
        b_max_lon.clone(),
    );
    let lat_contains = mask_and(a_min_lat.le(&b_min_lat), b_max_lat.le(&a_max_lat));
    let lon_contains = mask_and(a_min_lon.le(&b_min_lon), b_max_lon.le(&a_max_lon));

    mask_or(b_empty, mask_and(lat_contains, lon_contains))
}

/// Center latitude/longitude of a bounding box.
#[inline]
#[must_use]
pub fn bbox_center_latlon<T: Real>(min_lat: T, min_lon: T, max_lat: T, max_lon: T) -> (T, T) {
    (
        (min_lat + max_lat) * T::half(),
        (min_lon + max_lon) * T::half(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bbox_contains_latlon_is_mask_valued() {
        assert!(bbox_contains_latlon::<f64>(
            10.0, 20.0, 15.0, 25.0, 12.0, 21.0
        ));
        assert!(bbox_contains_latlon::<f64>(
            10.0, 20.0, 15.0, 25.0, 10.0, 25.0
        ));
        assert!(!bbox_contains_latlon::<f64>(
            10.0, 20.0, 15.0, 25.0, 16.0, 21.0
        ));
        assert!(!bbox_contains_latlon::<f64>(
            10.0,
            20.0,
            15.0,
            25.0,
            f64::NAN,
            21.0
        ));
    }

    #[test]
    fn bbox_pair_predicates_match_closed_rectangles() {
        assert!(!bbox_is_empty::<f64>(10.0, 20.0, 15.0, 25.0));
        assert!(bbox_is_empty::<f64>(15.0, 25.0, 10.0, 20.0));
        assert!(bboxes_intersect::<f64>(
            10.0, 20.0, 15.0, 25.0, 15.0, 25.0, 20.0, 30.0
        ));
        assert!(!bboxes_intersect::<f64>(
            10.0, 20.0, 15.0, 25.0, 16.0, 25.0, 20.0, 30.0
        ));
        assert!(bbox_contains_bbox::<f64>(
            10.0, 20.0, 15.0, 25.0, 11.0, 21.0, 12.0, 22.0
        ));
        assert!(bbox_contains_bbox::<f64>(
            10.0, 20.0, 15.0, 25.0, 12.0, 22.0, 11.0, 21.0
        ));
        assert!(!bbox_contains_bbox::<f64>(
            10.0, 20.0, 15.0, 25.0, 9.0, 21.0, 12.0, 22.0
        ));
    }

    #[test]
    fn bbox_value_kernels_run_for_f32() {
        assert_eq!(bbox_center_latlon(10.0_f32, 20.0, 16.0, 24.0), (13.0, 22.0));
        assert!(bbox_contains_latlon::<f32>(
            10.0, 20.0, 15.0, 25.0, 12.0, 21.0
        ));
    }
}
