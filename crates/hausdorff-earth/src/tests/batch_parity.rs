//! Bulk kernels against their scalar peers, elementwise.
//!
//! `LonLatBatch` exists to be fast, which means it is a second implementation of
//! things that already work — and a second implementation is a chance to disagree.
//! The failure mode is specific and nasty: the scalar path is what the tests and
//! the API docs describe, and the bulk path is what production actually runs. A
//! divergence is therefore invisible to exactly the tests most likely to exist.
//!
//! So the oracle here is **the scalar function itself**, applied pointwise. That
//! is not a weak oracle: the scalar geodesics are independently pinned against
//! GeographicLib identities in `hausdorff-base`, so agreement with them inherits
//! that validation rather than asserting the two paths are merely self-consistent.
//!
//! Agreement is required to be **bit-exact**, not approximate. The bulk kernels
//! are supposed to call the same formula over a column, not a faster
//! approximation of it, so any difference at all means the two paths have diverged
//! in a way worth knowing about. The one place that would be too strict is a
//! deliberately different summation order, and there is none here.
//!
//! The batch is also where the *degenerate* inputs matter most, because a column
//! is built from whatever a file contained: empty batches, single points, `NaN`
//! rows, out-of-range latitudes, and antimeridian spans all have to behave the
//! same as they would one at a time.

use crate::earth::{BBox, GeoPoint, LonLatBatch};
use hausdorff_base::earth::haversine_m;

fn g(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

/// A spread of points chosen to cover the seams rather than the easy middle:
/// both poles, both sides of the antimeridian, the equator, and the prime
/// meridian.
fn spread() -> Vec<GeoPoint> {
    vec![
        g(0.0, 0.0),
        g(90.0, 0.0),
        g(-90.0, 0.0),
        g(0.0, 180.0),
        g(0.0, -180.0),
        g(10.0, 179.999),
        g(10.0, -179.999),
        g(-45.5, 73.25),
        g(89.9, -179.9),
        g(-89.9, 179.9),
        g(51.5074, -0.1278),
        g(-33.8688, 151.2093),
        g(35.6762, 139.6503),
        g(-0.0, -0.0),
    ]
}

/// Targets to measure against, including ones coincident with batch members.
fn targets() -> Vec<GeoPoint> {
    vec![
        g(0.0, 0.0),
        g(90.0, 0.0),
        g(-90.0, 180.0),
        g(45.0, -90.0),
        g(10.0, 179.999),
        g(-0.0, 0.0),
    ]
}

#[test]
fn haversine_bulk_matches_scalar_bit_for_bit() {
    let points = spread();
    let batch = LonLatBatch::from_points(&points);
    for target in targets() {
        let bulk = batch.haversine_to_m(target);
        assert_eq!(bulk.len(), points.len(), "one output per input row");
        for (i, p) in points.iter().enumerate() {
            let scalar = haversine_m(*p, target);
            assert_eq!(
                bulk[i].to_bits(),
                scalar.to_bits(),
                "row {i} ({p:?}) to {target:?}: bulk {} vs scalar {}",
                bulk[i],
                scalar
            );
        }
        // And the unit-typed form must be the same number, just wrapped.
        let typed = batch.haversine_to(target);
        for (i, t) in typed.iter().enumerate() {
            assert_eq!(
                t.value().to_bits(),
                bulk[i].to_bits(),
                "row {i}: the typed wrapper must not change the value"
            );
        }
    }
}

#[test]
fn within_radius_agrees_with_comparing_the_scalar_distance() {
    // The radius predicate is allowed to be a cheaper formula than a full distance
    // followed by a comparison — that is the point of having it. But cheaper must
    // not mean different: the verdict has to match for every row, including rows
    // exactly on the boundary.
    let points = spread();
    let batch = LonLatBatch::from_points(&points);
    for target in targets() {
        let distances = batch.haversine_to_m(target);
        for radius in [
            0.0,
            1.0,
            1_000.0,
            1_000_000.0,
            10_000_000.0,
            20_100_000.0,
            f64::INFINITY,
        ] {
            let mask = batch.within_radius_m(target, radius);
            assert_eq!(mask.len(), points.len());
            for (i, inside) in mask.iter().enumerate() {
                let expected = distances[i] <= radius;
                assert_eq!(
                    *inside, expected,
                    "row {i} ({:?}) at radius {radius}: distance {} says {expected}, mask says {inside}",
                    points[i], distances[i]
                );
            }
        }
        // A negative radius can contain nothing, including a coincident point.
        for inside in batch.within_radius_m(target, -1.0) {
            assert!(!inside, "a negative radius must contain nothing");
        }
    }
}

#[test]
fn bbox_mask_agrees_with_scalar_containment() {
    let points = spread();
    let batch = LonLatBatch::from_points(&points);
    let boxes = [
        BBox::new(-10.0, -10.0, 10.0, 10.0),
        BBox::new(-90.0, -180.0, 90.0, 180.0),
        BBox::new(89.0, -180.0, 90.0, 180.0),
        BBox::new(0.0, 0.0, 0.0, 0.0),
        BBox::new(-1.0, 179.0, 11.0, 180.0),
    ];
    for bbox in boxes {
        let mask = batch.bbox_contains_mask(bbox);
        assert_eq!(mask.len(), points.len());
        for (i, inside) in mask.iter().enumerate() {
            assert_eq!(
                *inside,
                bbox.contains(points[i]),
                "row {i} ({:?}) against {bbox:?}",
                points[i]
            );
        }
    }
}

#[test]
fn h3_bulk_encoding_matches_one_at_a_time() {
    use hausdorff_index::index::h3::H3Cell;
    let points = spread();
    let batch = LonLatBatch::from_points(&points);
    for resolution in [0u8, 1, 7, 15] {
        let bulk = batch.h3_cells(resolution);
        assert_eq!(bulk.len(), points.len());
        for (i, p) in points.iter().enumerate() {
            let scalar = H3Cell::for_point(p.lat, p.lon, resolution);
            assert_eq!(
                bulk[i], scalar,
                "row {i} ({p:?}) at resolution {resolution}"
            );
        }
    }
    // An out-of-range resolution must fail for every row rather than panicking or
    // producing a bogus cell.
    for resolution in [16u8, 32, 255] {
        for cell in batch.h3_cells(resolution) {
            assert!(cell.is_none(), "resolution {resolution} is not valid");
        }
    }
}

#[test]
fn validity_mask_matches_the_scalar_range_checks() {
    let points = vec![
        g(0.0, 0.0),
        g(90.0, 180.0),
        g(-90.0, -180.0),
        g(91.0, 0.0),
        g(-91.0, 0.0),
        g(0.0, 181.0),
        g(0.0, -181.0),
        g(f64::NAN, 0.0),
        g(0.0, f64::NAN),
        g(f64::INFINITY, 0.0),
    ];
    let batch = LonLatBatch::from_points(&points);
    let mask = batch.valid_lonlat_mask();
    assert_eq!(mask.len(), points.len());
    for (i, valid) in mask.iter().enumerate() {
        let p = points[i];
        let expected = (-90.0..=90.0).contains(&p.lat) && (-180.0..=180.0).contains(&p.lon);
        assert_eq!(
            *valid, expected,
            "row {i} ({p:?}): lat and lon in range is {expected}"
        );
    }
}

#[test]
fn an_empty_batch_produces_empty_outputs_and_no_bbox() {
    let batch = LonLatBatch::from_points(&[]);
    assert!(batch.is_empty());
    assert_eq!(batch.len(), 0);
    assert!(batch.haversine_to_m(g(0.0, 0.0)).is_empty());
    assert!(batch.within_radius_m(g(0.0, 0.0), 1.0).is_empty());
    assert!(batch.valid_lonlat_mask().is_empty());
    assert!(batch.h3_cells(7).is_empty());
    assert!(
        batch.bounding_box().is_none(),
        "an empty batch has no bounding box, rather than a degenerate one at the origin"
    );
}

#[test]
fn a_single_point_batch_has_a_degenerate_bbox_containing_it() {
    for p in spread() {
        let batch = LonLatBatch::from_points(std::slice::from_ref(&p));
        let bbox = batch.bounding_box().expect("one point has a bounding box");
        assert!(
            bbox.contains(p),
            "the bounding box of {p:?} must contain it, got {bbox:?}"
        );
        assert_eq!(
            batch.haversine_to_m(p)[0],
            0.0,
            "distance to itself is zero"
        );
    }
}

#[test]
fn bounding_box_contains_every_row() {
    // The invariant that makes a bbox usable as a prune: no false negatives.
    let points = spread();
    let batch = LonLatBatch::from_points(&points);
    let bbox = batch.bounding_box().expect("non-empty");
    for (i, p) in points.iter().enumerate() {
        assert!(
            bbox.contains(*p),
            "row {i} ({p:?}) is outside its own batch's bounding box {bbox:?}"
        );
    }
    // Every row is inside, so the containment mask must be all true — the bulk
    // path and the bbox must agree about the batch that produced it.
    for (i, inside) in batch.bbox_contains_mask(bbox).iter().enumerate() {
        assert!(*inside, "row {i} excluded by its own batch's bounding box");
    }
}

#[test]
fn nan_rows_do_not_contaminate_their_neighbours() {
    // A columnar kernel that accumulated across rows — or vectorized with a
    // horizontal reduction — could let one bad row poison the ones beside it. The
    // clean rows must be unaffected.
    let points = vec![
        g(0.0, 0.0),
        g(f64::NAN, f64::NAN),
        g(10.0, 10.0),
        g(0.0, f64::NAN),
        g(-20.0, 30.0),
    ];
    let batch = LonLatBatch::from_points(&points);
    let target = g(5.0, 5.0);
    let bulk = batch.haversine_to_m(target);
    for (i, p) in points.iter().enumerate() {
        let scalar = haversine_m(*p, target);
        if scalar.is_nan() {
            assert!(bulk[i].is_nan(), "row {i} should stay NaN");
        } else {
            assert_eq!(
                bulk[i].to_bits(),
                scalar.to_bits(),
                "row {i} ({p:?}) was contaminated by a neighbouring NaN row"
            );
        }
    }
    // A NaN row cannot be within any radius.
    let mask = batch.within_radius_m(target, 1e9);
    assert!(!mask[1], "a NaN row is not within any radius");
    assert!(!mask[3], "a half-NaN row is not within any radius");
    assert!(mask[0] && mask[2] && mask[4], "clean rows still qualify");
}

#[test]
fn from_columns_rejects_mismatched_lengths() {
    assert!(LonLatBatch::from_columns(vec![1.0, 2.0], vec![3.0]).is_none());
    assert!(LonLatBatch::from_columns(vec![], vec![1.0]).is_none());
    assert!(LonLatBatch::from_columns(vec![1.0], vec![]).is_none());
    assert!(LonLatBatch::from_columns(vec![], vec![]).is_some());
    // And the column order is (lons, lats) — easy to transpose, so pinned.
    let b = LonLatBatch::from_columns(vec![10.0], vec![20.0]).expect("valid");
    assert_eq!(b.point(0).lon, 10.0, "first column is longitude");
    assert_eq!(b.point(0).lat, 20.0, "second column is latitude");
}

#[test]
fn large_batches_stay_consistent_with_the_scalar_path() {
    // A size past any plausible SIMD lane count and not a multiple of one, so a
    // vectorized kernel's scalar tail is exercised rather than skipped.
    let mut points = Vec::new();
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    for _ in 0..1037 {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let lat = ((seed >> 11) as f64 / (1u64 << 53) as f64) * 180.0 - 90.0;
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let lon = ((seed >> 11) as f64 / (1u64 << 53) as f64) * 360.0 - 180.0;
        points.push(g(lat, lon));
    }
    let batch = LonLatBatch::from_points(&points);
    let target = g(48.8566, 2.3522);
    let bulk = batch.haversine_to_m(target);
    for (i, p) in points.iter().enumerate() {
        assert_eq!(
            bulk[i].to_bits(),
            haversine_m(*p, target).to_bits(),
            "row {i} of {} diverged",
            points.len()
        );
    }
}
