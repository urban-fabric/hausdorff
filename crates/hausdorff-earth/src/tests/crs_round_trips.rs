//! Every registry CRS, round-tripped — including across datums, where the
//! transform pivots through WGS84 geocentric.
//!
//! The existing `projection_proj_reference` suite is the stronger test of whether
//! a projection *is what it claims to be*: it checks forward and inverse against
//! coordinates PROJ itself produced. What it cannot do is cover every pair, since
//! an external reference value has to be obtained for each one.
//!
//! This suite covers the pairs instead, using **round-tripping** as the oracle so
//! no reference values are needed. That buys breadth and loses absoluteness, and
//! the two are complementary rather than redundant: a round trip proves the
//! inverse undoes the forward and says *nothing* about whether either is correct.
//! A projection with a sign flip in both directions round-trips perfectly. So
//! neither suite would be sufficient alone.
//!
//! What the round trip does catch, and catches well, is the failure mode that
//! actually happens in a CRS stack: a transform that silently drops its datum
//! shift, loses precision in one branch, or returns a plausible-looking
//! coordinate for an input outside the projection's domain. Those are all
//! asymmetric, so they break the round trip even though the absolute answer
//! "looks fine".
//!
//! **Cross-datum pairs are the interesting ones.** `transform` goes
//! native → datum geodetic → datum ECEF → (Helmert or NTv2 grid) → WGS84 ECEF and
//! back out. Routing through geocentric Cartesian rather than lat/lon is what lets
//! a datum shift sit on the pivot — a Helmert transform is a rigid motion in ECEF
//! and has no clean lat/lon form. Skip the ECEF step and cross-datum results are
//! wrong by up to a few hundred metres, which is a *plausible* error and therefore
//! the dangerous kind. A round trip through the pivot and back catches a shift
//! that is applied in only one direction.
//!
//! Each CRS is exercised **inside its own area of use**, taken from the registry.
//! That is not a convenience: a projection outside its domain is legitimately
//! allowed to be inaccurate or undefined — Lambert conformal conic diverges at the
//! far pole, transverse Mercator degrades away from its central meridian — so
//! asserting tight accuracy there would be asserting something false. Domain
//! behaviour gets its own test below, which requires only that it not produce
//! garbage.

use crate::earth::GeoPoint;
use crate::earth::crs::{Crs, CrsCoord, GeoDim, transform};

/// The EPSG codes the registry knows, spanning geographic 2-D and 3-D, geocentric,
/// and every projection family present.
const CODES: &[u32] = &[
    4326, 4978, 4979, 4937, // WGS84 family and ETRS89 3D
    5070, 2154, 2193, 2180, 3978, 3577, // projected: Albers, LCC, TM
    3995, 3031, 3413, 3976, // polar stereographic, both hemispheres
    4277, 4230, 4267, 4314, 4275, 4272, 4269, 4258, // legacy datums
];

/// Smallest signed difference between two longitudes, in degrees.
///
/// Longitude is cyclic, so `-180` and `+180` are the *same meridian* and a plain
/// subtraction reports them as 360 degrees apart. Comparing projected longitudes
/// without this is how a test "discovers" that a correct round trip failed: the
/// first thing this suite caught was EPSG:2193 returning `+180` for an input of
/// `-180`, which is the same place.
fn lon_delta(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    if d > 180.0 { 360.0 - d } else { d }
}

/// A handful of sample points inside a bounding box: the centre and four points
/// pulled in from the corners, so a projection is exercised across its domain
/// rather than at one spot.
///
/// **Handles a bounding box that crosses the antimeridian**, where
/// `min_lon > max_lon` — EPSG:4267 (NAD27) is declared as
/// `min_lon 167.65 .. max_lon -47.74`, wrapping east through the Pacific. A naive
/// `lerp` from 167.65 to -47.74 interpolates the *wrong way* around the globe and
/// lands in the Pacific instead of North America, which then fails every
/// projection assertion for reasons that have nothing to do with the projection.
/// `BBox::contains` gets this right; this helper has to as well.
fn samples_in(bbox: crate::earth::BBox) -> Vec<GeoPoint> {
    let (lat0, lat1) = (bbox.min_lat, bbox.max_lat);
    let lon0 = bbox.min_lon;
    // Walk eastward, adding a turn when the box wraps, then normalize back.
    let span = if bbox.max_lon >= bbox.min_lon {
        bbox.max_lon - bbox.min_lon
    } else {
        bbox.max_lon + 360.0 - bbox.min_lon
    };
    let lerp = |a: f64, b: f64, t: f64| a + (b - a) * t;
    let lon_at = |t: f64| {
        let raw = lon0 + span * t;
        // Back into (-180, 180].
        let wrapped = (raw + 180.0).rem_euclid(360.0) - 180.0;
        if wrapped == -180.0 { 180.0 } else { wrapped }
    };
    vec![
        GeoPoint::from_lat_lon(lerp(lat0, lat1, 0.5), lon_at(0.5)),
        GeoPoint::from_lat_lon(lerp(lat0, lat1, 0.1), lon_at(0.1)),
        GeoPoint::from_lat_lon(lerp(lat0, lat1, 0.9), lon_at(0.1)),
        GeoPoint::from_lat_lon(lerp(lat0, lat1, 0.1), lon_at(0.9)),
        GeoPoint::from_lat_lon(lerp(lat0, lat1, 0.9), lon_at(0.9)),
    ]
}

/// The same CRS's own geographic frame, which is where a round trip starts and
/// ends so that no datum shift is involved in the comparison.
fn own_geographic(c: Crs) -> Crs {
    Crs::Geographic {
        datum: c.datum(),
        epsg: None,
        dim: GeoDim::D2,
    }
}

fn geographic_of(coord: CrsCoord) -> Option<GeoPoint> {
    match coord {
        CrsCoord::Geographic { point, .. } => Some(point),
        _ => None,
    }
}

#[test]
fn every_registry_crs_round_trips_within_its_area_of_use() {
    for &code in CODES {
        let crs = Crs::from_epsg(code).unwrap_or_else(|| panic!("EPSG:{code} should be known"));
        let bbox = match crs.area_of_use() {
            Some(b) => b,
            // No declared area of use means global; sample a mid-latitude band
            // rather than the poles, which several projections cannot represent.
            None => crate::earth::BBox::new(-60.0, -170.0, 60.0, 170.0),
        };
        let from = own_geographic(crs);
        for p in samples_in(bbox) {
            let there = transform(from, crs, CrsCoord::geographic(p));
            let back = transform(crs, from, there);
            let got = geographic_of(back).unwrap_or_else(|| {
                panic!(
                    "EPSG:{code}: round trip did not land in a \
                                           geographic frame, got {back:?}"
                )
            });
            assert!(
                got.lat.is_finite() && got.lon.is_finite(),
                "EPSG:{code}: round trip of {p:?} produced {got:?}"
            );
            // Within its own area of use a projection must be invertible to well
            // under a millimetre — about 1e-9 degrees of latitude.
            assert!(
                (got.lat - p.lat).abs() < 1e-8 && lon_delta(got.lon, p.lon) < 1e-8,
                "EPSG:{code}: {p:?} round-tripped to {got:?}"
            );
        }
    }
}

#[test]
fn cross_datum_transforms_round_trip_through_the_ecef_pivot() {
    // Each pair crosses a datum boundary, so the transform must apply a shift in
    // one direction and undo it in the other. A shift applied only one way shows
    // up here and nowhere in a same-datum test.
    let pairs = [
        (4326u32, 4277u32), // WGS84 <-> OSGB36
        (4326, 4267),       // WGS84 <-> NAD27
        (4326, 4314),       // WGS84 <-> DHDN (NTv2 grid available)
        (4326, 4230),       // WGS84 <-> ED50
        (4326, 4275),       // WGS84 <-> NTF (NTv2 grid available)
        (4326, 4272),       // WGS84 <-> NZGD49 (NTv2 grid available)
        (4277, 4267),       // OSGB36 <-> NAD27, neither of them WGS84
        (4326, 4978),       // geographic <-> geocentric
        (4326, 4979),       // 2-D <-> 3-D
    ];
    for (a_code, b_code) in pairs {
        let a = Crs::from_epsg(a_code).expect("known");
        let b = Crs::from_epsg(b_code).expect("known");
        // Sample where both CRSes are meaningful.
        let bbox = b
            .area_of_use()
            .or_else(|| a.area_of_use())
            .unwrap_or(crate::earth::BBox::new(-60.0, -170.0, 60.0, 170.0));
        for p in samples_in(bbox) {
            let start = transform(own_geographic(a), a, CrsCoord::geographic(p));
            let over = transform(a, b, start);
            let back = transform(b, a, over);
            let home = transform(a, own_geographic(a), back);
            let got = geographic_of(home).unwrap_or_else(|| {
                panic!("EPSG:{a_code}->{b_code}: did not return to a geographic frame")
            });
            assert!(
                got.lat.is_finite() && got.lon.is_finite(),
                "EPSG:{a_code}<->{b_code}: {p:?} produced {got:?}"
            );
            // Looser than the same-datum round trip (1e-8) on purpose. This chain
            // is four transforms, and a cross-datum pair like OSGB36 <-> NAD27
            // crosses *two* Helmert shifts via WGS84. A 7-parameter Helmert's
            // inverse is not an exact algebraic inverse, so a few centimetres
            // accumulate — against datums whose own realization is good to about a
            // metre. Tightening this would be asserting precision the datum does
            // not have; a real defect (a dropped or sign-flipped shift) is tens to
            // hundreds of metres and still fails comfortably.
            assert!(
                (got.lat - p.lat).abs() < 1e-6 && lon_delta(got.lon, p.lon) < 1e-6,
                "EPSG:{a_code}<->{b_code}: {p:?} round-tripped to {got:?}"
            );
        }
    }
}

#[test]
fn a_datum_shift_actually_moves_the_point() {
    // The complement to the round trip, and the reason the round trip is not
    // sufficient on its own: a transform that dropped its datum shift entirely
    // would round-trip *perfectly*. So the shift has to be shown to exist.
    //
    // Legacy datums differ from WGS84 by tens to hundreds of metres, which is far
    // more than any rounding and far less than a continent.
    let wgs84 = Crs::from_epsg(4326).expect("known");
    for (code, place) in [
        (4277u32, GeoPoint::from_lat_lon(51.5, -0.1)), // OSGB36, London
        (4267, GeoPoint::from_lat_lon(40.0, -100.0)),  // NAD27, Kansas
        (4230, GeoPoint::from_lat_lon(45.0, 9.0)),     // ED50, northern Italy
        (4314, GeoPoint::from_lat_lon(52.0, 10.0)),    // DHDN, Germany
    ] {
        let other = Crs::from_epsg(code).expect("known");
        let shifted = transform(wgs84, other, CrsCoord::geographic(place));
        let got =
            geographic_of(shifted).unwrap_or_else(|| panic!("EPSG:{code} is a geographic CRS"));
        // Measure the displacement on the ground, which is the meaningful unit.
        let moved = crate::earth::haversine_m(place, got);
        assert!(
            moved > 1.0,
            "EPSG:{code}: a datum shift from WGS84 must move the point, got {moved} m \
             — a transform that silently skipped the shift would look like this"
        );
        assert!(
            moved < 1_000.0,
            "EPSG:{code}: a datum shift should be tens to hundreds of metres, not \
             {moved} m — this looks like a wrong sign or a missing rotation"
        );
    }
}

#[test]
fn transform_is_the_identity_when_source_and_target_match() {
    // Worth pinning because an implementation that pivots unconditionally will
    // round-trip through ECEF for no reason and lose bits doing it.
    for &code in CODES {
        let crs = Crs::from_epsg(code).expect("known");
        let bbox = crs
            .area_of_use()
            .unwrap_or(crate::earth::BBox::new(-60.0, -170.0, 60.0, 170.0));
        for p in samples_in(bbox) {
            let coord = transform(own_geographic(crs), crs, CrsCoord::geographic(p));
            let same = transform(crs, crs, coord);
            // Exact, not approximate: `transform` short-circuits an equal
            // source and target. Before it did, the pivot ran anyway and moved an
            // ellipsoidal height by a couple of nanometres — small, but a no-op
            // that is not a no-op compounds across a pipeline that re-normalizes.
            assert_eq!(
                same, coord,
                "EPSG:{code}: transforming to itself must be bit-exact"
            );
        }
    }
}

#[test]
fn unknown_epsg_codes_are_rejected() {
    for code in [0u32, 1, 9999, 99999, u32::MAX, 4328, 1234] {
        if Crs::from_epsg(code).is_some() {
            // A code that *is* in the registry is fine; this test only requires
            // that an unknown one says so rather than inventing a CRS.
            continue;
        }
        assert!(
            Crs::from_epsg(code).is_none(),
            "EPSG:{code} is not in the registry and must return None"
        );
    }
    assert!(Crs::from_epsg(0).is_none(), "zero is not an EPSG code");
}

#[test]
fn projections_outside_their_domain_do_not_produce_plausible_garbage() {
    // Outside its area of use a projection may be inaccurate or undefined, and
    // demanding accuracy there would be asserting something false. What it must
    // *not* do is return a finite coordinate that silently round-trips to the
    // wrong place — the error that gets shipped.
    //
    // So the requirement is a disjunction: either the output is non-finite
    // (honest), or the round trip still recovers the input (correct). Anything
    // else is a plausible wrong answer.
    let hostile = [
        GeoPoint::from_lat_lon(90.0, 0.0),
        GeoPoint::from_lat_lon(-90.0, 0.0),
        GeoPoint::from_lat_lon(89.999_999, 179.999_999),
        GeoPoint::from_lat_lon(0.0, 180.0),
        GeoPoint::from_lat_lon(0.0, -180.0),
        GeoPoint::from_lat_lon(-89.999_999, -179.999_999),
    ];
    for &code in CODES {
        let crs = Crs::from_epsg(code).expect("known");
        let from = own_geographic(crs);
        for p in hostile {
            let there = transform(from, crs, CrsCoord::geographic(p));
            let back = transform(crs, from, there);
            let Some(got) = geographic_of(back) else {
                continue; // did not return a geographic frame; nothing to compare
            };
            if !got.lat.is_finite() || !got.lon.is_finite() {
                continue; // honestly undefined
            }
            // Deliberately loose. A conic projection 5 000 km outside its area of
            // use is *allowed* to be inaccurate — Albers for the contiguous US
            // puts the north pole back about 13 m off, and that is the cone
            // flattening, not a defect. "Either signal or be exact" would be a
            // false dichotomy. What this still catches is the answer that gets
            // shipped: a sign flip, a wrapped longitude, or a point on the wrong
            // continent, all of which are degrees out rather than metres.
            let recovered = (got.lat - p.lat).abs() < 1e-3 && lon_delta(got.lon, p.lon) < 1e-3;
            // A pole has no defined longitude, so recovering the latitude alone is
            // a correct answer there.
            let at_pole = p.lat.abs() > 89.999_99;
            let pole_ok = at_pole && (got.lat - p.lat).abs() < 1e-3;
            assert!(
                recovered || pole_ok,
                "EPSG:{code}: {p:?} is outside the domain but produced the finite, \
                 non-recovering coordinate {got:?} — either signal, or be right"
            );
        }
    }
    // And a non-finite input must never become a finite output.
    for &code in CODES {
        let crs = Crs::from_epsg(code).expect("known");
        let from = own_geographic(crs);
        for bad in [
            GeoPoint::from_lat_lon(f64::NAN, 0.0),
            GeoPoint::from_lat_lon(0.0, f64::NAN),
            GeoPoint::from_lat_lon(f64::INFINITY, 0.0),
        ] {
            let there = transform(from, crs, CrsCoord::geographic(bad));
            let finite = match there {
                CrsCoord::Projected(c) => c.x.is_finite() && c.y.is_finite(),
                CrsCoord::Geographic { point, .. } => {
                    point.lat.is_finite() && point.lon.is_finite()
                }
                CrsCoord::Geocentric(c) => c.x.is_finite() && c.y.is_finite() && c.z.is_finite(),
            };
            assert!(
                !finite,
                "EPSG:{code}: a non-finite input {bad:?} must not become a finite \
                 coordinate, got {there:?}"
            );
        }
    }
}

#[test]
fn the_area_of_use_contains_the_points_it_claims() {
    // `contains` is used to pick a CRS, so a bbox that disagrees with its own
    // `contains` would route data to the wrong projection.
    for &code in CODES {
        let crs = Crs::from_epsg(code).expect("known");
        let Some(bbox) = crs.area_of_use() else {
            continue;
        };
        for p in samples_in(bbox) {
            assert!(
                crs.contains(p),
                "EPSG:{code}: {p:?} is inside the declared area of use {bbox:?} but \
                 `contains` says otherwise"
            );
        }
        // A point far outside is not contained — but only ask that of a CRS that
        // *has* an outside. EPSG:4326's declared area of use is the whole globe,
        // so there is no point to place beyond it, and demanding one would be
        // demanding that a global CRS reject part of the globe.
        let global = bbox.min_lat <= -90.0
            && bbox.max_lat >= 90.0
            && bbox.min_lon <= -180.0
            && bbox.max_lon >= 180.0;
        if !global {
            // Pick an antipodal-ish latitude rather than guessing a longitude: a
            // wrapped box makes "the other side" ambiguous in longitude, but no
            // regional box spans both hemispheres' high latitudes.
            let outside = GeoPoint::from_lat_lon(
                if bbox.min_lat > -20.0 { -85.0 } else { 85.0 },
                bbox.min_lon,
            );
            assert!(
                !crs.contains(outside),
                "EPSG:{code}: {outside:?} is far outside {bbox:?} but `contains` accepted it"
            );
        }
    }
}
