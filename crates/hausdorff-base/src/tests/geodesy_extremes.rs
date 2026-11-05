//! Geodesy where it is hardest: antipodes, poles, degenerate pairs, and invalid
//! input.
//!
//! The Karney solver exists *because* of these cases — Vincenty's method fails to
//! converge on nearly-antipodal pairs, and that convergence failure is the whole
//! reason to carry a more complicated algorithm. So the cases below are not
//! paranoia about unlikely inputs; they are the ones the algorithm was chosen for,
//! and they are where a port is most likely to have dropped something.
//!
//! **The oracles here are identities, not transcribed digits.** Three of them:
//!
//! - **Symmetry** — `inverse(a, b)` and `inverse(b, a)` must report the same
//!   distance. Nothing in the solver enforces this; the two calls take different
//!   branches through the azimuth seeding.
//! - **Round trip** — walking the reported distance along the reported initial
//!   bearing must land back on `b`. This ties `inverse` and `direct` together, and
//!   neither can be wrong alone without the other agreeing to be wrong the same
//!   way.
//! - **Independent closed form** — the pole-to-pole distance is half the meridian
//!   circumference, and the pole-to-equator distance is a quarter of it, so those
//!   two must be exactly 2:1 regardless of what either number is.
//!
//! An identity cannot be transcribed wrongly, which is what makes it worth more
//! here than a reference value copied out of a table.
//!
//! One result worth not "fixing": antipodal points **on the equator** are
//! 20 003 931 m apart, not the 20 037 508 m of half the equatorial circumference,
//! and the initial bearing is due north. On an oblate ellipsoid the meridian is
//! shorter than the equator, so the shortest path between equatorial antipodes
//! goes over the pole. GeographicLib reports the same, and a test that "corrected"
//! this to the equatorial figure would be asserting that the Earth is a sphere.

use crate::earth::{Bearing, Ellipsoid, GeoPoint, destination_m, haversine_m};

fn g(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_lat_lon(lat, lon)
}

/// Every pair that has a defined answer, paired with a label for failure messages.
fn hard_pairs() -> Vec<(&'static str, GeoPoint, GeoPoint)> {
    vec![
        ("identical points", g(45.0, 10.0), g(45.0, 10.0)),
        ("antipodes on the equator", g(0.0, 0.0), g(0.0, 180.0)),
        ("antipodes off the equator", g(30.0, 20.0), g(-30.0, -160.0)),
        ("nearly antipodal", g(0.0, 0.0), g(0.5, 179.7)),
        (
            "nearly antipodal, tighter",
            g(0.0, 0.0),
            g(0.0, 179.999_999_9),
        ),
        ("nearly antipodal, high lat", g(89.9, 0.0), g(-89.9, 180.0)),
        ("pole to pole", g(90.0, 0.0), g(-90.0, 0.0)),
        ("pole to pole across lon", g(90.0, 37.0), g(-90.0, -122.0)),
        ("north pole to equator", g(90.0, 0.0), g(0.0, 0.0)),
        ("both at the north pole", g(90.0, 0.0), g(90.0, 170.0)),
        ("both at the south pole", g(-90.0, 15.0), g(-90.0, -15.0)),
        ("one metre apart", g(45.0, 10.0), g(45.000_009, 10.0)),
        (
            "sub-micrometre apart",
            g(45.0, 10.0),
            g(45.0, 10.000_000_000_001),
        ),
        ("equatorial, short", g(0.0, 0.0), g(0.0, 0.001)),
        ("equatorial, long", g(0.0, -70.0), g(0.0, 70.0)),
        ("meridional", g(-40.0, 5.0), g(40.0, 5.0)),
        (
            "straddling the antimeridian",
            g(10.0, 179.9),
            g(10.0, -179.9),
        ),
        ("across the pole", g(89.0, 0.0), g(89.0, 180.0)),
    ]
}

#[test]
fn inverse_is_symmetric_on_every_hard_pair() {
    let e = Ellipsoid::WGS84;
    for (name, a, b) in hard_pairs() {
        let fwd = e.inverse(a, b).distance_m;
        let rev = e.inverse(b, a).distance_m;
        assert!(
            fwd.is_finite(),
            "{name}: distance must be finite, got {fwd}"
        );
        assert!(
            fwd >= 0.0,
            "{name}: distance must be non-negative, got {fwd}"
        );
        // Exact equality is a stronger claim than it looks: the two directions seed
        // the azimuth search differently. Allowing a nanometre keeps it honest
        // without letting a real asymmetry through.
        assert!(
            (fwd - rev).abs() <= 1e-9,
            "{name}: asymmetric, {fwd} vs {rev}"
        );
    }
}

#[test]
fn inverse_and_direct_round_trip_on_every_hard_pair() {
    let e = Ellipsoid::WGS84;
    for (name, a, b) in hard_pairs() {
        let inv = e.inverse(a, b);
        let landed = e.destination_m(a, inv.initial_bearing, inv.distance_m);
        let miss = e.inverse(landed, b).distance_m;
        // A micrometre over distances up to 20 000 km is ~5e-14 relative, i.e. at
        // the f64 floor. At the poles the bearing is genuinely under-determined, so
        // the landing point can differ in longitude while being the same place;
        // measuring the *miss distance* rather than comparing coordinates is what
        // makes the assertion well-posed there.
        assert!(
            miss <= 1e-6,
            "{name}: round trip missed by {miss} m (d = {}, azi1 = {})",
            inv.distance_m,
            inv.initial_bearing.degrees()
        );
    }
}

#[test]
fn degenerate_pairs_have_zero_distance() {
    let e = Ellipsoid::WGS84;
    for (name, a, b) in [
        ("identical", g(45.0, 10.0), g(45.0, 10.0)),
        (
            "same north pole, different lon",
            g(90.0, 0.0),
            g(90.0, 170.0),
        ),
        (
            "same south pole, different lon",
            g(-90.0, 15.0),
            g(-90.0, -15.0),
        ),
    ] {
        let d = e.distance_m(a, b);
        assert!(d < 1e-9, "{name}: expected ~0, got {d}");
        assert!(
            haversine_m(a, b) < 1e-9,
            "{name}: haversine should agree at zero"
        );
    }
}

#[test]
fn pole_distances_match_the_meridian_closed_form() {
    let e = Ellipsoid::WGS84;
    let half = e.distance_m(g(90.0, 0.0), g(-90.0, 0.0));
    let quarter = e.distance_m(g(90.0, 0.0), g(0.0, 0.0));
    // Independent of the absolute values: a half meridian is exactly two quarters.
    assert!(
        (half - 2.0 * quarter).abs() < 1e-6,
        "half meridian {half} should be twice the quarter {quarter}"
    );
    // And the equatorial antipodal distance must be the polar route, which on an
    // oblate ellipsoid is *shorter* than half the equatorial circumference.
    let equatorial_antipodes = e.distance_m(g(0.0, 0.0), g(0.0, 180.0));
    let half_equator = std::f64::consts::PI * 6_378_137.0;
    assert!(
        equatorial_antipodes < half_equator,
        "on an oblate ellipsoid the polar route ({equatorial_antipodes}) must beat \
         the equatorial one ({half_equator})"
    );
    assert!(
        (equatorial_antipodes - half).abs() < 1.0,
        "equatorial antipodes should travel a half meridian: {equatorial_antipodes} vs {half}"
    );
}

#[test]
fn geodesic_distance_grows_monotonically_along_a_meridian() {
    // A sweep rather than a point check: the inverse solver must not have a
    // non-monotone patch where two branches of the azimuth search meet.
    let e = Ellipsoid::WGS84;
    let origin = g(0.0, 0.0);
    let mut prev = -1.0;
    let mut i = 0;
    while i <= 1800 {
        let lat = f64::from(i) / 20.0 - 90.0;
        let d = e.distance_m(origin, g(lat, 0.0));
        if lat <= 0.0 {
            // Southbound: distance shrinks toward the origin.
            assert!(d.is_finite(), "distance to lat {lat} is not finite");
        } else {
            assert!(
                d >= prev - 1e-6,
                "distance is not monotone crossing lat {lat}: {prev} then {d}"
            );
        }
        if lat >= 0.0 {
            prev = d;
        }
        i += 1;
    }
}

#[test]
fn longitude_is_cyclic_but_latitude_is_not() {
    let e = Ellipsoid::WGS84;
    // Longitude wraps: any value is a real angle, so adding whole turns changes
    // nothing. This is the convention, not an accident, so it is pinned.
    let base = e.distance_m(g(10.0, 20.0), g(-5.0, 100.0));
    for turns in [1.0, -1.0, 3.0, -1000.0] {
        let shifted = e.distance_m(g(10.0, 20.0 + 360.0 * turns), g(-5.0, 100.0));
        assert!(
            (shifted - base).abs() < 1e-6,
            "longitude should be cyclic: {turns} turns moved {base} to {shifted}"
        );
    }
    // Latitude does not wrap. GeographicLib's convention is NaN for out of range,
    // and the ellipsoidal path follows it.
    for bad_lat in [91.0, -91.0, 180.0, 1e9] {
        let d = e.distance_m(g(bad_lat, 0.0), g(0.0, 0.0));
        assert!(
            d.is_nan(),
            "latitude {bad_lat} is out of range; the ellipsoidal path should be NaN, got {d}"
        );
    }
}

#[test]
fn haversine_does_not_validate_latitude_but_karney_does() {
    // The two metrics genuinely differ here, and the difference is documented on
    // `haversine_m`. Pinned so that nobody "fixes" one to match the other without
    // deciding which convention the library wants.
    let e = Ellipsoid::WGS84;
    let out_of_range = g(91.0, 0.0);
    let origin = g(0.0, 0.0);
    assert!(
        e.distance_m(out_of_range, origin).is_nan(),
        "the ellipsoidal path signals an out-of-range latitude"
    );
    assert!(
        haversine_m(out_of_range, origin).is_finite(),
        "the spherical path deliberately does not check, and returns a finite number"
    );
}

#[test]
fn nan_coordinates_propagate_rather_than_panicking() {
    let e = Ellipsoid::WGS84;
    for (name, a, b) in [
        ("NaN lat", g(f64::NAN, 0.0), g(0.0, 0.0)),
        ("NaN lon", g(0.0, f64::NAN), g(0.0, 0.0)),
        ("inf lat", g(f64::INFINITY, 0.0), g(0.0, 0.0)),
        ("inf lon", g(0.0, f64::INFINITY), g(0.0, 0.0)),
        (
            "NaN both sides",
            g(f64::NAN, f64::NAN),
            g(f64::NAN, f64::NAN),
        ),
    ] {
        assert!(
            e.distance_m(a, b).is_nan(),
            "{name}: a non-finite coordinate should produce NaN, not a number"
        );
        assert!(haversine_m(a, b).is_nan(), "{name}: haversine likewise");
    }
}

#[test]
fn direct_survives_hostile_distances_and_bearings() {
    let e = Ellipsoid::WGS84;
    let origin = g(45.0, 10.0);
    // A non-finite distance has no defined endpoint, so NaN is the answer; the
    // requirement is that it does not panic or return a plausible-looking point.
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let p = e.destination_m(origin, Bearing::from_degrees(45.0), bad);
        assert!(
            p.lat.is_nan() && p.lon.is_nan(),
            "distance {bad} should give NaN, got {p:?}"
        );
    }
    // A negative distance is meaningful: it walks backwards along the geodesic, so
    // it must land diametrically opposite the forward walk.
    let fwd = e.destination_m(origin, Bearing::from_degrees(45.0), 1_000_000.0);
    let back = e.destination_m(origin, Bearing::from_degrees(45.0), -1_000_000.0);
    let rev = e.destination_m(origin, Bearing::from_degrees(225.0), 1_000_000.0);
    assert!(
        e.distance_m(back, rev) < 1e-6,
        "walking -d at 45 deg should equal walking +d at 225 deg: {back:?} vs {rev:?}"
    );
    assert!(e.distance_m(origin, fwd).is_finite());
    // Bearings are angles, so they wrap; a huge or negative bearing is valid.
    for turns in [1.0, -1.0, 2.0] {
        let p = e.destination_m(origin, Bearing::from_degrees(45.0 + 360.0 * turns), 5_000.0);
        assert!(
            e.distance_m(
                p,
                e.destination_m(origin, Bearing::from_degrees(45.0), 5_000.0)
            ) < 1e-6,
            "bearing should be cyclic over {turns} turns"
        );
    }
    // More than a full circumference is still a defined walk, just a wrapped one.
    let far = e.destination_m(origin, Bearing::from_degrees(90.0), 1e12);
    assert!(
        far.lat.is_finite() && far.lon.is_finite(),
        "a very long walk should still land somewhere real, got {far:?}"
    );
    assert!(
        (-90.0..=90.0).contains(&far.lat) && (-180.0..=180.0).contains(&far.lon),
        "a wrapped walk must still produce in-range coordinates, got {far:?}"
    );
}

#[test]
fn spherical_destination_stays_in_range_everywhere() {
    // The spherical `destination_m` clamps latitude into the asin domain and
    // normalizes longitude; sweep bearings and distances to confirm no escape.
    let origin = g(89.0, 170.0);
    let mut bearing_i = 0;
    while bearing_i < 360 {
        for dist in [1.0, 1e5, 1e6, 1e7, 2e7, 4e7] {
            let p = destination_m(origin, Bearing::from_degrees(f64::from(bearing_i)), dist);
            assert!(
                (-90.0..=90.0).contains(&p.lat),
                "latitude escaped at bearing {bearing_i}, distance {dist}: {p:?}"
            );
            assert!(
                (-180.0..=180.0).contains(&p.lon),
                "longitude escaped at bearing {bearing_i}, distance {dist}: {p:?}"
            );
        }
        bearing_i += 7;
    }
}
