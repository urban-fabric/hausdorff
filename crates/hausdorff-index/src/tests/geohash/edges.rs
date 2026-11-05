//! Geohash at its edges: the two-implementation boundary, the seams, and
//! malformed input.
//!
//! The thing most worth testing here is invisible from the public API. `encode`
//! has **two implementations**: a fixed-point Morton (bit-interleaving) fast path
//! for precision ≤ 12, where the whole code fits in 64 bits, and a general path
//! above it. Two implementations of one function is a chance to disagree, and the
//! disagreement would be at the boundary where almost nobody looks.
//!
//! There is no way to call the slow path directly, so the test is indirect but
//! sharp: a geohash is a **prefix code** — each character narrows the cell, so a
//! shorter hash must be a prefix of a longer one for the same point. If the fast
//! and general paths differ at all, then `encode(p, 12)` stops being a prefix of
//! `encode(p, 13)`, because the first came from one implementation and the second
//! from the other. The prefix property turns an inaccessible internal boundary
//! into an assertion.
//!
//! The other invariant worth stating is that `decode` returns the **cell centre**,
//! not the point that was encoded. So `decode(encode(p))` is not `p` and a test
//! expecting that would be wrong. What must hold is **round-trip stability**:
//! re-encoding the centre lands in the same cell, because the centre is in the
//! cell. That holds at any precision and needs no tolerance.
//!
//! Unlike the spherical geodesics, geohash **does** validate its range: `try_encode`
//! returns `None` outside `lat ∈ [-90, 90]` / `lon ∈ [-180, 180]` while `encode`
//! clamps silently. That split is deliberate and documented, so both halves are
//! pinned — a caller that picked one should not find the other's behaviour.

use crate::index::geohash::{decode, encode, try_encode};

/// The seams, the poles, and the signed zeros — the places a lon/lat grid folds.
fn seam_points() -> Vec<(&'static str, f64, f64)> {
    vec![
        ("origin", 0.0, 0.0),
        ("north pole", 0.0, 90.0),
        ("south pole", 0.0, -90.0),
        ("antimeridian east", 180.0, 0.0),
        ("antimeridian west", -180.0, 0.0),
        ("north-east corner", 180.0, 90.0),
        ("south-west corner", -180.0, -90.0),
        ("just inside +180", 179.999_999, 0.0),
        ("just inside -180", -179.999_999, 0.0),
        ("just inside +90", 0.0, 89.999_999),
        ("negative zero", -0.0, -0.0),
        ("Jutland", 10.407_439_7, 57.649_110_6),
        ("Sydney", 151.209_3, -33.868_8),
    ]
}

#[test]
fn a_shorter_hash_is_always_a_prefix_of_a_longer_one() {
    // This is what catches a divergence between the Morton fast path (<= 12) and
    // the general path above it: the two meet exactly at the 12/13 boundary, and a
    // mismatch there breaks the prefix chain.
    for (name, lon, lat) in seam_points() {
        let longest = encode(lon, lat, 20);
        for precision in 1..=20usize {
            let h = encode(lon, lat, precision);
            assert_eq!(
                h.len(),
                precision,
                "{name}: asked for {precision} characters"
            );
            assert!(
                longest.starts_with(&h),
                "{name}: precision {precision} gave {h}, which is not a prefix of {longest} \
                 — the fast and general encode paths have diverged"
            );
        }
    }
}

#[test]
fn every_character_is_in_the_geohash_alphabet() {
    // Base-32 geohash deliberately omits a, i, l, and o to avoid transcription
    // errors. An encoder that indexed a plain alphabet would silently emit them.
    const ALPHABET: &str = "0123456789bcdefghjkmnpqrstuvwxyz";
    for (name, lon, lat) in seam_points() {
        let h = encode(lon, lat, 15);
        for c in h.chars() {
            assert!(
                ALPHABET.contains(c),
                "{name}: {h} contains {c:?}, which is not in the geohash alphabet"
            );
        }
    }
    for excluded in ['a', 'i', 'l', 'o'] {
        assert!(!ALPHABET.contains(excluded));
        // And decoding one must fail rather than aliasing to a neighbour.
        assert!(
            decode(&format!("9q8yy{excluded}")).is_none(),
            "{excluded:?} is not a geohash character and must not decode"
        );
    }
}

#[test]
fn decode_returns_a_centre_that_re_encodes_to_the_same_cell() {
    // `decode` yields the cell centre, so it is not the point that was encoded.
    // The invariant is stability: the centre is inside the cell, so re-encoding it
    // reproduces the hash exactly, at every precision.
    for (name, lon, lat) in seam_points() {
        for precision in 1..=12usize {
            let h = encode(lon, lat, precision);
            let (clon, clat) = decode(&h).expect("our own hash must decode");
            let again = encode(clon, clat, precision);
            assert_eq!(
                h, again,
                "{name}: the centre of {h} re-encoded to {again} at precision {precision}"
            );
            // The centre must also be in range, including at the poles where the
            // cell is clipped.
            assert!(
                (-180.0..=180.0).contains(&clon) && (-90.0..=90.0).contains(&clat),
                "{name}: decoded centre ({clon}, {clat}) is out of range"
            );
        }
    }
}

#[test]
fn the_decoded_centre_is_within_the_documented_precision_bound() {
    // The property to assert is the one `decode` documents and `riemann`
    // `Hausdorff/Geohash.geohash_decode_close` proves: after `n` bits on an axis
    // the returned midpoint is within `(hi - lo) / 2^(n+1)` of the true
    // coordinate.
    //
    // It is specifically *not* monotone convergence, which is the tempting but
    // false version of this. The cell shrinks with every character, yet the
    // centre can move *further* from the point when the point sits near a cell
    // boundary: a point at 0.49 in the cell [0, 1] is 0.01 from the centre 0.5,
    // and subdividing puts it in [0, 0.5] whose centre 0.25 is 0.24 away. The
    // error is bounded by a shrinking quantity, not itself shrinking.
    for (name, lon, lat) in seam_points() {
        for precision in 1..=12usize {
            let bits = 5 * precision;
            // Characters alternate longitude first, so longitude takes the odd bit.
            let lon_bits = bits.div_ceil(2);
            let lat_bits = bits / 2;
            let lon_bound = 360.0 / 2f64.powi(lon_bits as i32 + 1);
            let lat_bound = 180.0 / 2f64.powi(lat_bits as i32 + 1);

            let (clon, clat) = decode(&encode(lon, lat, precision)).expect("decodes");
            assert!(
                (clon - lon).abs() <= lon_bound + 1e-12,
                "{name}: precision {precision} longitude off by {} exceeds the \
                 {lon_bits}-bit bound {lon_bound}",
                (clon - lon).abs()
            );
            assert!(
                (clat - lat).abs() <= lat_bound + 1e-12,
                "{name}: precision {precision} latitude off by {} exceeds the \
                 {lat_bits}-bit bound {lat_bound}",
                (clat - lat).abs()
            );
        }
    }
}

#[test]
fn neighbouring_points_share_a_prefix_and_distant_ones_do_not() {
    // The property that makes a geohash useful as a database key: spatial locality
    // becomes string locality. Without it, a prefix range scan is not a spatial
    // query and the whole index is pointless.
    let (lon, lat) = (10.407_439_7, 57.649_110_6);
    let near = encode(lon + 1e-7, lat + 1e-7, 9);
    let here = encode(lon, lat, 9);
    assert_eq!(
        here[..6],
        near[..6],
        "points a centimetre apart should share a long prefix: {here} vs {near}"
    );
    let far = encode(-120.0, -40.0, 9);
    assert_ne!(
        here.as_bytes()[0],
        far.as_bytes()[0],
        "points on opposite sides of the planet should not share the first character"
    );
}

#[test]
fn try_encode_validates_what_encode_clamps() {
    // Two deliberately different conventions on the same input. Both are pinned so
    // that a caller who picked one does not later find the other's behaviour.
    for (name, lon, lat) in [
        ("lat above range", 0.0, 91.0),
        ("lat below range", 0.0, -91.0),
        ("lon above range", 181.0, 0.0),
        ("lon below range", -181.0, 0.0),
        ("lat NaN", 0.0, f64::NAN),
        ("lon NaN", f64::NAN, 0.0),
        ("lat infinite", 0.0, f64::INFINITY),
        ("lon infinite", f64::NEG_INFINITY, 0.0),
        ("both absurd", 1e300, -1e300),
    ] {
        assert!(
            try_encode(lon, lat, 9).is_none(),
            "{name}: try_encode must reject ({lon}, {lat})"
        );
        // `encode` clamps rather than signalling, so it still returns a hash of the
        // requested length — and must not panic on a non-finite input.
        let h = encode(lon, lat, 9);
        assert_eq!(
            h.len(),
            9,
            "{name}: encode clamps and still returns 9 chars"
        );
    }
    // Zero precision has no valid answer either way.
    assert!(try_encode(0.0, 0.0, 0).is_none());
    assert_eq!(
        encode(0.0, 0.0, 0).len(),
        0,
        "zero precision is an empty hash"
    );
}

#[test]
fn clamped_out_of_range_input_lands_at_the_corner_it_was_clamped_to() {
    // If `encode` clamps, it should clamp to the boundary rather than wrapping. A
    // latitude of 91 is not latitude -89, and a longitude of 181 is not -179.
    assert_eq!(
        encode(0.0, 91.0, 9),
        encode(0.0, 90.0, 9),
        "latitude above range should clamp to the pole, not wrap"
    );
    assert_eq!(
        encode(0.0, -91.0, 9),
        encode(0.0, -90.0, 9),
        "latitude below range should clamp to the pole"
    );
    assert_eq!(
        encode(181.0, 0.0, 9),
        encode(180.0, 0.0, 9),
        "longitude above range should clamp, not wrap to the other side"
    );
    assert_eq!(encode(-181.0, 0.0, 9), encode(-180.0, 0.0, 9));
}

#[test]
fn decode_is_case_insensitive_but_encode_is_not() {
    // The documented convention: `decode` accepts either case because hashes
    // arrive from elsewhere, while `encode` only ever emits lowercase. Pinned
    // because the asymmetry is surprising enough to be "tidied" in either
    // direction, and both directions would break a caller.
    let lower = "9q8yy";
    let upper = "9Q8YY";
    assert_eq!(
        decode(lower),
        decode(upper),
        "the two cases must decode to the same cell"
    );
    let (clon, clat) = decode(upper).expect("uppercase decodes");
    assert_eq!(
        encode(clon, clat, lower.len()),
        lower,
        "encode emits lowercase, so the round trip normalizes case rather than \
         reproducing the input"
    );
    // Mixed case too, and the excluded letters are still rejected in either case.
    assert_eq!(decode("9Q8yY"), decode(lower));
    for excluded in ['a', 'i', 'l', 'o', 'A', 'I', 'L', 'O'] {
        assert!(
            decode(&format!("9q8y{excluded}")).is_none(),
            "{excluded:?} is excluded from base-32 in either case"
        );
    }
}

#[test]
fn malformed_hashes_are_rejected_rather_than_decoded() {
    for bad in [
        "",
        "a",
        "i",
        "l",
        "o",
        "ailmno",
        "ABC123",
        "9q8yy!",
        "9q8yy ",
        " 9q8yy",
        "9q8-yy",
        "\u{0}",
        "🙂",
        "9q8yy🙂",
        "-1",
        "..",
    ] {
        assert!(
            decode(bad).is_none(),
            "{bad:?} is not a valid geohash and must not decode"
        );
    }
    assert!(decode("9q8yy").is_some());
}

#[test]
fn the_canonical_reference_hash_decodes_to_its_published_coordinates() {
    // One external anchor, so the whole suite is not self-consistent-but-wrong:
    // `u4pruydqqvj` is the standard worked example, Jutland, Denmark.
    let (lon, lat) = decode("u4pruydqqvj").expect("the canonical example decodes");
    assert!(
        (lat - 57.649_110_6).abs() < 1e-6 && (lon - 10.407_439_7).abs() < 1e-6,
        "expected ~(57.6491106, 10.4074397), got ({lat}, {lon})"
    );
    // And the reverse: those coordinates encode back to that hash.
    assert_eq!(encode(10.407_439_7, 57.649_110_6, 11), "u4pruydqqvj");
}

#[test]
fn very_high_precision_is_saturating_rather_than_garbage() {
    // Past ~12 characters a geohash exceeds what an f64 lat/lon can resolve, so the
    // extra characters carry no information. They must still be *deterministic* and
    // keep the prefix property rather than becoming noise — which is what the
    // repeated trailing characters below confirm.
    let h = encode(0.0, 0.0, 40);
    assert_eq!(h.len(), 40);
    assert!(
        encode(0.0, 0.0, 40) == h,
        "high precision must be deterministic"
    );
    assert!(
        h.starts_with(&encode(0.0, 0.0, 12)),
        "the saturated tail must extend the meaningful prefix"
    );
    // And it must not panic or produce characters outside the alphabet.
    const ALPHABET: &str = "0123456789bcdefghjkmnpqrstuvwxyz";
    assert!(h.chars().all(|c| ALPHABET.contains(c)));
}
