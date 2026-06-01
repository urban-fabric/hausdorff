//! `geohash` — native geohash encoding/decoding of lon/lat points (Niemeyer's
//! public-domain base-32 algorithm), so the workspace need not depend on the
//! external `geohash` crate. A geohash is a string spatial index: each character
//! adds 5 bits, alternating longitude/latitude, narrowing a lon/lat cell.

const BASE32: &[u8; 32] = b"0123456789bcdefghjkmnpqrstuvwxyz";

/// Encode `(lon, lat)` in degrees to a geohash of `precision` characters.
#[must_use]
pub fn encode(lon: f64, lat: f64, precision: usize) -> String {
    let mut lat_range = (-90.0_f64, 90.0_f64);
    let mut lon_range = (-180.0_f64, 180.0_f64);
    let mut hash = String::with_capacity(precision);
    let mut even = true;
    let mut bit = 0u8;
    let mut ch = 0u8;

    while hash.len() < precision {
        if even {
            let mid = midpoint(lon_range);
            if lon >= mid {
                ch |= 1 << (4 - bit);
                lon_range.0 = mid;
            } else {
                lon_range.1 = mid;
            }
        } else {
            let mid = midpoint(lat_range);
            if lat >= mid {
                ch |= 1 << (4 - bit);
                lat_range.0 = mid;
            } else {
                lat_range.1 = mid;
            }
        }
        even = !even;
        if bit < 4 {
            bit += 1;
        } else {
            hash.push(BASE32[ch as usize] as char);
            bit = 0;
            ch = 0;
        }
    }
    hash
}

/// Decode a geohash to the `(lon, lat)` of its cell center, or `None` on an
/// invalid character.
#[must_use]
pub fn decode(hash: &str) -> Option<(f64, f64)> {
    let mut lat_range = (-90.0_f64, 90.0_f64);
    let mut lon_range = (-180.0_f64, 180.0_f64);
    let mut even = true;

    for c in hash.bytes() {
        let idx = BASE32.iter().position(|&b| b == c.to_ascii_lowercase())?;
        for bit in (0..5).rev() {
            let on = (idx >> bit) & 1 == 1;
            if even {
                let mid = midpoint(lon_range);
                if on {
                    lon_range.0 = mid;
                } else {
                    lon_range.1 = mid;
                }
            } else {
                let mid = midpoint(lat_range);
                if on {
                    lat_range.0 = mid;
                } else {
                    lat_range.1 = mid;
                }
            }
            even = !even;
        }
    }
    Some((midpoint(lon_range), midpoint(lat_range)))
}

#[inline]
fn midpoint(r: (f64, f64)) -> f64 {
    (r.0 + r.1) / 2.0
}
