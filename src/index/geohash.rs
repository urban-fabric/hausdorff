//! `geohash` — native geohash encoding/decoding of lon/lat points (Niemeyer's
//! public-domain base-32 algorithm), so the workspace need not depend on the
//! external `geohash` crate. A geohash is a string spatial index: each character
//! adds 5 bits, alternating longitude/latitude, narrowing a lon/lat cell.

const BASE32: &[u8; 32] = b"0123456789bcdefghjkmnpqrstuvwxyz";

/// Encode `(lon, lat)` in degrees to a geohash of `precision` characters,
/// returning `None` for a zero `precision`, non-finite coordinates, or
/// coordinates outside `lon ∈ [-180, 180]` / `lat ∈ [-90, 90]`. The checked
/// form of [`encode`], which clamps such inputs silently.
#[must_use]
pub fn try_encode(lon: f64, lat: f64, precision: usize) -> Option<String> {
    if precision == 0
        || !lon.is_finite()
        || !lat.is_finite()
        || !(-180.0..=180.0).contains(&lon)
        || !(-90.0..=90.0).contains(&lat)
    {
        return None;
    }
    Some(encode(lon, lat, precision))
}

/// Largest precision the fixed-point (Morton) fast path can serve: 12 chars =
/// 60 bits ≤ the 64-bit interleaved code.
const MORTON_MAX_PRECISION: usize = 12;

/// Encode `(lon, lat)` in degrees to a geohash of `precision` characters.
///
/// Up to 12 chars this bit-interleaves 32-bit fixed-point coordinates
/// (Morton/Z-order) — a few bit ops instead of `precision × 5` bisection steps.
/// Longer hashes fall back to a bisection reference path.
#[must_use]
pub fn encode(lon: f64, lat: f64, precision: usize) -> String {
    if precision == 0 {
        return String::new();
    }
    if precision <= MORTON_MAX_PRECISION {
        return encode_morton(lon, lat, precision);
    }
    encode_bisect(lon, lat, precision)
}

/// Spread the 32 bits of `v` into the even bit positions of a `u64`.
#[inline]
fn spread(v: u32) -> u64 {
    let mut x = u64::from(v);
    x = (x | (x << 16)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x << 8)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x << 2)) & 0x3333_3333_3333_3333;
    x = (x | (x << 1)) & 0x5555_5555_5555_5555;
    x
}

/// Map a fraction in `[0, 1]` to 32-bit fixed point — the geohash bisection bits,
/// MSB first.
#[inline]
fn quantize(frac: f64) -> u32 {
    (frac * 4_294_967_296.0).clamp(0.0, 4_294_967_295.0) as u32
}

fn encode_morton(lon: f64, lat: f64, precision: usize) -> String {
    let lon_q = quantize((lon + 180.0) / 360.0);
    let lat_q = quantize((lat + 90.0) / 180.0);
    let morton = (spread(lon_q) << 1) | spread(lat_q);
    let mut hash = String::with_capacity(precision);
    for j in 0..precision {
        let symbol = ((morton >> (64 - 5 * (j + 1))) & 0x1F) as usize;
        hash.push(BASE32[symbol] as char);
    }
    hash
}

/// The bisection reference encoder: narrows the lon/lat cell five bits per char.
fn encode_bisect(lon: f64, lat: f64, precision: usize) -> String {
    let mut lat_range = (-90.0_f64, 90.0_f64);
    let mut lon_range = (-180.0_f64, 180.0_f64);
    let mut hash = String::with_capacity(precision);
    let mut even = true;
    for _ in 0..precision {
        let mut ch = 0u8;
        for bit in 0..5u8 {
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
        }
        hash.push(BASE32[ch as usize] as char);
    }
    hash
}

/// Decode a geohash to the `(lon, lat)` of its cell center, or `None` for an
/// empty string or an invalid character.
#[must_use]
pub fn decode(hash: &str) -> Option<(f64, f64)> {
    if hash.is_empty() {
        return None;
    }
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
