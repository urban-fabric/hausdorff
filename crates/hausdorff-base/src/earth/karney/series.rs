//! The Fourier-series machinery: the expansion order, the ellipsoid-independent
//! `A₁`/`A₂`/`C₁`/`C₁'`/`C₂` coefficient generators, and the raw `A₃`/`C₃`/`C₄`
//! polynomial coefficient tables (in the third flattening) the solver pre-derives
//! per ellipsoid.

use super::math::{polyval, sq};

/// Order of the series expansions (GeographicLib's compile-time default).
pub(super) const GEODESIC_ORDER: usize = 6;
/// The `C1`/`C1'`/`C2` Fourier coefficient arrays are 1-indexed: slot 0 unused.
pub(super) const CARR_SIZE: usize = GEODESIC_ORDER + 1;
/// Number of entries in the `C3` coefficient triangle.
pub(super) const NC3X: usize = 15;
/// Number of entries in the `C4` (area series) coefficient triangle.
pub(super) const NC4X: usize = 21;

/// `A₁(ε) − 1`: the mean-to-true distance scale of the series.
pub(super) fn a1m1f(eps: f64) -> f64 {
    const COEFF: [f64; 5] = [1.0, 4.0, 64.0, 0.0, 256.0];
    let m = GEODESIC_ORDER / 2;
    let t = polyval(m, &COEFF, sq(eps)) / COEFF[m + 1];
    (t + eps) / (1.0 - eps)
}

/// The `C₁` Fourier coefficients (distance series).
pub(super) fn c1f(eps: f64, c: &mut [f64; CARR_SIZE]) {
    const COEFF: [f64; 18] = [
        -1.0, 6.0, -16.0, 32.0, -9.0, 64.0, -128.0, 2048.0, 9.0, -16.0, 768.0, 3.0, -5.0, 512.0,
        -7.0, 1280.0, -7.0, 2048.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    // Indexed loop as in the reference port (the iterator form benched slower
    // there).
    #[allow(clippy::needless_range_loop)]
    for l in 1..=GEODESIC_ORDER {
        let m = (GEODESIC_ORDER - l) / 2;
        c[l] = d * polyval(m, &COEFF[o..], eps2) / COEFF[o + m + 1];
        o += m + 2;
        d *= eps;
    }
}

/// The `C₁'` Fourier coefficients (the inverse of the distance series, for the
/// direct problem).
pub(super) fn c1pf(eps: f64, c: &mut [f64; CARR_SIZE]) {
    const COEFF: [f64; 18] = [
        205.0, -432.0, 768.0, 1536.0, 4005.0, -4736.0, 3840.0, 12288.0, -225.0, 116.0, 384.0,
        -7173.0, 2695.0, 7680.0, 3467.0, 7680.0, 38081.0, 61440.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    // Indexed loop as in the reference port (the iterator form benched slower
    // there).
    #[allow(clippy::needless_range_loop)]
    for l in 1..=GEODESIC_ORDER {
        let m = (GEODESIC_ORDER - l) / 2;
        c[l] = d * polyval(m, &COEFF[o..], eps2) / COEFF[o + m + 1];
        o += m + 2;
        d *= eps;
    }
}

/// `A₂(ε) − 1`: the scale of the reduced-length series.
pub(super) fn a2m1f(eps: f64) -> f64 {
    const COEFF: [f64; 5] = [-11.0, -28.0, -192.0, 0.0, 256.0];
    let m = GEODESIC_ORDER / 2;
    let t = polyval(m, &COEFF, sq(eps)) / COEFF[m + 1];
    (t - eps) / (1.0 + eps)
}

/// The `C₂` Fourier coefficients (reduced-length series).
pub(super) fn c2f(eps: f64, c: &mut [f64; CARR_SIZE]) {
    const COEFF: [f64; 18] = [
        1.0, 2.0, 16.0, 32.0, 35.0, 64.0, 384.0, 2048.0, 15.0, 80.0, 768.0, 7.0, 35.0, 512.0, 63.0,
        1280.0, 77.0, 2048.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    // Indexed loop as in the reference port (the iterator form benched slower
    // there).
    #[allow(clippy::needless_range_loop)]
    for l in 1..=GEODESIC_ORDER {
        let m = (GEODESIC_ORDER - l) / 2;
        c[l] = d * polyval(m, &COEFF[o..], eps2) / COEFF[o + m + 1];
        o += m + 2;
        d *= eps;
    }
}

/// Coefficients of the `A₃` polynomial in the third flattening `n`.
pub(super) const COEFF_A3: [f64; 18] = [
    -3.0, 128.0, -2.0, -3.0, 64.0, -1.0, -3.0, -1.0, 16.0, 3.0, -1.0, -2.0, 8.0, 1.0, -1.0, 2.0,
    1.0, 1.0,
];

/// Coefficients of the `C₃` polynomials in the third flattening `n`.
pub(super) const COEFF_C3: [f64; 45] = [
    3.0, 128.0, 2.0, 5.0, 128.0, -1.0, 3.0, 3.0, 64.0, -1.0, 0.0, 1.0, 8.0, -1.0, 1.0, 4.0, 5.0,
    256.0, 1.0, 3.0, 128.0, -3.0, -2.0, 3.0, 64.0, 1.0, -3.0, 2.0, 32.0, 7.0, 512.0, -10.0, 9.0,
    384.0, 5.0, -9.0, 5.0, 192.0, 7.0, 512.0, -14.0, 7.0, 512.0, 21.0, 2560.0,
];

/// Coefficients of the `C₄` polynomials (area series) in the third flattening.
pub(super) const COEFF_C4: [f64; 77] = [
    97.0, 15015.0, 1088.0, 156.0, 45045.0, -224.0, -4784.0, 1573.0, 45045.0, -10656.0, 14144.0,
    -4576.0, -858.0, 45045.0, 64.0, 624.0, -4576.0, 6864.0, -3003.0, 15015.0, 100.0, 208.0, 572.0,
    3432.0, -12012.0, 30030.0, 45045.0, 1.0, 9009.0, -2944.0, 468.0, 135135.0, 5792.0, 1040.0,
    -1287.0, 135135.0, 5952.0, -11648.0, 9152.0, -2574.0, 135135.0, -64.0, -624.0, 4576.0, -6864.0,
    3003.0, 135135.0, 8.0, 10725.0, 1856.0, -936.0, 225225.0, -8448.0, 4992.0, -1144.0, 225225.0,
    -1440.0, 4160.0, -4576.0, 1716.0, 225225.0, -136.0, 63063.0, 1024.0, -208.0, 105105.0, 3584.0,
    -3328.0, 1144.0, 315315.0, -128.0, 135135.0, -2560.0, 832.0, 405405.0, 128.0, 99099.0,
];
