//! Spherical cell measurement: the exact great-circle area of a cell boundary
//! ring (via Cagnoli's formula with Kahan summation), great-circle point
//! distance, and H3's published average-area / average-edge-length tables.

/// Authalic Earth radius (km) used by H3 for spherical area/length conversions.
pub(super) const EARTH_RADIUS_KM: f64 = 6371.007180918475;

pub(super) const AVG_HEXAGON_AREA_KM2: [f64; 16] = [
    4.357449416078383e+06,
    6.097884417941332e+05,
    8.680_178_039_899_72e4,
    1.239343465508816e+04,
    1.770347654491307e+03,
    2.529038581819449e+02,
    3.612906216441245e+01,
    5.161293359717191e+00,
    7.373275975944177e-01,
    1.053325134272067e-01,
    1.504750190766435e-02,
    2.149643129451879e-03,
    3.070_918_756_316_06e-4,
    4.387026794728296e-05,
    6.267181135324313e-06,
    8.953_115_907_605_79e-7,
];

pub(super) const AVG_HEXAGON_AREA_M2: [f64; 16] = [
    4.357_449_416_078_39e12,
    6.097884417941339e+11,
    8.680178039899731e+10,
    1.239343465508818e+10,
    1.770347654491309e+09,
    2.529038581819452e+08,
    3.612_906_216_441_25e7,
    5.161293359717198e+06,
    7.373275975944188e+05,
    1.053325134272069e+05,
    1.504750190766437e+04,
    2.149643129451882e+03,
    3.070918756316063e+02,
    4.387026794728301e+01,
    6.267181135324322e+00,
    8.953115907605802e-01,
];

pub(super) const AVG_HEXAGON_EDGE_LENGTH_KM: [f64; 16] = [
    1281.256011,
    483.0568391,
    182.5129565,
    68.97922179,
    26.07175968,
    9.854090990,
    3.724532667,
    1.406475763,
    0.531414010,
    0.200786148,
    0.075863783,
    0.028663897,
    0.010830188,
    0.004092010,
    0.001546100,
    0.000584169,
];

pub(super) const AVG_HEXAGON_EDGE_LENGTH_M: [f64; 16] = [
    1281256.011,
    483056.8391,
    182512.9565,
    68979.22179,
    26071.75968,
    9854.090990,
    3724.532667,
    1406.475763,
    531.4140101,
    200.7861476,
    75.86378287,
    28.66389748,
    10.83018784,
    4.092010473,
    1.546099657,
    0.584168630,
];

/// Great-circle distance in radians between two `(lat, lng)` points (radians),
/// via the haversine formula.
pub(super) fn distance_rads(a: (f64, f64), b: (f64, f64)) -> f64 {
    let sin_lat = ((b.0 - a.0) * 0.5).sin();
    let sin_lng = ((b.1 - a.1) * 0.5).sin();
    let value = sin_lat.mul_add(sin_lat, a.0.cos() * b.0.cos() * sin_lng * sin_lng);
    2.0 * value.sqrt().atan2((1.0 - value).sqrt())
}

/// Exact spherical area (rad²) enclosed by a `(lat, lng)` boundary ring,
/// summing Cagnoli's per-edge contribution with Kahan compensation.
pub(super) fn linear_ring_area_rads2(ring: &[(f64, f64)]) -> f64 {
    let mut adder = FloatAdder::default();
    for (a, b) in ring
        .iter()
        .zip(ring.iter().cycle().skip(1))
        .take(ring.len())
    {
        adder.add(cagnoli((a.1, a.0), (b.1, b.0)));
    }

    let area = adder.value();
    if area < 0.0 {
        area + 4.0 * std::f64::consts::PI
    } else {
        area
    }
}

fn cagnoli(a: (f64, f64), b: (f64, f64)) -> f64 {
    let a_lat = a.1.mul_add(0.5, std::f64::consts::PI * 0.25);
    let b_lat = b.1.mul_add(0.5, std::f64::consts::PI * 0.25);
    let sin_a = a_lat.sin() * b_lat.sin();
    let cos_a = a_lat.cos() * b_lat.cos();

    let delta = b.0 - a.0;
    -2.0 * (sin_a * delta.sin()).atan2(sin_a * delta.cos() + cos_a)
}

#[derive(Default)]
struct FloatAdder {
    sum: f64,
    correction: f64,
}

impl FloatAdder {
    fn add(&mut self, value: f64) {
        let y = value - self.correction;
        let t = self.sum + y;
        self.correction = (t - self.sum) - y;
        self.sum = t;
    }

    fn value(self) -> f64 {
        self.sum + self.correction
    }
}
