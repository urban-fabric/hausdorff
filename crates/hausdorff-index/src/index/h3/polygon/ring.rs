use super::H3PolygonError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct BBox {
    min_lat: f64,
    max_lat: f64,
    min_lng: f64,
    max_lng: f64,
}

impl BBox {
    pub(super) fn from_rings<'a>(
        rings: impl IntoIterator<Item = &'a Vec<(f64, f64)>>,
    ) -> Result<Self, H3PolygonError> {
        let mut points = Vec::new();
        for ring in rings {
            points.extend(ring.iter().copied());
        }
        Self::from_points(points).ok_or(H3PolygonError::RingTooShort)
    }

    pub(super) fn from_points(points: impl IntoIterator<Item = (f64, f64)>) -> Option<Self> {
        let mut bbox = None::<Self>;
        for (lat, lng) in points {
            bbox = Some(match bbox {
                Some(current) => Self {
                    min_lat: current.min_lat.min(lat),
                    max_lat: current.max_lat.max(lat),
                    min_lng: current.min_lng.min(lng),
                    max_lng: current.max_lng.max(lng),
                },
                None => Self {
                    min_lat: lat,
                    max_lat: lat,
                    min_lng: lng,
                    max_lng: lng,
                },
            });
        }
        bbox
    }

    pub(super) fn contains(self, point: (f64, f64)) -> bool {
        point.0 >= self.min_lat
            && point.0 <= self.max_lat
            && point.1 >= self.min_lng
            && point.1 <= self.max_lng
    }

    pub(super) fn intersects(self, other: &Self) -> bool {
        self.min_lat <= other.max_lat
            && self.max_lat >= other.min_lat
            && self.min_lng <= other.max_lng
            && self.max_lng >= other.min_lng
    }

    pub(super) fn center(self) -> (f64, f64) {
        (
            (self.min_lat + self.max_lat) * 0.5,
            (self.min_lng + self.max_lng) * 0.5,
        )
    }
}

pub(super) fn normalize_ring(mut ring: Vec<(f64, f64)>) -> Result<Vec<(f64, f64)>, H3PolygonError> {
    if ring.first() == ring.last() {
        ring.pop();
    }
    if ring.len() < 3 {
        return Err(H3PolygonError::RingTooShort);
    }
    for &(lat, lng) in &ring {
        if !lat.is_finite()
            || !lng.is_finite()
            || !(-90.0..=90.0).contains(&lat)
            || !(-180.0..=180.0).contains(&lng)
        {
            return Err(H3PolygonError::InvalidCoordinate);
        }
    }
    Ok(ring)
}

/// True iff some ring edge spans more than 180° of longitude, meaning the ring
/// takes the short way across the antimeridian (mirrors the H3/h3o heuristic).
pub(super) fn ring_is_transmeridian(ring: &[(f64, f64)]) -> bool {
    ring_edges(ring).any(|(a, b)| (a.1 - b.1).abs() > 180.0)
}

/// Shifts west-hemisphere longitudes by `+360` so an antimeridian-crossing
/// ring becomes continuous.
pub(super) fn shift_ring_east(ring: &mut [(f64, f64)]) {
    for point in ring {
        if point.1 < 0.0 {
            point.1 += 360.0;
        }
    }
}

/// Ring edges as consecutive point pairs, closing back to the first point.
pub(super) fn ring_edges(ring: &[(f64, f64)]) -> impl Iterator<Item = ((f64, f64), (f64, f64))> {
    ring.iter()
        .copied()
        .zip(ring.iter().copied().cycle().skip(1))
        .take(ring.len())
}

/// Twice the signed area of triangle `abc` in the lat/lng plane.
fn orientation(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0).mul_add(c.1 - a.1, -((b.1 - a.1) * (c.0 - a.0)))
}

/// True iff segments `ab` and `cd` cross at an interior point of both —
/// endpoint touches and collinear overlaps do not count, so a shared boundary
/// does not break full containment.
pub(super) fn segments_cross_properly(
    a: (f64, f64),
    b: (f64, f64),
    c: (f64, f64),
    d: (f64, f64),
) -> bool {
    let d1 = orientation(c, d, a);
    let d2 = orientation(c, d, b);
    let d3 = orientation(a, b, c);
    let d4 = orientation(a, b, d);
    ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
}

/// True iff segments `ab` and `cd` intersect at all, including endpoint
/// touches and collinear overlaps.
pub(super) fn segments_intersect(
    a: (f64, f64),
    b: (f64, f64),
    c: (f64, f64),
    d: (f64, f64),
) -> bool {
    if segments_cross_properly(a, b, c, d) {
        return true;
    }
    let d1 = orientation(c, d, a);
    let d2 = orientation(c, d, b);
    let d3 = orientation(a, b, c);
    let d4 = orientation(a, b, d);
    (d1 == 0.0 && collinear_point_on_segment(c, d, a))
        || (d2 == 0.0 && collinear_point_on_segment(c, d, b))
        || (d3 == 0.0 && collinear_point_on_segment(a, b, c))
        || (d4 == 0.0 && collinear_point_on_segment(a, b, d))
}

/// Bounds check for a point already known to be collinear with segment `ab`.
fn collinear_point_on_segment(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> bool {
    p.0 >= a.0.min(b.0) && p.0 <= a.0.max(b.0) && p.1 >= a.1.min(b.1) && p.1 <= a.1.max(b.1)
}

pub(super) fn ring_contains_point(ring: &[(f64, f64)], point: (f64, f64)) -> bool {
    let (lat, lng) = point;
    let mut inside = false;

    for (&a, &b) in ring
        .iter()
        .zip(ring.iter().cycle().skip(1))
        .take(ring.len())
    {
        if point_on_segment(a, b, point) {
            return true;
        }

        let (lat_a, lng_a) = a;
        let (lat_b, lng_b) = b;
        let crosses = (lat_a > lat) != (lat_b > lat);
        if crosses {
            let lng_at_lat = (lng_b - lng_a).mul_add((lat - lat_a) / (lat_b - lat_a), lng_a);
            if lng < lng_at_lat {
                inside = !inside;
            }
        }
    }

    inside
}

fn point_on_segment(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> bool {
    const EPS: f64 = 1e-12;
    let cross = (p.0 - a.0).mul_add(b.1 - a.1, -((p.1 - a.1) * (b.0 - a.0)));
    if cross.abs() > EPS {
        return false;
    }
    let min_lat = a.0.min(b.0) - EPS;
    let max_lat = a.0.max(b.0) + EPS;
    let min_lng = a.1.min(b.1) - EPS;
    let max_lng = a.1.max(b.1) + EPS;
    p.0 >= min_lat && p.0 <= max_lat && p.1 >= min_lng && p.1 <= max_lng
}

pub(super) fn linear_ring_area_rads2(ring: &[(f64, f64)]) -> f64 {
    let area = euclid_numeric::compensated::sum_iter(
        ring.iter()
            .zip(ring.iter().cycle().skip(1))
            .take(ring.len())
            .map(|(a, b)| {
                cagnoli(
                    (a.1.to_radians(), a.0.to_radians()),
                    (b.1.to_radians(), b.0.to_radians()),
                )
            }),
    );
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
