//! Polygon area on the round Earth, two ways — mirrors of `geo`'s traits.
//! Coordinates are interpreted as **lon/lat degrees** (`x = lon`, `y = lat`);
//! results are in **metres²**.
//!
//! - [`ChamberlainDuquetteArea`] — the cheap **spherical** approximation
//!   (JPL Publication 07-03), evaluated on a sphere of the equatorial Earth
//!   radius for exact agreement with geo. Off by up to ~0.5% vs the ellipsoid.
//! - [`GeodesicArea`] — the exact **ellipsoidal** (Karney) area and perimeter
//!   on WGS84, via [`hausdorff_base::earth::Ellipsoid::ring_perimeter_area_m`].

use hausdorff_base::earth::{Ellipsoid, GeoPoint};

use crate::space::r2::{Geometry, LineString, MultiPolygon, Polygon, Rect, Triangle};

/// Equatorial Earth radius (WGS84), the radius the Chamberlain–Duquette formula
/// is conventionally evaluated at — kept identical to geo for exact agreement.
pub const EQUATORIAL_EARTH_RADIUS_M: f64 = 6_378_137.0;

/// Spherical (Chamberlain–Duquette) polygon area, signed and unsigned.
pub trait ChamberlainDuquetteArea {
    /// Signed area in m²: sign follows the exterior ring's winding (holes
    /// subtract).
    fn chamberlain_duquette_signed_area(&self) -> f64;

    /// Absolute area in m².
    fn chamberlain_duquette_unsigned_area(&self) -> f64 {
        self.chamberlain_duquette_signed_area().abs()
    }
}

/// Signed spherical area contribution of a single (closed) ring.
fn ring_area(coords: &[LineStringCoord]) -> f64 {
    let n = coords.len();
    if n <= 2 {
        return 0.0;
    }
    let mut total = 0.0;
    for i in 0..n {
        let (lower, middle, upper) = if i == n - 2 {
            (n - 2, n - 1, 0)
        } else if i == n - 1 {
            (n - 1, 0, 1)
        } else {
            (i, i + 1, i + 2)
        };
        let p1 = coords[lower];
        let p2 = coords[middle];
        let p3 = coords[upper];
        total += (p3.0.to_radians() - p1.0.to_radians()) * p2.1.to_radians().sin();
    }
    total * EQUATORIAL_EARTH_RADIUS_M * EQUATORIAL_EARTH_RADIUS_M / -2.0
}

type LineStringCoord = (f64, f64);

fn ring_lonlat(ring: &LineString) -> Vec<LineStringCoord> {
    ring.0.iter().map(|c| (c.x, c.y)).collect()
}

impl ChamberlainDuquetteArea for Polygon {
    fn chamberlain_duquette_signed_area(&self) -> f64 {
        let exterior = ring_area(&ring_lonlat(self.exterior()));
        self.interiors().iter().fold(exterior, |total, hole| {
            total - ring_area(&ring_lonlat(hole))
        })
    }
}

impl ChamberlainDuquetteArea for MultiPolygon {
    fn chamberlain_duquette_signed_area(&self) -> f64 {
        self.iter()
            .map(ChamberlainDuquetteArea::chamberlain_duquette_signed_area)
            .sum()
    }
}

/// Exact ellipsoidal (Karney) area and perimeter on WGS84 — mirror of `geo`'s
/// `GeodesicArea`. Areas are in m², perimeters in metres.
///
/// The *signed* forms assume the simple-feature convention — exterior rings
/// counter-clockwise, holes clockwise — and return a negative area for a
/// clockwise exterior. The *unsigned* forms read every exterior as wound
/// counter-clockwise, which also makes them correct for polygons enclosing
/// more than half the planet (where the signed form wraps).
pub trait GeodesicArea {
    /// Signed geodesic area, m² (counter-clockwise exterior positive).
    fn geodesic_area_signed(&self) -> f64 {
        self.geodesic_perimeter_area_signed().1
    }

    /// Unsigned geodesic area, m².
    fn geodesic_area_unsigned(&self) -> f64 {
        self.geodesic_perimeter_area_unsigned().1
    }

    /// Geodesic perimeter, metres: every ring's boundary length (exterior
    /// plus holes), like [`Perimeter`](crate::space::r2::Perimeter).
    fn geodesic_perimeter(&self) -> f64 {
        self.geodesic_perimeter_area_signed().0
    }

    /// Perimeter and signed area in one pass.
    fn geodesic_perimeter_area_signed(&self) -> (f64, f64);

    /// Perimeter and unsigned area in one pass.
    fn geodesic_perimeter_area_unsigned(&self) -> (f64, f64);
}

/// Perimeter and area of `poly` under the geo convention: the exterior ring is
/// read as counter-clockwise-positive and the holes as clockwise-positive;
/// hole areas are subtracted by magnitude, following the exterior's sign.
fn polygon_perimeter_area(poly: &Polygon, signed: bool) -> (f64, f64) {
    let e = Ellipsoid::WGS84;
    let ring = |r: &LineString, ccw: bool| {
        // x = lon, y = lat → (lat, lon) for the geodesic kernel.
        let latlon: Vec<GeoPoint> =
            r.0.iter()
                .map(|c| GeoPoint::from_lon_lat(c.x, c.y))
                .collect();
        e.ring_perimeter_area_winding_m(&latlon, ccw, signed)
    };

    let (outer_perimeter, outer_area) = ring(poly.exterior(), true);
    let mut inner_perimeter = 0.0;
    let mut inner_area = 0.0;
    for hole in poly.interiors() {
        let (p, a) = ring(hole, false);
        inner_perimeter += p;
        inner_area += a.abs();
    }
    if outer_area < 0.0 && inner_area > 0.0 {
        inner_area = -inner_area;
    }
    (outer_perimeter + inner_perimeter, outer_area - inner_area)
}

impl GeodesicArea for Polygon {
    fn geodesic_perimeter_area_signed(&self) -> (f64, f64) {
        polygon_perimeter_area(self, true)
    }

    fn geodesic_perimeter_area_unsigned(&self) -> (f64, f64) {
        polygon_perimeter_area(self, false)
    }
}

impl GeodesicArea for MultiPolygon {
    fn geodesic_perimeter_area_signed(&self) -> (f64, f64) {
        self.iter()
            .map(GeodesicArea::geodesic_perimeter_area_signed)
            .fold((0.0, 0.0), |(p, a), (pi, ai)| (p + pi, a + ai))
    }

    fn geodesic_perimeter_area_unsigned(&self) -> (f64, f64) {
        self.iter()
            .map(GeodesicArea::geodesic_perimeter_area_unsigned)
            .fold((0.0, 0.0), |(p, a), (pi, ai)| (p + pi, a + ai))
    }
}

impl GeodesicArea for Rect {
    fn geodesic_perimeter_area_signed(&self) -> (f64, f64) {
        self.to_polygon().geodesic_perimeter_area_signed()
    }

    fn geodesic_perimeter_area_unsigned(&self) -> (f64, f64) {
        self.to_polygon().geodesic_perimeter_area_unsigned()
    }
}

impl GeodesicArea for Triangle {
    fn geodesic_perimeter_area_signed(&self) -> (f64, f64) {
        self.to_polygon().geodesic_perimeter_area_signed()
    }

    fn geodesic_perimeter_area_unsigned(&self) -> (f64, f64) {
        self.to_polygon().geodesic_perimeter_area_unsigned()
    }
}

/// Areal parts contribute; puntal and linear parts are zero; collections sum.
impl GeodesicArea for Geometry {
    fn geodesic_perimeter_area_signed(&self) -> (f64, f64) {
        geometry_perimeter_area(self, true)
    }

    fn geodesic_perimeter_area_unsigned(&self) -> (f64, f64) {
        geometry_perimeter_area(self, false)
    }
}

fn geometry_perimeter_area(g: &Geometry, signed: bool) -> (f64, f64) {
    let go = |x: &dyn GeodesicArea| {
        if signed {
            x.geodesic_perimeter_area_signed()
        } else {
            x.geodesic_perimeter_area_unsigned()
        }
    };
    match g {
        Geometry::Polygon(p) => go(p),
        Geometry::MultiPolygon(mp) => go(mp),
        Geometry::Rect(r) => go(r),
        Geometry::Triangle(t) => go(t),
        Geometry::GeometryCollection(gs) => gs
            .iter()
            .map(|sub| geometry_perimeter_area(sub, signed))
            .fold((0.0, 0.0), |(p, a), (pi, ai)| (p + pi, a + ai)),
        Geometry::Point(_)
        | Geometry::MultiPoint(_)
        | Geometry::Line(_)
        | Geometry::LineString(_)
        | Geometry::MultiLineString(_) => (0.0, 0.0),
    }
}
