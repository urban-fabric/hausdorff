//! Metric spaces and the measures over them — [`Euclidean`] (planar) and
//! [`Haversine`] (great-circle) — mirroring `geo`'s `line_measures`:
//! [`Distance`], [`Length`], [`Bearing`], [`Destination`].
//!
//! Under [`Haversine`], a [`Point`]/[`Coord`] is read as **`x = lon`, `y = lat`
//! degrees** (the geo convention) and the spherical work is delegated to the
//! crate's validated [`crate::earth`] kernels, so results stay consistent with the
//! rest of the crate.
//!
//! The ellipsoidal `Geodesic` (Karney) metric space is a future addition — it
//! needs a geographiclib-style port rather than the spherical kernels here.

use crate::earth::{Bearing as CompassBearing, GeoPoint, destination_m, haversine_m};

use super::{Coord, Geometry, LineString, MultiLineString, MultiPolygon, Point, Polygon};

/// The planar Euclidean metric space.
#[derive(Clone, Copy, Debug, Default)]
pub struct Euclidean;

/// The spherical great-circle metric space (coordinates read as `x = lon`,
/// `y = lat`), on a sphere of [`crate::earth::EARTH_RADIUS_M`].
#[derive(Clone, Copy, Debug, Default)]
pub struct Haversine;

/// A `(x = lon, y = lat)` planar point as a spherical [`GeoPoint`].
#[inline]
fn as_geo(c: Coord) -> GeoPoint {
    GeoPoint::from_lat_lon(c.y, c.x)
}

/// Distance between two points under a metric space (the metres unit for
/// [`Haversine`], the plane's own unit for [`Euclidean`]).
pub trait Distance<A, B> {
    fn distance(&self, a: A, b: B) -> f64;
}

impl Distance<Coord, Coord> for Euclidean {
    fn distance(&self, a: Coord, b: Coord) -> f64 {
        let d = a - b;
        (d.x * d.x + d.y * d.y).sqrt()
    }
}

impl Distance<Point, Point> for Euclidean {
    fn distance(&self, a: Point, b: Point) -> f64 {
        self.distance(a.0, b.0)
    }
}

impl Distance<Coord, Coord> for Haversine {
    fn distance(&self, a: Coord, b: Coord) -> f64 {
        haversine_m(as_geo(a), as_geo(b))
    }
}

impl Distance<Point, Point> for Haversine {
    fn distance(&self, a: Point, b: Point) -> f64 {
        haversine_m(as_geo(a.0), as_geo(b.0))
    }
}

impl Distance<GeoPoint, GeoPoint> for Haversine {
    fn distance(&self, a: GeoPoint, b: GeoPoint) -> f64 {
        haversine_m(a, b)
    }
}

/// Total length of a linear geometry under a metric space.
pub trait Length<G> {
    fn length(&self, geometry: &G) -> f64;
}

impl Length<LineString> for Euclidean {
    fn length(&self, ls: &LineString) -> f64 {
        ls.0.windows(2).map(|w| self.distance(w[0], w[1])).sum()
    }
}

impl Length<LineString> for Haversine {
    fn length(&self, ls: &LineString) -> f64 {
        ls.0.windows(2).map(|w| self.distance(w[0], w[1])).sum()
    }
}

impl<M> Length<MultiLineString> for M
where
    M: Length<LineString>,
{
    fn length(&self, mls: &MultiLineString) -> f64 {
        mls.iter().map(|ls| self.length(ls)).sum()
    }
}

/// Total length of the **curve** parts of a geometry under metric `M` — the 1-D
/// measure of its linear content (PostGIS `ST_Length`, matching `geo`'s `Length`).
/// Areal and puntal geometries contribute `0`; a collection sums its linear
/// members. A polygon's ring boundary is *not* counted here — that is
/// [`Perimeter`]. (`length + perimeter` is the JTS `getLength` total.)
impl<M> Length<Geometry> for M
where
    M: Length<LineString>,
{
    fn length(&self, g: &Geometry) -> f64 {
        geometry_length(self, g)
    }
}

/// Linear length of `g` under metric `m` (curve parts only). A free function
/// (recursing on itself for a collection) so every inner call is the unambiguous
/// `Length<LineString>`, not the `Length<Geometry>` impl being defined.
fn geometry_length<M: Length<LineString>>(m: &M, g: &Geometry) -> f64 {
    match g {
        Geometry::Line(l) => m.length(&LineString::new(vec![l.start, l.end])),
        Geometry::LineString(ls) => m.length(ls),
        Geometry::MultiLineString(mls) => mls.iter().map(|ls| m.length(ls)).sum(),
        Geometry::GeometryCollection(gs) => gs.iter().map(|sub| geometry_length(m, sub)).sum(),
        // points and areal geometries have no curve length
        Geometry::Point(_)
        | Geometry::MultiPoint(_)
        | Geometry::Polygon(_)
        | Geometry::MultiPolygon(_)
        | Geometry::Rect(_)
        | Geometry::Triangle(_) => 0.0,
    }
}

/// Total boundary length of the **areal** parts of a geometry under a metric —
/// a polygon's perimeter (exterior + holes), summed over the parts (PostGIS
/// `ST_Perimeter`). Linear and puntal geometries contribute `0`. The areal
/// counterpart of [`Length`]; their sum is the JTS `getLength` total.
pub trait Perimeter<G> {
    fn perimeter(&self, geometry: &G) -> f64;
}

impl<M> Perimeter<Polygon> for M
where
    M: Length<LineString>,
{
    fn perimeter(&self, p: &Polygon) -> f64 {
        polygon_perimeter(self, p)
    }
}

impl<M> Perimeter<MultiPolygon> for M
where
    M: Length<LineString>,
{
    fn perimeter(&self, mp: &MultiPolygon) -> f64 {
        mp.iter().map(|p| polygon_perimeter(self, p)).sum()
    }
}

impl<M> Perimeter<Geometry> for M
where
    M: Length<LineString>,
{
    fn perimeter(&self, g: &Geometry) -> f64 {
        geometry_perimeter(self, g)
    }
}

/// Areal boundary length of `g` under metric `m` (perimeter of its polygonal
/// parts). A free function (recursing on itself for a collection) so every inner
/// call resolves unambiguously.
fn geometry_perimeter<M: Length<LineString>>(m: &M, g: &Geometry) -> f64 {
    match g {
        Geometry::Polygon(p) => polygon_perimeter(m, p),
        Geometry::MultiPolygon(mp) => mp.iter().map(|p| polygon_perimeter(m, p)).sum(),
        Geometry::Rect(r) => polygon_perimeter(m, &r.to_polygon()),
        Geometry::Triangle(t) => polygon_perimeter(m, &t.to_polygon()),
        Geometry::GeometryCollection(gs) => gs.iter().map(|sub| geometry_perimeter(m, sub)).sum(),
        // points and curves have no areal boundary
        Geometry::Point(_)
        | Geometry::MultiPoint(_)
        | Geometry::Line(_)
        | Geometry::LineString(_)
        | Geometry::MultiLineString(_) => 0.0,
    }
}

/// A polygon's perimeter under metric `M`: its exterior ring plus every hole.
fn polygon_perimeter<M: Length<LineString>>(m: &M, p: &Polygon) -> f64 {
    m.length(p.exterior()) + p.interiors().iter().map(|h| m.length(h)).sum::<f64>()
}

/// Initial bearing in **degrees** (`0 = North`, clockwise) from `a` to `b`.
pub trait Bearing<A, B> {
    fn bearing(&self, a: A, b: B) -> f64;
}

impl Bearing<Point, Point> for Haversine {
    fn bearing(&self, a: Point, b: Point) -> f64 {
        CompassBearing::initial(as_geo(a.0), as_geo(b.0)).degrees()
    }
}

impl Bearing<Coord, Coord> for Haversine {
    fn bearing(&self, a: Coord, b: Coord) -> f64 {
        CompassBearing::initial(as_geo(a), as_geo(b)).degrees()
    }
}

/// The forward geodesic: the point reached from `origin` by travelling
/// `distance` along the initial compass `bearing_degrees`.
pub trait Destination {
    fn destination(&self, origin: Point, bearing_degrees: f64, distance: f64) -> Point;
}

impl Destination for Haversine {
    fn destination(&self, origin: Point, bearing_degrees: f64, distance: f64) -> Point {
        let dest = destination_m(
            as_geo(origin.0),
            CompassBearing::from_degrees(bearing_degrees),
            distance,
        );
        Point::new(dest.lon, dest.lat) // back to x = lon, y = lat
    }
}
