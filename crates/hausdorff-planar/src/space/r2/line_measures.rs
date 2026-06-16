//! Metric spaces and the measures over them — [`Euclidean`] (planar),
//! [`Haversine`] (great-circle), [`Rhumb`] (loxodrome), and [`Geodesic`]
//! (WGS84 ellipsoidal, Karney) — mirroring `geo`'s `line_measures`:
//! [`Distance`], [`Length`], [`Bearing`], [`Destination`],
//! [`InterpolatePoint`].
//!
//! Under the geographic metrics, a [`Point`]/[`Coord`] is read as **`x = lon`,
//! `y = lat` degrees** (the geo convention) and the spherical/ellipsoidal work
//! is delegated to the crate's validated [`hausdorff_base::earth`] kernels, so results
//! stay consistent with the rest of the crate.

use hausdorff_base::earth::{
    Bearing as CompassBearing, Ellipsoid, GeoPoint, destination_m, haversine_m, rhumb_bearing,
    rhumb_destination_m, rhumb_distance_m, rhumb_point_at_ratio,
};

use super::{Coord, Geometry, LineString, MultiLineString, MultiPolygon, Point, Polygon};

/// The planar Euclidean metric space.
#[derive(Clone, Copy, Debug, Default)]
pub struct Euclidean;

/// The spherical great-circle metric space (coordinates read as `x = lon`,
/// `y = lat`), on a sphere of [`hausdorff_base::earth::EARTH_RADIUS_M`].
#[derive(Clone, Copy, Debug, Default)]
pub struct Haversine;

/// The spherical rhumb-line (loxodrome) metric space (coordinates read as
/// `x = lon`, `y = lat`), on a sphere of [`hausdorff_base::earth::EARTH_RADIUS_M`]:
/// constant-bearing courses, straight lines on a Mercator chart.
#[derive(Clone, Copy, Debug, Default)]
pub struct Rhumb;
/// The WGS84 **ellipsoidal** geodesic metric space (coordinates read as
/// `x = lon`, `y = lat`) — Karney's algorithms via
/// [`Ellipsoid::WGS84`](hausdorff_base::earth::Ellipsoid::WGS84). Exact to round-off
/// where [`Haversine`]'s spherical model is off by up to ~0.5%; roughly an
/// order of magnitude slower.
#[derive(Clone, Copy, Debug, Default)]
pub struct Geodesic;

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

impl Distance<Coord, Coord> for Rhumb {
    fn distance(&self, a: Coord, b: Coord) -> f64 {
        rhumb_distance_m(as_geo(a), as_geo(b))
    }
}

impl Distance<Point, Point> for Rhumb {
    fn distance(&self, a: Point, b: Point) -> f64 {
        rhumb_distance_m(as_geo(a.0), as_geo(b.0))
    }
}

impl Distance<GeoPoint, GeoPoint> for Rhumb {
    fn distance(&self, a: GeoPoint, b: GeoPoint) -> f64 {
        rhumb_distance_m(a, b)
    }
}

impl Distance<Coord, Coord> for Geodesic {
    fn distance(&self, a: Coord, b: Coord) -> f64 {
        Ellipsoid::WGS84.distance_m(as_geo(a), as_geo(b))
    }
}

impl Distance<Point, Point> for Geodesic {
    fn distance(&self, a: Point, b: Point) -> f64 {
        self.distance(a.0, b.0)
    }
}

impl Distance<GeoPoint, GeoPoint> for Geodesic {
    fn distance(&self, a: GeoPoint, b: GeoPoint) -> f64 {
        Ellipsoid::WGS84.distance_m(a, b)
    }
}

/// Total length of a linear geometry under a metric space.
///
/// Proof-backed (Euclidean) by `riemann` `Hausdorff/Length`: the polyline length
/// `Σ‖vᵢ₊₁−vᵢ‖` is non-negative (`polyLen_nonneg`), a two-vertex line measures the
/// endpoint distance (`polyLen_pair`), and obeys the **chord bound** — the
/// straight-line endpoint distance is at most the polyline length
/// (`dist_le_polyLen`), so a refined or coarsened path can't drift past its
/// length. (The sphere instance is `Hausdorff/SpherePolyline`.)
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

impl Length<LineString> for Rhumb {
    fn length(&self, ls: &LineString) -> f64 {
        ls.0.windows(2).map(|w| self.distance(w[0], w[1])).sum()
    }
}

impl Length<LineString> for Geodesic {
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

impl Bearing<Point, Point> for Rhumb {
    fn bearing(&self, a: Point, b: Point) -> f64 {
        rhumb_bearing(as_geo(a.0), as_geo(b.0)).degrees()
    }
}

impl Bearing<Coord, Coord> for Rhumb {
    fn bearing(&self, a: Coord, b: Coord) -> f64 {
        rhumb_bearing(as_geo(a), as_geo(b)).degrees()
    }
}

impl Bearing<Point, Point> for Geodesic {
    fn bearing(&self, a: Point, b: Point) -> f64 {
        self.bearing(a.0, b.0)
    }
}

impl Bearing<Coord, Coord> for Geodesic {
    fn bearing(&self, a: Coord, b: Coord) -> f64 {
        Ellipsoid::WGS84
            .inverse(as_geo(a), as_geo(b))
            .initial_bearing
            .degrees()
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
        Point::new(dest.lon, dest.lat)
    }
}

impl Destination for Rhumb {
    fn destination(&self, origin: Point, bearing_degrees: f64, distance: f64) -> Point {
        let dest = rhumb_destination_m(
            as_geo(origin.0),
            CompassBearing::from_degrees(bearing_degrees),
            distance,
        );
        Point::new(dest.lon, dest.lat)
    }
}

impl Destination for Geodesic {
    fn destination(&self, origin: Point, bearing_degrees: f64, distance: f64) -> Point {
        let dest = Ellipsoid::WGS84.destination_m(
            as_geo(origin.0),
            CompassBearing::from_degrees(bearing_degrees),
            distance,
        );
        Point::new(dest.lon, dest.lat)
    }
}

/// Point interpolation along a metric space's "straight line" between two
/// endpoints — the basis for [`Densify`](super::densify::Densify).
pub trait InterpolatePoint {
    /// The point `ratio_from_start` of the way (by this metric's distance)
    /// from `start` to `end`: `0` is `start`, `1` is `end`.
    fn point_at_ratio_between(&self, start: Point, end: Point, ratio_from_start: f64) -> Point;

    /// The point `distance_from_start` along the line from `start` toward
    /// `end` (the metres unit for spherical metrics, the plane's own unit for
    /// [`Euclidean`]).
    fn point_at_distance_between(
        &self,
        start: Point,
        end: Point,
        distance_from_start: f64,
    ) -> Point;
}

impl InterpolatePoint for Euclidean {
    fn point_at_ratio_between(&self, start: Point, end: Point, ratio_from_start: f64) -> Point {
        let d = end.0 - start.0;
        Point::new(
            start.x() + d.x * ratio_from_start,
            start.y() + d.y * ratio_from_start,
        )
    }

    fn point_at_distance_between(
        &self,
        start: Point,
        end: Point,
        distance_from_start: f64,
    ) -> Point {
        let total = self.distance(start, end);
        self.point_at_ratio_between(start, end, distance_from_start / total)
    }
}

/// Great-circle interpolation via the spherical destination kernel.
///
/// Proof-backed by `riemann` `Hausdorff/Slerp.slerp_unit` / `slerp_angle` and
/// `Hausdorff/Destination.dest_unit` / `dest_dot`: the waypoint stays on the
/// sphere and is reached at the requested central angle from the origin.
impl InterpolatePoint for Haversine {
    fn point_at_ratio_between(&self, start: Point, end: Point, ratio_from_start: f64) -> Point {
        let start_geo = as_geo(start.0);
        let end_geo = as_geo(end.0);
        let total = haversine_m(start_geo, end_geo);
        let bearing = CompassBearing::initial(start_geo, end_geo);
        let p = destination_m(start_geo, bearing, total * ratio_from_start);
        Point::new(p.lon, p.lat)
    }

    fn point_at_distance_between(
        &self,
        start: Point,
        end: Point,
        distance_from_start: f64,
    ) -> Point {
        let start_geo = as_geo(start.0);
        let end_geo = as_geo(end.0);
        let bearing = CompassBearing::initial(start_geo, end_geo);
        let p = destination_m(start_geo, bearing, distance_from_start);
        Point::new(p.lon, p.lat)
    }
}

impl InterpolatePoint for Rhumb {
    fn point_at_ratio_between(&self, start: Point, end: Point, ratio_from_start: f64) -> Point {
        let p = rhumb_point_at_ratio(as_geo(start.0), as_geo(end.0), ratio_from_start);
        Point::new(p.lon, p.lat)
    }

    fn point_at_distance_between(
        &self,
        start: Point,
        end: Point,
        distance_from_start: f64,
    ) -> Point {
        let bearing = rhumb_bearing(as_geo(start.0), as_geo(end.0));
        let dest = rhumb_destination_m(as_geo(start.0), bearing, distance_from_start);
        Point::new(dest.lon, dest.lat)
    }
}
