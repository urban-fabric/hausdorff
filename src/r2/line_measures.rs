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

use super::{Coord, LineString, MultiLineString, Point};

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
