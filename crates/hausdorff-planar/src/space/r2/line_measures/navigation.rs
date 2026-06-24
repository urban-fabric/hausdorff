use super::*;

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
