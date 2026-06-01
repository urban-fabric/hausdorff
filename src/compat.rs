//! Compatibility boundary for external geometry crates and interchange formats.
//!
//! Downstream UrbanFabric crates should depend on `hausdorff` geometry types.
//! This module contains the temporary adapters to `geo`/`geo-types`/`geo-traits`
//! needed for WKB/WKT/GeoJSON and topology/overlay operations that have not yet
//! been ported to native `hausdorff` kernels.

use geo::{BooleanOps, Contains, GeodesicArea, Intersects, Relate};

use crate::r2::{
    Coord, Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    Rect, Triangle,
};

/// Encode a lon/lat coordinate as a geohash string (native; see
/// [`crate::geohash`]).
#[must_use]
pub fn encode_geohash(lon: f64, lat: f64, precision: usize) -> Option<String> {
    Some(crate::geohash::encode(lon, lat, precision))
}

/// Topological predicates delegated through the compatibility boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopologyPredicate {
    /// The geometries share at least one point.
    Intersects,
    /// The left geometry contains the right.
    Contains,
    /// The left geometry is within the right.
    Within,
    /// Every point of the right lies in the left.
    Covers,
    /// Every point of the left lies in the right.
    CoveredBy,
    /// The geometries have no point in common.
    Disjoint,
    /// The geometries touch only at their boundaries.
    Touches,
    /// The geometries cross.
    Crosses,
    /// The geometries are topologically equal.
    Equals,
    /// The geometries overlap at the same dimension.
    Overlaps,
}

// WKB / WKT / GeoJSON now route through hausdorff's own native codecs in
// [`crate::r2_format`]; the external `wkb`/`wkt`/`geojson` crates are kept only
// as dev-dependency conformance references.

/// Decode one WKB value into a hausdorff geometry.
#[must_use]
pub fn decode_wkb(bytes: &[u8]) -> Option<Geometry> {
    crate::r2_format::wkb::from_wkb(bytes).ok()
}

/// Encode one hausdorff geometry as little-endian WKB.
#[must_use]
pub fn encode_wkb(geom: &Geometry) -> Option<Vec<u8>> {
    Some(crate::r2_format::wkb::to_wkb(geom))
}

/// Convert a hausdorff geometry to WKT text.
#[must_use]
pub fn geometry_to_wkt(geom: &Geometry) -> String {
    crate::r2_format::wkt::to_wkt(geom)
}

/// Parse WKT text into a hausdorff geometry.
#[must_use]
pub fn wkt_to_geometry(text: &str) -> Option<Geometry> {
    crate::r2_format::wkt::from_wkt(text).ok()
}

/// Convert a hausdorff geometry to a GeoJSON geometry string.
#[must_use]
pub fn geometry_to_geojson(geom: &Geometry) -> String {
    crate::r2_format::geojson::to_geojson(geom)
}

/// Parse a GeoJSON geometry string into a hausdorff geometry.
#[must_use]
pub fn geojson_to_geometry(text: &str) -> Option<Geometry> {
    crate::r2_format::geojson::from_geojson(text).ok()
}

/// Geodesic area in square metres, delegated through the compatibility boundary.
#[must_use]
pub fn geodesic_area_unsigned_m2(geom: &Geometry) -> f64 {
    to_geo_geometry(geom).geodesic_area_unsigned()
}

/// Evaluate a topology predicate through the compatibility boundary.
#[must_use]
pub fn topology_predicate(left: &Geometry, right: &Geometry, predicate: TopologyPredicate) -> bool {
    let left = to_geo_geometry(left);
    let right = to_geo_geometry(right);
    match predicate {
        TopologyPredicate::Intersects => left.intersects(&right),
        TopologyPredicate::Contains => left.contains(&right),
        _ => {
            let matrix = left.relate(&right);
            match predicate {
                TopologyPredicate::Within => matrix.is_within(),
                TopologyPredicate::Covers => matrix.is_covers(),
                TopologyPredicate::CoveredBy => matrix.is_coveredby(),
                TopologyPredicate::Disjoint => matrix.is_disjoint(),
                TopologyPredicate::Touches => matrix.is_touches(),
                TopologyPredicate::Crosses => matrix.is_crosses(),
                TopologyPredicate::Equals => matrix.is_equal_topo(),
                TopologyPredicate::Overlaps => matrix.is_overlaps(),
                TopologyPredicate::Intersects | TopologyPredicate::Contains => unreachable!(),
            }
        }
    }
}

/// Pairwise polygon intersection.
#[must_use]
pub fn intersection(left: &MultiPolygon, right: &MultiPolygon) -> MultiPolygon {
    from_geo_multi_polygon(&to_geo_multi_polygon(left).intersection(&to_geo_multi_polygon(right)))
}

/// Pairwise polygon union.
#[must_use]
pub fn union(left: &MultiPolygon, right: &MultiPolygon) -> MultiPolygon {
    from_geo_multi_polygon(&to_geo_multi_polygon(left).union(&to_geo_multi_polygon(right)))
}

/// Pairwise polygon difference.
#[must_use]
pub fn difference(left: &MultiPolygon, right: &MultiPolygon) -> MultiPolygon {
    from_geo_multi_polygon(&to_geo_multi_polygon(left).difference(&to_geo_multi_polygon(right)))
}

/// Pairwise polygon symmetric difference.
#[must_use]
pub fn xor(left: &MultiPolygon, right: &MultiPolygon) -> MultiPolygon {
    from_geo_multi_polygon(&to_geo_multi_polygon(left).xor(&to_geo_multi_polygon(right)))
}

fn to_geo_coord(coord: Coord) -> geo_types::Coord<f64> {
    geo_types::Coord {
        x: coord.x,
        y: coord.y,
    }
}

fn from_geo_coord(coord: geo_types::Coord<f64>) -> Coord {
    Coord::new(coord.x, coord.y)
}

fn to_geo_line_string(line_string: &LineString) -> geo_types::LineString<f64> {
    geo_types::LineString(line_string.0.iter().copied().map(to_geo_coord).collect())
}

fn from_geo_line_string(line_string: geo_types::LineString<f64>) -> LineString {
    LineString::new(line_string.0.into_iter().map(from_geo_coord).collect())
}

fn to_geo_polygon(polygon: &Polygon) -> geo_types::Polygon<f64> {
    geo_types::Polygon::new(
        to_geo_line_string(polygon.exterior()),
        polygon.interiors().iter().map(to_geo_line_string).collect(),
    )
}

fn from_geo_polygon(polygon: geo_types::Polygon<f64>) -> Polygon {
    let (exterior, interiors) = polygon.into_inner();
    Polygon::new(
        from_geo_line_string(exterior),
        interiors.into_iter().map(from_geo_line_string).collect(),
    )
}

fn to_geo_multi_polygon(multi_polygon: &MultiPolygon) -> geo_types::MultiPolygon<f64> {
    geo_types::MultiPolygon(multi_polygon.iter().map(to_geo_polygon).collect())
}

fn from_geo_multi_polygon(multi_polygon: &geo_types::MultiPolygon<f64>) -> MultiPolygon {
    MultiPolygon::new(
        multi_polygon
            .0
            .iter()
            .cloned()
            .map(from_geo_polygon)
            .collect(),
    )
}

/// Convert a hausdorff geometry to its `geo_types` equivalent (the bridge to the
/// external geo ecosystem; used by conformance tests and remaining geo-backed ops).
pub fn to_geo_geometry(geom: &Geometry) -> geo_types::Geometry<f64> {
    match geom {
        Geometry::Point(point) => {
            geo_types::Geometry::Point(geo_types::Point::new(point.x(), point.y()))
        }
        Geometry::Line(line) => geo_types::Geometry::Line(geo_types::Line::new(
            to_geo_coord(line.start),
            to_geo_coord(line.end),
        )),
        Geometry::LineString(line_string) => {
            geo_types::Geometry::LineString(to_geo_line_string(line_string))
        }
        Geometry::Polygon(polygon) => geo_types::Geometry::Polygon(to_geo_polygon(polygon)),
        Geometry::MultiPoint(multi_point) => {
            geo_types::Geometry::MultiPoint(geo_types::MultiPoint(
                multi_point
                    .iter()
                    .map(|point| geo_types::Point::new(point.x(), point.y()))
                    .collect(),
            ))
        }
        Geometry::MultiLineString(multi_line_string) => geo_types::Geometry::MultiLineString(
            geo_types::MultiLineString(multi_line_string.iter().map(to_geo_line_string).collect()),
        ),
        Geometry::MultiPolygon(multi_polygon) => {
            geo_types::Geometry::MultiPolygon(to_geo_multi_polygon(multi_polygon))
        }
        Geometry::Rect(rect) => geo_types::Geometry::Rect(geo_types::Rect::new(
            to_geo_coord(rect.min()),
            to_geo_coord(rect.max()),
        )),
        Geometry::Triangle(triangle) => geo_types::Geometry::Triangle(geo_types::Triangle::new(
            to_geo_coord(triangle.0),
            to_geo_coord(triangle.1),
            to_geo_coord(triangle.2),
        )),
        Geometry::GeometryCollection(collection) => geo_types::Geometry::GeometryCollection(
            geo_types::GeometryCollection(collection.iter().map(to_geo_geometry).collect()),
        ),
    }
}

/// Convert a `geo_types` geometry to its hausdorff equivalent.
pub fn from_geo_geometry(geom: geo_types::Geometry<f64>) -> Geometry {
    match geom {
        geo_types::Geometry::Point(point) => Geometry::Point(Point::new(point.x(), point.y())),
        geo_types::Geometry::Line(line) => Geometry::Line(Line::new(
            from_geo_coord(line.start),
            from_geo_coord(line.end),
        )),
        geo_types::Geometry::LineString(line_string) => {
            Geometry::LineString(from_geo_line_string(line_string))
        }
        geo_types::Geometry::Polygon(polygon) => Geometry::Polygon(from_geo_polygon(polygon)),
        geo_types::Geometry::MultiPoint(multi_point) => Geometry::MultiPoint(MultiPoint::new(
            multi_point
                .0
                .into_iter()
                .map(|point| Point::new(point.x(), point.y()))
                .collect(),
        )),
        geo_types::Geometry::MultiLineString(multi_line_string) => {
            Geometry::MultiLineString(MultiLineString::new(
                multi_line_string
                    .0
                    .into_iter()
                    .map(from_geo_line_string)
                    .collect(),
            ))
        }
        geo_types::Geometry::MultiPolygon(multi_polygon) => Geometry::MultiPolygon(
            MultiPolygon::new(multi_polygon.0.into_iter().map(from_geo_polygon).collect()),
        ),
        geo_types::Geometry::Rect(rect) => Geometry::Rect(Rect::new(
            from_geo_coord(rect.min()),
            from_geo_coord(rect.max()),
        )),
        geo_types::Geometry::Triangle(triangle) => Geometry::Triangle(Triangle::new(
            from_geo_coord(triangle.v1()),
            from_geo_coord(triangle.v2()),
            from_geo_coord(triangle.v3()),
        )),
        geo_types::Geometry::GeometryCollection(collection) => {
            Geometry::GeometryCollection(collection.0.into_iter().map(from_geo_geometry).collect())
        }
    }
}
