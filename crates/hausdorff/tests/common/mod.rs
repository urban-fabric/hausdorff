//! Test-only bridge to the `geo`/`geo_types` conformance oracle.
//!
//! hausdorff has **no runtime dependency** on external geometry crates — every
//! kernel/codec is native. This helper exists purely so the format and topology
//! conformance suites can cross-check the native implementations against `geo`.
//! It was previously `hausdorff::compat::{to_geo_geometry, from_geo_geometry}`;
//! moving it here keeps the `geo_types` conversion entirely in test scope and
//! `geo-types` a dev-dependency only.
#![allow(dead_code)]

use hausdorff::{
    Coord, Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    Rect, Triangle,
};

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

/// Convert a hausdorff geometry to its `geo_types` equivalent (the bridge to the
/// external geo ecosystem; used only by conformance tests).
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
