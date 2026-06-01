//! Compatibility boundary for external geometry crates and interchange formats.
//!
//! Downstream UrbanFabric crates depend on `hausdorff` geometry types and call
//! through this boundary. Almost everything here is now backed by native
//! `hausdorff` kernels — formats ([`crate::r2_format`]), geohash, overlay
//! ([`crate::r2::bool_ops`]), `Intersects`/`Disjoint`, and geodesic area.
//!
//! The **only** remaining use of the heavy `geo` crate is the DE-9IM topology
//! predicates beyond intersection (`Contains`/`Within`/`Covers`/`CoveredBy`/
//! `Touches`/`Crosses`/`Equals`/`Overlaps`) in [`topology_predicate`]. These are
//! exposed by galois's Polars API but are not consumed by any leaf application,
//! so a full native `Relate` matrix is intentionally out of scope (we don't need
//! the full `geo`). The lightweight `geo-types` bridge ([`to_geo_geometry`] /
//! [`from_geo_geometry`]) is kept deliberately for the WKB/shapefile I/O edges.

use geo::{Contains, Relate};

use crate::r2::{
    BooleanOps, Coord, Geometry, Line, LineString, MultiLineString, MultiPoint, MultiPolygon,
    Point, Polygon, Rect, Triangle,
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

/// Geodesic area in square metres (native spherical Chamberlain–Duquette; see
/// [`crate::r2::ChamberlainDuquetteArea`]). This is the spherical approximation,
/// not geo's ellipsoidal Karney area — within ~0.5%, which the workspace does
/// not depend on to tighter tolerance.
#[must_use]
pub fn geodesic_area_unsigned_m2(geom: &Geometry) -> f64 {
    use crate::r2::ChamberlainDuquetteArea;
    geom.chamberlain_duquette_unsigned_area()
}

/// Evaluate a topology predicate through the compatibility boundary.
///
/// `Intersects`/`Disjoint` are native ([`crate::r2::Intersects`]); the remaining
/// DE-9IM predicates still delegate to `geo`'s `Relate` until a native relate
/// matrix lands.
#[must_use]
pub fn topology_predicate(left: &Geometry, right: &Geometry, predicate: TopologyPredicate) -> bool {
    use crate::r2::Intersects as _;
    match predicate {
        TopologyPredicate::Intersects => return left.intersects(right),
        TopologyPredicate::Disjoint => return !left.intersects(right),
        _ => {}
    }
    let l = to_geo_geometry(left);
    let r = to_geo_geometry(right);
    match predicate {
        TopologyPredicate::Contains => l.contains(&r),
        TopologyPredicate::Within => l.relate(&r).is_within(),
        TopologyPredicate::Covers => l.relate(&r).is_covers(),
        TopologyPredicate::CoveredBy => l.relate(&r).is_coveredby(),
        TopologyPredicate::Touches => l.relate(&r).is_touches(),
        TopologyPredicate::Crosses => l.relate(&r).is_crosses(),
        TopologyPredicate::Equals => l.relate(&r).is_equal_topo(),
        TopologyPredicate::Overlaps => l.relate(&r).is_overlaps(),
        TopologyPredicate::Intersects | TopologyPredicate::Disjoint => unreachable!(),
    }
}

/// Pairwise polygon intersection (native Martinez–Rueda; see
/// [`crate::r2::bool_ops`]).
#[must_use]
pub fn intersection(left: &MultiPolygon, right: &MultiPolygon) -> MultiPolygon {
    left.intersection(right)
}

/// Pairwise polygon union (native).
#[must_use]
pub fn union(left: &MultiPolygon, right: &MultiPolygon) -> MultiPolygon {
    left.union(right)
}

/// Pairwise polygon difference (native).
#[must_use]
pub fn difference(left: &MultiPolygon, right: &MultiPolygon) -> MultiPolygon {
    left.difference(right)
}

/// Pairwise polygon symmetric difference (native).
#[must_use]
pub fn xor(left: &MultiPolygon, right: &MultiPolygon) -> MultiPolygon {
    left.xor(right)
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
