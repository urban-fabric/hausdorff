//! `format::geojson` — native GeoJSON *geometry* codec over
//! [`crate::space::r2::Geometry`].
//!
//! The JSON layer uses `serde_json` (a generic JSON library, not a geo crate);
//! the GeoJSON geometry schema (Point/LineString/Polygon/Multi*/
//! GeometryCollection, RFC 7946) is handled here, so the `geojson` crate is not
//! needed. Coordinates beyond the first two (elevation) are ignored on read.

use serde_json::{Map, Value};

use crate::space::r2::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};

/// A GeoJSON parse error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeoJsonError(pub String);

impl std::fmt::Display for GeoJsonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "invalid GeoJSON: {}", self.0)
    }
}
impl std::error::Error for GeoJsonError {}

// ---- writer -----------------------------------------------------------------

/// Serialize a geometry to a GeoJSON geometry string.
#[must_use]
pub fn to_geojson(geom: &Geometry) -> String {
    geom_to_value(geom).to_string()
}

fn coord_v(c: Coord) -> Value {
    Value::Array(vec![Value::from(c.x), Value::from(c.y)])
}

fn coords_v(coords: &[Coord]) -> Value {
    Value::Array(coords.iter().map(|c| coord_v(*c)).collect())
}

fn poly_v(p: &Polygon) -> Value {
    let mut rings = vec![coords_v(&p.exterior().0)];
    rings.extend(p.interiors().iter().map(|h| coords_v(&h.0)));
    Value::Array(rings)
}

fn object(type_name: &str, key: &str, value: Value) -> Value {
    let mut m = Map::new();
    m.insert("type".into(), Value::String(type_name.into()));
    m.insert(key.into(), value);
    Value::Object(m)
}

fn geom_to_value(geom: &Geometry) -> Value {
    match geom {
        Geometry::Point(p) => object("Point", "coordinates", coord_v(p.0)),
        Geometry::Line(l) => object(
            "LineString",
            "coordinates",
            Value::Array(vec![coord_v(l.start), coord_v(l.end)]),
        ),
        Geometry::LineString(ls) => object("LineString", "coordinates", coords_v(&ls.0)),
        Geometry::Polygon(p) => object("Polygon", "coordinates", poly_v(p)),
        Geometry::MultiPoint(mp) => object(
            "MultiPoint",
            "coordinates",
            Value::Array(mp.iter().map(|p| coord_v(p.0)).collect()),
        ),
        Geometry::MultiLineString(mls) => object(
            "MultiLineString",
            "coordinates",
            Value::Array(mls.iter().map(|ls| coords_v(&ls.0)).collect()),
        ),
        Geometry::MultiPolygon(mp) => object(
            "MultiPolygon",
            "coordinates",
            Value::Array(mp.iter().map(poly_v).collect()),
        ),
        Geometry::Rect(r) => geom_to_value(&Geometry::Polygon(r.to_polygon())),
        Geometry::Triangle(t) => geom_to_value(&Geometry::Polygon(t.to_polygon())),
        Geometry::GeometryCollection(gs) => object(
            "GeometryCollection",
            "geometries",
            Value::Array(gs.iter().map(geom_to_value).collect()),
        ),
    }
}

// ---- reader -----------------------------------------------------------------

/// Parse a GeoJSON geometry string into a geometry.
///
/// # Errors
/// Returns [`GeoJsonError`] if the text is not a valid GeoJSON geometry.
pub fn from_geojson(text: &str) -> Result<Geometry, GeoJsonError> {
    let v: Value = serde_json::from_str(text).map_err(|e| GeoJsonError(e.to_string()))?;
    value_to_geom(&v)
}

fn err<T>(msg: &str) -> Result<T, GeoJsonError> {
    Err(GeoJsonError(msg.into()))
}

fn parse_coord(v: &Value) -> Result<Coord, GeoJsonError> {
    let a = v.as_array().ok_or_else(|| GeoJsonError("coordinate not an array".into()))?;
    match (a.first().and_then(Value::as_f64), a.get(1).and_then(Value::as_f64)) {
        (Some(x), Some(y)) => Ok(Coord::new(x, y)),
        _ => err("coordinate needs two numbers"),
    }
}

fn parse_coords(v: &Value) -> Result<Vec<Coord>, GeoJsonError> {
    v.as_array()
        .ok_or_else(|| GeoJsonError("expected coordinate list".into()))?
        .iter()
        .map(parse_coord)
        .collect()
}

fn parse_rings(v: &Value) -> Result<Vec<Vec<Coord>>, GeoJsonError> {
    v.as_array()
        .ok_or_else(|| GeoJsonError("expected ring list".into()))?
        .iter()
        .map(parse_coords)
        .collect()
}

fn parse_polygon(v: &Value) -> Result<Polygon, GeoJsonError> {
    let mut rings = parse_rings(v)?.into_iter();
    let exterior = rings.next().map_or_else(|| LineString::new(vec![]), LineString::new);
    Ok(Polygon::new(exterior, rings.map(LineString::new).collect()))
}

fn value_to_geom(v: &Value) -> Result<Geometry, GeoJsonError> {
    let type_name = v
        .get("type")
        .and_then(Value::as_str)
        .ok_or_else(|| GeoJsonError("missing \"type\"".into()))?;
    let coords = || {
        v.get("coordinates")
            .ok_or_else(|| GeoJsonError("missing \"coordinates\"".into()))
    };
    match type_name {
        "Point" => Ok(Geometry::Point(Point(parse_coord(coords()?)?))),
        "LineString" => Ok(Geometry::LineString(LineString::new(parse_coords(coords()?)?))),
        "Polygon" => Ok(Geometry::Polygon(parse_polygon(coords()?)?)),
        "MultiPoint" => Ok(Geometry::MultiPoint(MultiPoint::new(
            parse_coords(coords()?)?.into_iter().map(Point).collect(),
        ))),
        "MultiLineString" => Ok(Geometry::MultiLineString(MultiLineString::new(
            parse_rings(coords()?)?.into_iter().map(LineString::new).collect(),
        ))),
        "MultiPolygon" => {
            let polys = coords()?
                .as_array()
                .ok_or_else(|| GeoJsonError("MultiPolygon coordinates not an array".into()))?
                .iter()
                .map(parse_polygon)
                .collect::<Result<Vec<_>, _>>()?;
            Ok(Geometry::MultiPolygon(MultiPolygon::new(polys)))
        }
        "GeometryCollection" => {
            let gs = v
                .get("geometries")
                .and_then(Value::as_array)
                .ok_or_else(|| GeoJsonError("missing \"geometries\"".into()))?
                .iter()
                .map(value_to_geom)
                .collect::<Result<Vec<_>, _>>()?;
            Ok(Geometry::GeometryCollection(gs))
        }
        other => err(&format!("unknown geometry type {other:?}")),
    }
}
