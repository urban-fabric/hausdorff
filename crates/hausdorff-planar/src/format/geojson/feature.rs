//! `Feature` and `FeatureCollection` — the GeoJSON documents that actually exist.
//!
//! The geometry codec beside this one handles the seven geometry objects, which is
//! what the *spec* is mostly about but not what files contain. Essentially every
//! `.geojson` in the wild is a `FeatureCollection`: a list of features, each
//! pairing one geometry with a bag of properties. Without this module hausdorff
//! could parse the GeoJSON nobody has and not the GeoJSON everybody has.
//!
//! ## Four details RFC 7946 is specific about, and parsers routinely get wrong
//!
//! - **`geometry` may be `null`.** A feature with attributes and no location is
//!   explicitly legal (§3.2), so the field is `Option<Geometry>` rather than being
//!   rejected or silently replaced with an empty geometry. A reader that treats
//!   null as an error fails on perfectly valid exports.
//! - **`id` may be a string *or* a number** (§3.2), so it is kept as a
//!   `serde_json::Value` rather than coerced. Coercing to `String` would make
//!   `1` and `"1"` the same feature; coercing to a number would drop UUIDs.
//! - **Foreign members are allowed anywhere** (§6.1) and must survive a round
//!   trip. Tools attach CRS hints, styling, and vendor metadata at the top level,
//!   and a codec that drops unknown keys silently destroys them.
//! - **`bbox` is kept verbatim and never recomputed.** It is advisory, it may
//!   disagree with the geometry, and it is 4 or 6 numbers depending on
//!   dimensionality. Recomputing it on write would be a different document than
//!   the one that was read; validating it would reject files other tools accept.
//!   So it is carried, not interpreted.
//!
//! ## Why `from_geojson_any` exists
//!
//! At a file boundary you usually do not know which of the three document kinds
//! you have — a bare geometry, a feature, or a collection — and the answer is in
//! the `type` member you have not parsed yet. Making the caller guess means
//! writing the try-each-in-turn ladder at every call site. [`GeoJson`] is that
//! ladder, once.

use serde_json::{Map, Value};

use super::{GeoJsonError, err, geom_to_value, value_to_geom};
use crate::space::r2::Geometry;

/// One GeoJSON feature: a geometry, its properties, and the optional members the
/// spec attaches to it.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Feature {
    /// The geometry, or `None` for an explicit `"geometry": null`.
    pub geometry: Option<Geometry>,
    /// The `properties` object. An absent or null `properties` reads as empty.
    pub properties: Map<String, Value>,
    /// The `id`, a string or a number, kept uncoerced.
    pub id: Option<Value>,
    /// The `bbox` as written: carried verbatim, never recomputed or validated.
    pub bbox: Option<Vec<f64>>,
    /// Members that are not part of the spec, preserved so a round trip is lossless.
    pub foreign: Map<String, Value>,
}

/// A GeoJSON `FeatureCollection`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct FeatureCollection {
    pub features: Vec<Feature>,
    /// The collection-level `bbox`, carried verbatim.
    pub bbox: Option<Vec<f64>>,
    /// Non-spec members at the collection level, preserved.
    pub foreign: Map<String, Value>,
}

/// Any of the three kinds of GeoJSON document.
///
/// The `type` member decides which, and it has to be read to find out — so this
/// is what [`from_geojson_any`] returns rather than making a caller try each
/// parser in turn.
#[derive(Clone, Debug, PartialEq)]
pub enum GeoJson {
    Geometry(Geometry),
    Feature(Feature),
    FeatureCollection(FeatureCollection),
}

/// Members this codec understands, and therefore must not copy into `foreign`.
const FEATURE_KNOWN: &[&str] = &["type", "geometry", "properties", "id", "bbox"];
const COLLECTION_KNOWN: &[&str] = &["type", "features", "bbox"];

fn collect_foreign(obj: &Map<String, Value>, known: &[&str]) -> Map<String, Value> {
    obj.iter()
        .filter(|(k, _)| !known.contains(&k.as_str()))
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect()
}

/// Parse a `bbox` member: an array of numbers, or `None` if absent.
///
/// A malformed `bbox` is an error rather than a silent `None`. It is advisory, but
/// a `bbox` of `["a", "b"]` means the document is not what it claims to be, and
/// quietly discarding it would hide that.
fn parse_bbox(v: Option<&Value>) -> Result<Option<Vec<f64>>, GeoJsonError> {
    match v {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Array(a)) => {
            let numbers: Option<Vec<f64>> = a.iter().map(Value::as_f64).collect();
            match numbers {
                Some(b) => Ok(Some(b)),
                None => err("bbox must contain only numbers"),
            }
        }
        Some(_) => err("bbox must be an array"),
    }
}

fn type_of(obj: &Map<String, Value>) -> Result<&str, GeoJsonError> {
    obj.get("type")
        .and_then(Value::as_str)
        .ok_or_else(|| GeoJsonError("missing or non-string \"type\" member".into()))
}

fn as_object(v: &Value) -> Result<&Map<String, Value>, GeoJsonError> {
    v.as_object()
        .ok_or_else(|| GeoJsonError("expected a JSON object".into()))
}

fn value_to_feature(v: &Value) -> Result<Feature, GeoJsonError> {
    let obj = as_object(v)?;
    if type_of(obj)? != "Feature" {
        return err("not a Feature");
    }
    // An explicit null geometry is legal; a *missing* geometry member is not, but
    // it is common enough in real exports that treating it as null is the useful
    // reading rather than a parse failure.
    let geometry = match obj.get("geometry") {
        None | Some(Value::Null) => None,
        Some(g) => Some(value_to_geom(g)?),
    };
    let properties = match obj.get("properties") {
        None | Some(Value::Null) => Map::new(),
        Some(Value::Object(p)) => p.clone(),
        Some(_) => return err("\"properties\" must be an object or null"),
    };
    let id = match obj.get("id") {
        None | Some(Value::Null) => None,
        Some(x @ (Value::String(_) | Value::Number(_))) => Some(x.clone()),
        Some(_) => return err("\"id\" must be a string or a number"),
    };
    Ok(Feature {
        geometry,
        properties,
        id,
        bbox: parse_bbox(obj.get("bbox"))?,
        foreign: collect_foreign(obj, FEATURE_KNOWN),
    })
}

fn value_to_collection(v: &Value) -> Result<FeatureCollection, GeoJsonError> {
    let obj = as_object(v)?;
    if type_of(obj)? != "FeatureCollection" {
        return err("not a FeatureCollection");
    }
    let features = match obj.get("features") {
        Some(Value::Array(a)) => a
            .iter()
            .map(value_to_feature)
            .collect::<Result<Vec<_>, _>>()?,
        // The spec requires the member; an empty collection is written as `[]`.
        _ => return err("\"features\" must be an array"),
    };
    Ok(FeatureCollection {
        features,
        bbox: parse_bbox(obj.get("bbox"))?,
        foreign: collect_foreign(obj, COLLECTION_KNOWN),
    })
}

fn feature_to_value(f: &Feature) -> Value {
    let mut obj = Map::new();
    obj.insert("type".into(), Value::from("Feature"));
    // Emitted even when `None`, as an explicit null: the member is required, and
    // omitting it produces a document some readers reject.
    obj.insert(
        "geometry".into(),
        f.geometry.as_ref().map_or(Value::Null, geom_to_value),
    );
    obj.insert("properties".into(), Value::Object(f.properties.clone()));
    if let Some(id) = &f.id {
        obj.insert("id".into(), id.clone());
    }
    if let Some(bbox) = &f.bbox {
        obj.insert("bbox".into(), Value::from(bbox.clone()));
    }
    for (k, v) in &f.foreign {
        obj.insert(k.clone(), v.clone());
    }
    Value::Object(obj)
}

fn collection_to_value(c: &FeatureCollection) -> Value {
    let mut obj = Map::new();
    obj.insert("type".into(), Value::from("FeatureCollection"));
    obj.insert(
        "features".into(),
        Value::Array(c.features.iter().map(feature_to_value).collect()),
    );
    if let Some(bbox) = &c.bbox {
        obj.insert("bbox".into(), Value::from(bbox.clone()));
    }
    for (k, v) in &c.foreign {
        obj.insert(k.clone(), v.clone());
    }
    Value::Object(obj)
}

impl Feature {
    /// A feature carrying just a geometry.
    #[must_use]
    pub fn from_geometry(geometry: Geometry) -> Self {
        Self {
            geometry: Some(geometry),
            ..Self::default()
        }
    }

    /// Parse one GeoJSON `Feature`.
    ///
    /// # Errors
    /// Returns [`GeoJsonError`] if the text is not a well-formed `Feature`.
    pub fn from_geojson(text: &str) -> Result<Self, GeoJsonError> {
        let v: Value = serde_json::from_str(text).map_err(|e| GeoJsonError(e.to_string()))?;
        value_to_feature(&v)
    }

    /// Serialize as a GeoJSON `Feature`.
    #[must_use]
    pub fn to_geojson(&self) -> String {
        feature_to_value(self).to_string()
    }
}

impl FeatureCollection {
    /// A collection of features with no collection-level members.
    #[must_use]
    pub fn new(features: Vec<Feature>) -> Self {
        Self {
            features,
            ..Self::default()
        }
    }

    /// Parse a GeoJSON `FeatureCollection`.
    ///
    /// # Errors
    /// Returns [`GeoJsonError`] if the text is not a well-formed collection.
    pub fn from_geojson(text: &str) -> Result<Self, GeoJsonError> {
        let v: Value = serde_json::from_str(text).map_err(|e| GeoJsonError(e.to_string()))?;
        value_to_collection(&v)
    }

    /// Serialize as a GeoJSON `FeatureCollection`.
    #[must_use]
    pub fn to_geojson(&self) -> String {
        collection_to_value(self).to_string()
    }

    /// Every non-null geometry, in order — the common reason to open the file.
    pub fn geometries(&self) -> impl Iterator<Item = &Geometry> + '_ {
        self.features.iter().filter_map(|f| f.geometry.as_ref())
    }
}

/// Parse any GeoJSON document, dispatching on its `type` member.
///
/// The entry point for a file boundary, where the kind of document is not known
/// until it has been read.
///
/// # Errors
/// Returns [`GeoJsonError`] if the text is not valid JSON, has no string `type`,
/// or that `type` is not a GeoJSON object kind.
pub fn from_geojson_any(text: &str) -> Result<GeoJson, GeoJsonError> {
    let v: Value = serde_json::from_str(text).map_err(|e| GeoJsonError(e.to_string()))?;
    let obj = as_object(&v)?;
    match type_of(obj)? {
        "Feature" => Ok(GeoJson::Feature(value_to_feature(&v)?)),
        "FeatureCollection" => Ok(GeoJson::FeatureCollection(value_to_collection(&v)?)),
        "Point" | "LineString" | "Polygon" | "MultiPoint" | "MultiLineString" | "MultiPolygon"
        | "GeometryCollection" => Ok(GeoJson::Geometry(value_to_geom(&v)?)),
        other => Err(GeoJsonError(format!("unknown GeoJSON type {other:?}"))),
    }
}

impl GeoJson {
    /// Serialize back to GeoJSON.
    #[must_use]
    pub fn to_geojson(&self) -> String {
        match self {
            GeoJson::Geometry(g) => geom_to_value(g).to_string(),
            GeoJson::Feature(f) => f.to_geojson(),
            GeoJson::FeatureCollection(c) => c.to_geojson(),
        }
    }

    /// Every geometry the document carries, whichever kind it is — so a caller
    /// that only wants geometry need not match on the kind at all.
    #[must_use]
    pub fn geometries(&self) -> Vec<&Geometry> {
        match self {
            GeoJson::Geometry(g) => vec![g],
            GeoJson::Feature(f) => f.geometry.iter().collect(),
            GeoJson::FeatureCollection(c) => c.geometries().collect(),
        }
    }
}
