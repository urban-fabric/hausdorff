//! Parsing and emitting **PROJJSON** (the JSON encoding of a CRS) ↔ [`Crs`].
//!
//! PROJJSON is the same CRS model as WKT2 in JSON: a `type`
//! (`GeographicCRS`/`GeodeticCRS`/`ProjectedCRS`), a `datum`/`datum_ensemble` with
//! an `ellipsoid`, a `coordinate_system`, and — for a projected CRS — a `base_crs`
//! and a `conversion` (`method.name` + `parameters[]`). The method/parameter names
//! are the EPSG ones shared with WKT2, so the projection mapping is reused from
//! [`super::method`]. A small dependency-free JSON reader ([`Json`]) parses the
//! document (matching the hand-rolled WKT/PROJ-string parsers); `serde_json` is
//! intentionally not pulled into this crate.

use super::method::{self, ParamUnit};
use super::{Crs, Datum, DatumShift, GeoDim, Helmert, datum};
use hausdorff_base::earth::Ellipsoid;

impl Crs {
    /// Parse a **PROJJSON** string into a [`Crs`], or `None` if malformed or using
    /// a CRS/projection this crate does not model. A recognised `id`/`base_crs` id
    /// yields the registry datum (with its shift to WGS84); the EPSG code is kept.
    #[must_use]
    pub fn from_projjson(s: &str) -> Option<Crs> {
        let j = Json::parse(s)?;
        match j.get("type")?.as_str()? {
            "ProjectedCRS" | "DerivedProjectedCRS" => {
                let base = j.get("base_crs")?;
                let datum = datum_of(base)?;
                let conv = j.get("conversion")?;
                let method = conv.get("method")?.get("name")?.as_str()?;
                let params = conv.get("parameters")?.as_array()?;
                let projection = method::projection_from_method(method, datum.ellipsoid, |name| {
                    params
                        .iter()
                        .find(|p| p.get("name").and_then(Json::as_str) == Some(name))
                        .and_then(|p| p.get("value"))
                        .and_then(Json::as_f64)
                })?;
                Some(Crs::Projected {
                    datum,
                    projection,
                    epsg: epsg_of(&j),
                })
            }
            "GeographicCRS" | "GeodeticCRS" => {
                let datum = datum_of(&j)?;
                let cs = j.get("coordinate_system")?;
                let n = cs.get("axis")?.as_array()?.len();
                match cs.get("subtype")?.as_str()? {
                    "ellipsoidal" => Some(Crs::Geographic {
                        datum,
                        epsg: epsg_of(&j),
                        dim: if n >= 3 { GeoDim::D3 } else { GeoDim::D2 },
                    }),
                    "Cartesian" => Some(Crs::Geocentric {
                        datum,
                        epsg: epsg_of(&j),
                    }),
                    _ => None,
                }
            }
            _ => None,
        }
    }

    /// Serialize this CRS to a **PROJJSON** string — the inverse of
    /// [`from_projjson`](Self::from_projjson), valid for any reader. As with
    /// [`to_wkt2`](Self::to_wkt2) the datum shift is carried through this crate's
    /// reader via the `base_crs` EPSG id, not embedded as a transformation.
    #[must_use]
    pub fn to_projjson(&self) -> String {
        match *self {
            Crs::Geographic { datum, epsg, dim } => {
                let axes = if dim == GeoDim::D3 {
                    format!("{AXIS_LAT},{AXIS_LON},{AXIS_H}")
                } else {
                    format!("{AXIS_LAT},{AXIS_LON}")
                };
                let (name, dname, ename) = geographic_names(epsg);
                format!(
                    r#"{{"type":"GeographicCRS","name":"{name}","datum":{},"coordinate_system":{{"subtype":"ellipsoidal","axis":[{axes}]}}{}}}"#,
                    datum_json(datum, dname, ename),
                    id_json(epsg)
                )
            }
            Crs::Geocentric { datum, epsg } => format!(
                r#"{{"type":"GeodeticCRS","name":"geocentric","datum":{},"coordinate_system":{{"subtype":"Cartesian","axis":[{AXIS_X},{AXIS_Y},{AXIS_Z}]}}{}}}"#,
                datum_json(datum, "datum", "ellipsoid"),
                id_json(epsg)
            ),
            Crs::Projected {
                datum,
                projection,
                epsg,
            } => {
                let base = format!(
                    r#"{{"type":"GeographicCRS","name":"base","datum":{},"coordinate_system":{{"subtype":"ellipsoidal","axis":[{AXIS_LAT},{AXIS_LON}]}}{}}}"#,
                    datum_json(datum, "datum", "ellipsoid"),
                    id_json(datum::base_geographic_code(datum))
                );
                format!(
                    r#"{{"type":"ProjectedCRS","name":"projected","base_crs":{base},"conversion":{},"coordinate_system":{{"subtype":"Cartesian","axis":[{AXIS_E},{AXIS_N}]}}{}}}"#,
                    conversion_json(projection),
                    id_json(epsg)
                )
            }
        }
    }

    /// PROJJSON for a standard EPSG `code` (`None` if the registry does not model
    /// it) — convenience over [`from_epsg`](Self::from_epsg) +
    /// [`to_projjson`](Self::to_projjson). Geometry producers (e.g. a GeoParquet
    /// writer) name a CRS by code rather than hand-authoring the document.
    #[must_use]
    pub fn projjson_from_epsg(code: u32) -> Option<String> {
        Some(Crs::from_epsg(code)?.to_projjson())
    }
}

// PROJJSON axis fragments (valid `name`/`direction`/`unit` for any reader).
const AXIS_LAT: &str =
    r#"{"name":"Latitude","abbreviation":"lat","direction":"north","unit":"degree"}"#;
const AXIS_LON: &str =
    r#"{"name":"Longitude","abbreviation":"lon","direction":"east","unit":"degree"}"#;
const AXIS_H: &str =
    r#"{"name":"Ellipsoidal height","abbreviation":"h","direction":"up","unit":"metre"}"#;
const AXIS_X: &str =
    r#"{"name":"Geocentric X","abbreviation":"X","direction":"geocentricX","unit":"metre"}"#;
const AXIS_Y: &str =
    r#"{"name":"Geocentric Y","abbreviation":"Y","direction":"geocentricY","unit":"metre"}"#;
const AXIS_Z: &str =
    r#"{"name":"Geocentric Z","abbreviation":"Z","direction":"geocentricZ","unit":"metre"}"#;
const AXIS_E: &str = r#"{"name":"Easting","abbreviation":"E","direction":"east","unit":"metre"}"#;
const AXIS_N: &str = r#"{"name":"Northing","abbreviation":"N","direction":"north","unit":"metre"}"#;

/// `{"type":"GeodeticReferenceFrame","name":<name>,"ellipsoid":{"name":<ell>,…}}`.
/// `name`/`ell_name` are the human-readable datum/ellipsoid names; pass the generic
/// `"datum"`/`"ellipsoid"` placeholders when none is known.
fn datum_json(d: Datum, name: &str, ell_name: &str) -> String {
    let e = d.ellipsoid;
    let shape = if e.f == 0.0 {
        format!(r#""semi_major_axis":{},"semi_minor_axis":{}"#, e.a, e.a)
    } else {
        format!(
            r#""semi_major_axis":{},"inverse_flattening":{}"#,
            e.a,
            1.0 / e.f
        )
    };
    format!(
        r#"{{"type":"GeodeticReferenceFrame","name":"{name}","ellipsoid":{{"name":"{ell_name}",{shape}}}}}"#
    )
}

/// Human-readable `(crs, datum, ellipsoid)` names for well-known geographic EPSG
/// codes, so emitted PROJJSON carries real names (`"NAD83"`, …) rather than the
/// generic placeholders. Unknown codes fall back to the placeholders; readers key
/// on the EPSG `id` regardless, so this only improves legibility.
fn geographic_names(epsg: Option<u32>) -> (&'static str, &'static str, &'static str) {
    match epsg {
        Some(4326) => ("WGS 84", "World Geodetic System 1984", "WGS 84"),
        Some(4269) => ("NAD83", "North American Datum 1983", "GRS 1980"),
        Some(4258) => (
            "ETRS89",
            "European Terrestrial Reference System 1989",
            "GRS 1980",
        ),
        _ => ("geographic", "datum", "ellipsoid"),
    }
}

/// `,"id":{"authority":"EPSG","code":code}` or empty.
fn id_json(code: Option<u32>) -> String {
    code.map_or(String::new(), |c| {
        format!(r#","id":{{"authority":"EPSG","code":{c}}}"#)
    })
}

/// The `conversion` object: `method` (name + EPSG id) and `parameters[]`.
fn conversion_json(projection: super::ProjectionKind) -> String {
    let Some((method, mid, params)) = method::conversion_params(projection) else {
        return r#"{"name":"conversion","method":{"name":"Topocentric"}}"#.to_string();
    };
    let ps: Vec<String> = params
        .iter()
        .map(|&(name, value, ref unit)| {
            let u = match unit {
                ParamUnit::Angle => "degree",
                ParamUnit::Length => "metre",
                ParamUnit::Scale => "unity",
            };
            format!(r#"{{"name":"{name}","value":{value},"unit":"{u}"}}"#)
        })
        .collect();
    format!(
        r#"{{"name":"conversion","method":{{"name":"{method}","id":{{"authority":"EPSG","code":{mid}}}}},"parameters":[{}]}}"#,
        ps.join(",")
    )
}

/// The EPSG code from a node's `id` object (`{"authority":"EPSG","code":N}`).
fn epsg_of(node: &Json) -> Option<u32> {
    let id = node.get("id")?;
    if id.get("authority")?.as_str()?.eq_ignore_ascii_case("EPSG") {
        #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
        id.get("code")?.as_f64().map(|n| n as u32)
    } else {
        None
    }
}

/// The datum of a (base) CRS node: the **registry** datum when the node carries an
/// EPSG id this crate knows (keeping its shift to WGS84), else the bare
/// `ellipsoid` (identity shift). Mirrors the WKT2 reader.
fn datum_of(node: &Json) -> Option<Datum> {
    if let Some(code) = epsg_of(node) {
        if let Some(crs) = Crs::from_epsg(code) {
            return Some(crs.datum());
        }
    }
    let d = node.get("datum").or_else(|| node.get("datum_ensemble"))?;
    let ell = d.get("ellipsoid")?;
    let a = ell.get("semi_major_axis")?.as_f64()?;
    let f = if let Some(rf) = ell.get("inverse_flattening").and_then(Json::as_f64) {
        if rf == 0.0 { 0.0 } else { 1.0 / rf }
    } else if let Some(b) = ell.get("semi_minor_axis").and_then(Json::as_f64) {
        (a - b) / a
    } else {
        0.0
    };
    Some(Datum {
        ellipsoid: Ellipsoid { a, f },
        to_wgs84: DatumShift::Helmert(Helmert::IDENTITY),
    })
}

mod parser;
use parser::Json;
