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
                format!(
                    r#"{{"type":"GeographicCRS","name":"geographic","datum":{},"coordinate_system":{{"subtype":"ellipsoidal","axis":[{axes}]}}{}}}"#,
                    datum_json(datum),
                    id_json(epsg)
                )
            }
            Crs::Geocentric { datum, epsg } => format!(
                r#"{{"type":"GeodeticCRS","name":"geocentric","datum":{},"coordinate_system":{{"subtype":"Cartesian","axis":[{AXIS_X},{AXIS_Y},{AXIS_Z}]}}{}}}"#,
                datum_json(datum),
                id_json(epsg)
            ),
            Crs::Projected {
                datum,
                projection,
                epsg,
            } => {
                let base = format!(
                    r#"{{"type":"GeographicCRS","name":"base","datum":{},"coordinate_system":{{"subtype":"ellipsoidal","axis":[{AXIS_LAT},{AXIS_LON}]}}{}}}"#,
                    datum_json(datum),
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
}

// PROJJSON axis fragments (valid `name`/`direction`/`unit` for any reader).
const AXIS_LAT: &str = r#"{"name":"Latitude","abbreviation":"lat","direction":"north","unit":"degree"}"#;
const AXIS_LON: &str = r#"{"name":"Longitude","abbreviation":"lon","direction":"east","unit":"degree"}"#;
const AXIS_H: &str = r#"{"name":"Ellipsoidal height","abbreviation":"h","direction":"up","unit":"metre"}"#;
const AXIS_X: &str = r#"{"name":"Geocentric X","abbreviation":"X","direction":"geocentricX","unit":"metre"}"#;
const AXIS_Y: &str = r#"{"name":"Geocentric Y","abbreviation":"Y","direction":"geocentricY","unit":"metre"}"#;
const AXIS_Z: &str = r#"{"name":"Geocentric Z","abbreviation":"Z","direction":"geocentricZ","unit":"metre"}"#;
const AXIS_E: &str = r#"{"name":"Easting","abbreviation":"E","direction":"east","unit":"metre"}"#;
const AXIS_N: &str = r#"{"name":"Northing","abbreviation":"N","direction":"north","unit":"metre"}"#;

/// `{"type":"GeodeticReferenceFrame","name":"datum","ellipsoid":{…}}`.
fn datum_json(d: Datum) -> String {
    let e = d.ellipsoid;
    let shape = if e.f == 0.0 {
        format!(r#""semi_major_axis":{},"semi_minor_axis":{}"#, e.a, e.a)
    } else {
        format!(r#""semi_major_axis":{},"inverse_flattening":{}"#, e.a, 1.0 / e.f)
    };
    format!(
        r#"{{"type":"GeodeticReferenceFrame","name":"datum","ellipsoid":{{"name":"ellipsoid",{shape}}}}}"#
    )
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

/// A minimal parsed JSON value — enough to read PROJJSON, no external crate.
/// `Null`/`Bool` complete the JSON model (PROJJSON has e.g. `deprecated` booleans)
/// even though this reader never inspects their contents.
enum Json {
    Null,
    #[allow(dead_code)]
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    /// Parse a JSON document, or `None` if malformed.
    fn parse(s: &str) -> Option<Json> {
        let mut p = JsonParser {
            b: s.as_bytes(),
            i: 0,
        };
        let v = p.value()?;
        p.ws();
        // Trailing garbage is tolerated (PROJJSON documents are single values).
        Some(v)
    }

    /// The value for `key` in an object.
    fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(m) => m.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(n) => Some(*n),
            _ => None,
        }
    }

    fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Arr(a) => Some(a),
            _ => None,
        }
    }
}

/// A byte-cursor recursive-descent JSON parser.
struct JsonParser<'a> {
    b: &'a [u8],
    i: usize,
}

impl JsonParser<'_> {
    fn ws(&mut self) {
        while matches!(self.b.get(self.i), Some(b' ' | b'\n' | b'\t' | b'\r')) {
            self.i += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.b.get(self.i).copied()
    }

    fn value(&mut self) -> Option<Json> {
        self.ws();
        match self.peek()? {
            b'{' => self.object(),
            b'[' => self.array(),
            b'"' => Some(Json::Str(self.string()?)),
            b't' => self.lit("true", Json::Bool(true)),
            b'f' => self.lit("false", Json::Bool(false)),
            b'n' => self.lit("null", Json::Null),
            _ => self.number(),
        }
    }

    fn lit(&mut self, word: &str, val: Json) -> Option<Json> {
        if self.b[self.i..].starts_with(word.as_bytes()) {
            self.i += word.len();
            Some(val)
        } else {
            None
        }
    }

    fn object(&mut self) -> Option<Json> {
        self.i += 1; // '{'
        let mut m = Vec::new();
        self.ws();
        if self.peek() == Some(b'}') {
            self.i += 1;
            return Some(Json::Obj(m));
        }
        loop {
            self.ws();
            let key = self.string()?;
            self.ws();
            if self.peek()? != b':' {
                return None;
            }
            self.i += 1;
            m.push((key, self.value()?));
            self.ws();
            match self.peek()? {
                b',' => self.i += 1,
                b'}' => {
                    self.i += 1;
                    return Some(Json::Obj(m));
                }
                _ => return None,
            }
        }
    }

    fn array(&mut self) -> Option<Json> {
        self.i += 1; // '['
        let mut a = Vec::new();
        self.ws();
        if self.peek() == Some(b']') {
            self.i += 1;
            return Some(Json::Arr(a));
        }
        loop {
            a.push(self.value()?);
            self.ws();
            match self.peek()? {
                b',' => self.i += 1,
                b']' => {
                    self.i += 1;
                    return Some(Json::Arr(a));
                }
                _ => return None,
            }
        }
    }

    fn string(&mut self) -> Option<String> {
        if self.peek()? != b'"' {
            return None;
        }
        self.i += 1;
        let mut s = String::new();
        loop {
            let c = self.peek()?;
            self.i += 1;
            match c {
                b'"' => return Some(s),
                b'\\' => {
                    let e = self.peek()?;
                    self.i += 1;
                    match e {
                        b'"' => s.push('"'),
                        b'\\' => s.push('\\'),
                        b'/' => s.push('/'),
                        b'n' => s.push('\n'),
                        b't' => s.push('\t'),
                        b'r' => s.push('\r'),
                        b'b' => s.push('\u{8}'),
                        b'f' => s.push('\u{c}'),
                        b'u' => {
                            let hex = std::str::from_utf8(self.b.get(self.i..self.i + 4)?).ok()?;
                            let cp = u32::from_str_radix(hex, 16).ok()?;
                            self.i += 4;
                            s.push(char::from_u32(cp).unwrap_or('\u{fffd}'));
                        }
                        _ => return None,
                    }
                }
                // ASCII, or a raw UTF-8 byte we never inspect (PROJJSON is ASCII-escaped).
                _ => s.push(c as char),
            }
        }
    }

    fn number(&mut self) -> Option<Json> {
        let start = self.i;
        while matches!(self.peek(), Some(c)
            if c.is_ascii_digit() || matches!(c, b'-' | b'+' | b'.' | b'e' | b'E'))
        {
            self.i += 1;
        }
        std::str::from_utf8(&self.b[start..self.i])
            .ok()?
            .parse()
            .ok()
            .map(Json::Num)
    }
}
