//! Parsing **WKT2** (ISO 19162 well-known text) into a [`Crs`].
//!
//! [`Crs::from_wkt2`] parses the WKT2:2019 single-CRS forms — `GEOGCRS` /
//! `GEODCRS` (geographic or geocentric, distinguished by their `CS`) and `PROJCRS`
//! (a `CONVERSION` of a `BASEGEOGCRS`). WKT2 is a nested keyword tree
//! `KEYWORD["arg", number, NESTED[…], bareword]`; a small recursive-descent
//! parser builds that tree, then the CRS is read off it.
//!
//! The `CONVERSION`'s `METHOD` name selects the [`ProjectionKind`] and the
//! `PARAMETER[...]`s fill its fields (the EPSG parameter names — "Latitude of
//! natural origin", "Scale factor at natural origin", …). When the datum (or a
//! projected CRS's `BASEGEOGCRS`) carries an `ID["EPSG",NNNN]` this crate knows,
//! the **registry datum** is used — so its Helmert / NTv2-grid shift to WGS84 is
//! kept; otherwise the datum is built from the `ELLIPSOID` with an identity shift
//! (projection-only). Supported methods mirror [`Crs::from_proj_string`]:
//! Transverse Mercator, Lambert Conic Conformal (1SP/2SP), Albers Equal Area,
//! Polar Stereographic (variant B), and the pseudo-Mercator.

use super::method::{self, ParamUnit};
use super::{Crs, Datum, DatumShift, GeoDim, Helmert, ProjectionKind, datum};
use hausdorff_base::earth::Ellipsoid;

impl Crs {
    /// Parse a **WKT2** string into a [`Crs`], or `None` if it is malformed or uses
    /// a CRS/projection this crate does not model.
    ///
    /// The CRS keeps its `ID["EPSG",NNNN]` code (in [`epsg`](Self::epsg)) when
    /// present. A recognised datum/base code yields the registry datum (with its
    /// shift to WGS84); otherwise the datum is the bare `ELLIPSOID` (identity
    /// shift).
    #[must_use]
    pub fn from_wkt2(s: &str) -> Option<Crs> {
        let node = WktNode::parse(s)?;
        match node.keyword.to_ascii_uppercase().as_str() {
            "GEOGCRS" | "GEOGRAPHICCRS" | "GEODCRS" | "GEODETICCRS" => {
                let datum = datum_of(&node)?;
                let cs = node.child("CS")?;
                let cs_type = cs.str_arg(0)?;
                if cs_type.eq_ignore_ascii_case("ellipsoidal") {
                    let dim = if cs.num_arg(1) == Some(3.0) {
                        GeoDim::D3
                    } else {
                        GeoDim::D2
                    };
                    Some(Crs::Geographic {
                        datum,
                        epsg: node.epsg(),
                        dim,
                    })
                } else if cs_type.eq_ignore_ascii_case("Cartesian") {
                    Some(Crs::Geocentric {
                        datum,
                        epsg: node.epsg(),
                    })
                } else {
                    None
                }
            }
            "PROJCRS" | "PROJECTEDCRS" => {
                let base = node
                    .child("BASEGEOGCRS")
                    .or_else(|| node.child("BASEGEODCRS"))?;
                let datum = datum_of(base)?;
                let conv = node.child("CONVERSION")?;
                let projection = projection_of(conv, datum.ellipsoid)?;
                Some(Crs::Projected {
                    datum,
                    projection,
                    epsg: node.epsg(),
                })
            }
            _ => None,
        }
    }
}

impl Crs {
    /// Serialize this CRS to a **WKT2** string — the inverse of
    /// [`from_wkt2`](Self::from_wkt2). The output is valid WKT2 (units, `PRIMEM`,
    /// method/axis order), accepted by other readers (e.g. PROJ).
    ///
    /// When the datum is a registry datum its geographic `ID["EPSG",NNNN]` is
    /// emitted on the `BASEGEOGCRS`, so [`from_wkt2`](Self::from_wkt2) recovers the
    /// datum *with* its shift to WGS84. Note this datum shift is **not embedded**
    /// in the WKT itself — a self-contained shift would need a `BOUNDCRS` /
    /// `ABRIDGEDTRANSFORMATION` (not emitted), so an external reader that does not
    /// resolve the EPSG id treats the datum as shift-free. The *projection* is
    /// always exact. (For a fully self-contained datum shift, prefer
    /// [`to_proj_string`](Self::to_proj_string), which embeds `+towgs84`.)
    #[must_use]
    pub fn to_wkt2(&self) -> String {
        match *self {
            Crs::Geographic { datum, epsg, dim } => {
                let (n, axes) = if dim == GeoDim::D3 {
                    (
                        3,
                        format!(
                            r#"AXIS["latitude",north,ORDER[1],{DEG}],AXIS["longitude",east,ORDER[2],{DEG}],AXIS["ellipsoidal height",up,ORDER[3],{METRE}]"#
                        ),
                    )
                } else {
                    (
                        2,
                        format!(
                            r#"AXIS["latitude",north,ORDER[1],{DEG}],AXIS["longitude",east,ORDER[2],{DEG}]"#
                        ),
                    )
                };
                format!(
                    r#"GEOGCRS["geographic",{},{PRIMEM},CS[ellipsoidal,{n}],{axes}{}]"#,
                    datum_node(datum),
                    id_node(epsg)
                )
            }
            Crs::Geocentric { datum, epsg } => format!(
                r#"GEODCRS["geocentric",{},{PRIMEM},CS[Cartesian,3],AXIS["(X)",geocentricX,ORDER[1],{METRE}],AXIS["(Y)",geocentricY,ORDER[2],{METRE}],AXIS["(Z)",geocentricZ,ORDER[3],{METRE}]{}]"#,
                datum_node(datum),
                id_node(epsg)
            ),
            Crs::Projected {
                datum,
                projection,
                epsg,
            } => {
                let base = format!(
                    r#"BASEGEOGCRS["base",{},{PRIMEM}{}]"#,
                    datum_node(datum),
                    id_node(datum::base_geographic_code(datum))
                );
                format!(
                    r#"PROJCRS["projected",{base},{},CS[Cartesian,2],AXIS["(E)",east,ORDER[1],{METRE}],AXIS["(N)",north,ORDER[2],{METRE}]{}]"#,
                    conversion_node(projection),
                    id_node(epsg)
                )
            }
        }
    }
}

/// WKT2 unit / prime-meridian fragments shared by the serializer.
const DEG: &str = r#"ANGLEUNIT["degree",0.0174532925199433]"#;
const METRE: &str = r#"LENGTHUNIT["metre",1]"#;
const PRIMEM: &str = r#"PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]]"#;

/// `ELLIPSOID["…",a,1/f,LENGTHUNIT["metre",1]]` (inverse flattening; 0 for a sphere).
fn ellipsoid_node(e: Ellipsoid) -> String {
    let rf = if e.f == 0.0 { 0.0 } else { 1.0 / e.f };
    format!(
        r#"ELLIPSOID["ellipsoid",{},{},LENGTHUNIT["metre",1]]"#,
        e.a, rf
    )
}

/// `DATUM["…",ELLIPSOID[…]]`.
fn datum_node(d: Datum) -> String {
    format!(r#"DATUM["datum",{}]"#, ellipsoid_node(d.ellipsoid))
}

/// `,ID["EPSG",code]` or empty.
fn id_node(code: Option<u32>) -> String {
    code.map_or(String::new(), |c| format!(r#",ID["EPSG",{c}]"#))
}

/// The `CONVERSION["…",METHOD[…],PARAMETER[…]…]` for a projection. Parameters
/// carry their EPSG units (angle / length / scale) so the WKT2 is valid for any
/// reader, not just this crate's parser.
fn conversion_node(projection: ProjectionKind) -> String {
    // Topocentric has no standard WKT2 method — best-effort, not re-parseable.
    let Some((method, mid, params)) = method::conversion_params(projection) else {
        return r#"CONVERSION["conversion",METHOD["Topocentric"]]"#.to_string();
    };
    let body: Vec<String> = params
        .iter()
        .map(|&(name, value, ref unit)| {
            let u = match unit {
                ParamUnit::Angle => DEG,
                ParamUnit::Length => METRE,
                ParamUnit::Scale => r#"SCALEUNIT["unity",1]"#,
            };
            format!(r#"PARAMETER["{name}",{value},{u}]"#)
        })
        .collect();
    format!(
        r#"CONVERSION["conversion",METHOD["{method}",ID["EPSG",{mid}]],{}]"#,
        body.join(",")
    )
}

/// The datum of a geographic/geocentric (or base) node: the **registry** datum
/// when the node carries an EPSG code this crate knows (keeping its shift to
/// WGS84), else built from the `ELLIPSOID` with an identity shift.
fn datum_of(node: &WktNode) -> Option<Datum> {
    if let Some(code) = node.epsg() {
        if let Some(crs) = Crs::from_epsg(code) {
            return Some(crs.datum());
        }
    }
    let ell = node
        .find_deep("ELLIPSOID")
        .or_else(|| node.find_deep("SPHEROID"))?;
    let a = ell.num_arg(1)?;
    let rf = ell.num_arg(2)?;
    let f = if rf == 0.0 { 0.0 } else { 1.0 / rf };
    Some(Datum {
        ellipsoid: Ellipsoid { a, f },
        to_wgs84: DatumShift::Helmert(Helmert::IDENTITY),
    })
}

/// Build the [`ProjectionKind`] from a `CONVERSION` node's `METHOD` + `PARAMETER`s.
fn projection_of(conv: &WktNode, ellipsoid: Ellipsoid) -> Option<ProjectionKind> {
    let method = conv
        .child("METHOD")
        .or_else(|| conv.child("PROJECTION"))?
        .str_arg(0)?;
    method::projection_from_method(method, ellipsoid, |name| {
        conv.children("PARAMETER")
            .find(|p| p.str_arg(0).is_some_and(|n| n.eq_ignore_ascii_case(name)))
            .and_then(|p| p.num_arg(1))
    })
}

/// A parsed WKT node `KEYWORD[args]`. An arg is a quoted string, a number, a
/// bareword (kept as a [`WktArg::Str`]), or a nested node.
struct WktNode {
    keyword: String,
    args: Vec<WktArg>,
}

enum WktArg {
    Str(String),
    Num(f64),
    Node(WktNode),
}

impl WktNode {
    /// Parse a full WKT string into its root node.
    fn parse(s: &str) -> Option<WktNode> {
        let mut p = WktParser {
            b: s.as_bytes(),
            i: 0,
        };
        p.skip_ws();
        let node = p.node()?;
        Some(node)
    }

    /// The first direct child node with this keyword (case-insensitive).
    fn child(&self, keyword: &str) -> Option<&WktNode> {
        self.args.iter().find_map(|a| match a {
            WktArg::Node(n) if n.keyword.eq_ignore_ascii_case(keyword) => Some(n),
            _ => None,
        })
    }

    /// All direct child nodes with this keyword (case-insensitive).
    fn children<'a>(&'a self, keyword: &'a str) -> impl Iterator<Item = &'a WktNode> {
        self.args.iter().filter_map(move |a| match a {
            WktArg::Node(n) if n.keyword.eq_ignore_ascii_case(keyword) => Some(n),
            _ => None,
        })
    }

    /// The first node with this keyword anywhere in the subtree (this node first).
    fn find_deep(&self, keyword: &str) -> Option<&WktNode> {
        if self.keyword.eq_ignore_ascii_case(keyword) {
            return Some(self);
        }
        self.args.iter().find_map(|a| match a {
            WktArg::Node(n) => n.find_deep(keyword),
            _ => None,
        })
    }

    /// The `i`-th positional arg as a string (quoted string or bareword).
    fn str_arg(&self, i: usize) -> Option<&str> {
        match self.args.get(i)? {
            WktArg::Str(s) => Some(s),
            _ => None,
        }
    }

    /// The `i`-th positional arg as a number.
    fn num_arg(&self, i: usize) -> Option<f64> {
        match self.args.get(i)? {
            WktArg::Num(n) => Some(*n),
            _ => None,
        }
    }

    /// The EPSG code from this node's direct `ID["EPSG",code]` (or `AUTHORITY`).
    fn epsg(&self) -> Option<u32> {
        self.children("ID")
            .chain(self.children("AUTHORITY"))
            .find(|id| {
                id.str_arg(0)
                    .is_some_and(|a| a.eq_ignore_ascii_case("EPSG"))
            })
            .and_then(|id| id.num_arg(1))
            .map(|n| n as u32)
    }
}

/// A byte-cursor recursive-descent parser over a WKT string.
struct WktParser<'a> {
    b: &'a [u8],
    i: usize,
}

impl WktParser<'_> {
    fn skip_ws(&mut self) {
        while matches!(self.b.get(self.i), Some(b' ' | b'\n' | b'\t' | b'\r')) {
            self.i += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.b.get(self.i).copied()
    }

    /// A keyword / bareword token: `[A-Za-z0-9_]+`.
    fn word(&mut self) -> Option<String> {
        let start = self.i;
        while matches!(self.peek(), Some(c) if c.is_ascii_alphanumeric() || c == b'_') {
            self.i += 1;
        }
        if self.i == start {
            return None;
        }
        Some(String::from_utf8_lossy(&self.b[start..self.i]).into_owned())
    }

    /// A node `KEYWORD[args]` — assumes the keyword starts at the cursor.
    fn node(&mut self) -> Option<WktNode> {
        let keyword = self.word()?;
        self.skip_ws();
        if self.peek() != Some(b'[') {
            return None;
        }
        let args = self.args()?;
        Some(WktNode { keyword, args })
    }

    /// The bracketed, comma-separated argument list — assumes the cursor is at `[`.
    fn args(&mut self) -> Option<Vec<WktArg>> {
        self.i += 1; // consume '['
        let mut args = Vec::new();
        self.skip_ws();
        if self.peek() == Some(b']') {
            self.i += 1;
            return Some(args);
        }
        loop {
            self.skip_ws();
            args.push(self.arg()?);
            self.skip_ws();
            match self.peek()? {
                b',' => self.i += 1,
                b']' => {
                    self.i += 1;
                    return Some(args);
                }
                _ => return None,
            }
        }
    }

    /// One argument: quoted string, number, nested node, or bareword.
    fn arg(&mut self) -> Option<WktArg> {
        match self.peek()? {
            b'"' => Some(WktArg::Str(self.quoted())),
            c if c == b'-' || c == b'+' || c == b'.' || c.is_ascii_digit() => {
                let start = self.i;
                while matches!(self.peek(), Some(c)
                    if c.is_ascii_digit() || matches!(c, b'.' | b'-' | b'+' | b'e' | b'E'))
                {
                    self.i += 1;
                }
                std::str::from_utf8(&self.b[start..self.i])
                    .ok()?
                    .parse()
                    .ok()
                    .map(WktArg::Num)
            }
            c if c.is_ascii_alphabetic() => {
                let word = self.word()?;
                self.skip_ws();
                if self.peek() == Some(b'[') {
                    Some(WktArg::Node(WktNode {
                        keyword: word,
                        args: self.args()?,
                    }))
                } else {
                    Some(WktArg::Str(word)) // bareword (e.g. `ellipsoidal`, `north`)
                }
            }
            _ => None,
        }
    }

    /// A double-quoted string, with WKT's `""` → `"` escape. Assumes the cursor
    /// is at the opening quote.
    fn quoted(&mut self) -> String {
        self.i += 1; // opening quote
        let mut s = String::new();
        while let Some(c) = self.peek() {
            self.i += 1;
            if c == b'"' {
                if self.peek() == Some(b'"') {
                    self.i += 1;
                    s.push('"');
                } else {
                    break;
                }
            } else {
                s.push(c as char);
            }
        }
        s
    }
}
