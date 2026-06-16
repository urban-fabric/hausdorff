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

use super::{Crs, Datum, DatumShift, GeoDim, Helmert, ProjectionKind};
use crate::earth::projection::{
    AlbersEqualArea, LambertConformalConic, PolarStereographic, TransverseMercator, WebMercator,
};
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
                let base = node.child("BASEGEOGCRS").or_else(|| node.child("BASEGEODCRS"))?;
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

/// The datum of a geographic/geocentric (or base) node: the **registry** datum
/// when the node carries an EPSG code this crate knows (keeping its shift to
/// WGS84), else built from the `ELLIPSOID` with an identity shift.
fn datum_of(node: &WktNode) -> Option<Datum> {
    if let Some(code) = node.epsg() {
        if let Some(crs) = Crs::from_epsg(code) {
            return Some(crs.datum());
        }
    }
    let ell = node.find_deep("ELLIPSOID").or_else(|| node.find_deep("SPHEROID"))?;
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
    let method = conv.child("METHOD").or_else(|| conv.child("PROJECTION"))?;
    let param = |name: &str| -> Option<f64> {
        conv.children("PARAMETER")
            .find(|p| p.str_arg(0).is_some_and(|n| n.eq_ignore_ascii_case(name)))
            .and_then(|p| p.num_arg(1))
    };
    Some(match method.str_arg(0)?.to_ascii_lowercase().as_str() {
        "transverse mercator" => ProjectionKind::TransverseMercator(TransverseMercator {
            ellipsoid,
            central_meridian: param("Longitude of natural origin")?,
            lat_origin: param("Latitude of natural origin").unwrap_or(0.0),
            k0: param("Scale factor at natural origin").unwrap_or(1.0),
            false_easting: param("False easting").unwrap_or(0.0),
            false_northing: param("False northing").unwrap_or(0.0),
        }),
        "lambert conic conformal (2sp)" => {
            ProjectionKind::LambertConformalConic(LambertConformalConic {
                ellipsoid,
                lat1: param("Latitude of 1st standard parallel")?,
                lat2: param("Latitude of 2nd standard parallel")?,
                lat0: param("Latitude of false origin").unwrap_or(0.0),
                lon0: param("Longitude of false origin")?,
                false_easting: param("Easting at false origin").unwrap_or(0.0),
                false_northing: param("Northing at false origin").unwrap_or(0.0),
            })
        }
        "lambert conic conformal (1sp)" => {
            let lat0 = param("Latitude of natural origin")?;
            ProjectionKind::LambertConformalConic(LambertConformalConic {
                ellipsoid,
                lat1: lat0,
                lat2: lat0,
                lat0,
                lon0: param("Longitude of natural origin")?,
                false_easting: param("False easting").unwrap_or(0.0),
                false_northing: param("False northing").unwrap_or(0.0),
            })
        }
        "albers equal area" => ProjectionKind::AlbersEqualArea(AlbersEqualArea {
            ellipsoid,
            lat1: param("Latitude of 1st standard parallel")?,
            lat2: param("Latitude of 2nd standard parallel")?,
            lat0: param("Latitude of false origin").unwrap_or(0.0),
            lon0: param("Longitude of false origin")?,
            false_easting: param("Easting at false origin").unwrap_or(0.0),
            false_northing: param("Northing at false origin").unwrap_or(0.0),
        }),
        "polar stereographic (variant b)" => {
            let lat_ts = param("Latitude of standard parallel")?;
            ProjectionKind::PolarStereographic(PolarStereographic {
                ellipsoid,
                lat_ts,
                lon0: param("Longitude of origin").unwrap_or(0.0),
                north: lat_ts > 0.0,
                false_easting: param("False easting").unwrap_or(0.0),
                false_northing: param("False northing").unwrap_or(0.0),
            })
        }
        "popular visualisation pseudo mercator" | "mercator (variant a)"
        | "mercator (variant b)" => {
            ProjectionKind::WebMercator(WebMercator { radius: ellipsoid.a })
        }
        _ => return None,
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
            .find(|id| id.str_arg(0).is_some_and(|a| a.eq_ignore_ascii_case("EPSG")))
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
