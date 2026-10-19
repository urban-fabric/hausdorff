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
    if let Some(code) = node.epsg()
        && let Some(crs) = Crs::from_epsg(code)
    {
        return Some(crs.datum());
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

mod parser;
use parser::WktNode;
