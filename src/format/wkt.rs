//! `wkt` — native Well-Known Text reader/writer over [`crate::space::r2::Geometry`].
//!
//! A self-contained codec for the OGC WKT geometry grammar (POINT, LINESTRING,
//! POLYGON, MULTIPOINT, MULTILINESTRING, MULTIPOLYGON, GEOMETRYCOLLECTION, and
//! their `EMPTY` forms) — no `geo`/`wkt`-crate round-trip. `r2::Rect`,
//! `r2::Triangle`, and `r2::Line` have no WKT keyword, so they serialize as the
//! equivalent POLYGON / LINESTRING.
//!
//! **Empty geometries are lossy in type.** The model has no typed-empty variant,
//! so every `… EMPTY` (e.g. `POINT EMPTY`, `LINESTRING EMPTY`) parses to an empty
//! [`MultiPoint`], and re-serializing yields `MULTIPOINT EMPTY`. The geometry is
//! preserved as empty, but its original type is not. (A typed-empty
//! representation would require extending [`Geometry`].)

use crate::space::r2::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};

/// A WKT parse error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WktError(pub String);

impl std::fmt::Display for WktError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "invalid WKT: {}", self.0)
    }
}
impl std::error::Error for WktError {}

/// Serialize a geometry to WKT.
#[must_use]
pub fn to_wkt(geom: &Geometry) -> String {
    let mut s = String::new();
    write_geom(&mut s, geom);
    s
}

fn num(s: &mut String, v: f64) {
    use std::fmt::Write;
    let _ = write!(s, "{v}");
}

fn coord(s: &mut String, c: Coord) {
    num(s, c.x);
    s.push(' ');
    num(s, c.y);
}

fn coord_list(s: &mut String, coords: &[Coord]) {
    s.push('(');
    for (i, c) in coords.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        coord(s, *c);
    }
    s.push(')');
}

fn rings(s: &mut String, poly: &Polygon) {
    s.push('(');
    coord_list(s, &poly.exterior().0);
    for hole in poly.interiors() {
        s.push_str(", ");
        coord_list(s, &hole.0);
    }
    s.push(')');
}

fn write_geom(s: &mut String, geom: &Geometry) {
    match geom {
        Geometry::Point(p) => {
            s.push_str("POINT (");
            coord(s, p.0);
            s.push(')');
        }
        Geometry::Line(l) => write_geom(
            s,
            &Geometry::LineString(LineString::new(vec![l.start, l.end])),
        ),
        Geometry::LineString(ls) => {
            s.push_str("LINESTRING ");
            if ls.0.is_empty() {
                s.push_str("EMPTY");
            } else {
                coord_list(s, &ls.0);
            }
        }
        Geometry::Polygon(p) => {
            s.push_str("POLYGON ");
            if p.exterior().0.is_empty() {
                s.push_str("EMPTY");
            } else {
                rings(s, p);
            }
        }
        Geometry::MultiPoint(mp) => {
            s.push_str("MULTIPOINT ");
            if mp.0.is_empty() {
                s.push_str("EMPTY");
            } else {
                s.push('(');
                for (i, p) in mp.iter().enumerate() {
                    if i > 0 {
                        s.push_str(", ");
                    }
                    s.push('(');
                    coord(s, p.0);
                    s.push(')');
                }
                s.push(')');
            }
        }
        Geometry::MultiLineString(mls) => {
            s.push_str("MULTILINESTRING ");
            if mls.0.is_empty() {
                s.push_str("EMPTY");
            } else {
                s.push('(');
                for (i, ls) in mls.iter().enumerate() {
                    if i > 0 {
                        s.push_str(", ");
                    }
                    coord_list(s, &ls.0);
                }
                s.push(')');
            }
        }
        Geometry::MultiPolygon(mp) => {
            s.push_str("MULTIPOLYGON ");
            if mp.0.is_empty() {
                s.push_str("EMPTY");
            } else {
                s.push('(');
                for (i, p) in mp.iter().enumerate() {
                    if i > 0 {
                        s.push_str(", ");
                    }
                    rings(s, p);
                }
                s.push(')');
            }
        }
        Geometry::Rect(r) => write_geom(s, &Geometry::Polygon(r.to_polygon())),
        Geometry::Triangle(t) => write_geom(s, &Geometry::Polygon(t.to_polygon())),
        Geometry::GeometryCollection(gs) => {
            s.push_str("GEOMETRYCOLLECTION ");
            if gs.is_empty() {
                s.push_str("EMPTY");
            } else {
                s.push('(');
                for (i, g) in gs.iter().enumerate() {
                    if i > 0 {
                        s.push_str(", ");
                    }
                    write_geom(s, g);
                }
                s.push(')');
            }
        }
    }
}

/// Parse WKT text into a geometry.
///
/// # Errors
/// Returns [`WktError`] if the text is not valid WKT.
pub fn from_wkt(text: &str) -> Result<Geometry, WktError> {
    let mut p = Parser {
        s: text.as_bytes(),
        i: 0,
    };
    let g = p.geometry()?;
    p.skip_ws();
    if p.i != p.s.len() {
        return Err(WktError("trailing input".into()));
    }
    Ok(g)
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn skip_ws(&mut self) {
        while self.i < self.s.len() && self.s[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }

    fn keyword(&mut self) -> String {
        self.skip_ws();
        let start = self.i;
        while self.i < self.s.len() && self.s[self.i].is_ascii_alphabetic() {
            self.i += 1;
        }
        String::from_utf8_lossy(&self.s[start..self.i]).to_ascii_uppercase()
    }

    fn eat(&mut self, c: u8) -> bool {
        self.skip_ws();
        if self.i < self.s.len() && self.s[self.i] == c {
            self.i += 1;
            true
        } else {
            false
        }
    }

    fn expect(&mut self, c: u8) -> Result<(), WktError> {
        if self.eat(c) {
            Ok(())
        } else {
            Err(WktError(format!("expected '{}'", c as char)))
        }
    }

    /// Peek whether the next keyword (without consuming) is EMPTY; consume it if so.
    fn empty(&mut self) -> bool {
        let save = self.i;
        if self.keyword() == "EMPTY" {
            true
        } else {
            self.i = save;
            false
        }
    }

    fn number(&mut self) -> Result<f64, WktError> {
        self.skip_ws();
        let start = self.i;
        while self.i < self.s.len() {
            let b = self.s[self.i];
            if b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'.' | b'e' | b'E') {
                self.i += 1;
            } else {
                break;
            }
        }
        std::str::from_utf8(&self.s[start..self.i])
            .ok()
            .and_then(|t| t.parse::<f64>().ok())
            .ok_or_else(|| WktError("expected number".into()))
    }

    fn coord(&mut self) -> Result<Coord, WktError> {
        let x = self.number()?;
        let y = self.number()?;
        Ok(Coord::new(x, y))
    }

    /// `( c, c, ... )`
    fn coord_seq(&mut self) -> Result<Vec<Coord>, WktError> {
        self.expect(b'(')?;
        let mut v = vec![self.coord()?];
        while self.eat(b',') {
            v.push(self.coord()?);
        }
        self.expect(b')')?;
        Ok(v)
    }

    fn geometry(&mut self) -> Result<Geometry, WktError> {
        let kw = self.keyword();
        match kw.as_str() {
            "POINT" => {
                if self.empty() {
                    return Ok(Geometry::MultiPoint(MultiPoint::new(vec![])));
                }
                self.expect(b'(')?;
                let c = self.coord()?;
                self.expect(b')')?;
                Ok(Geometry::Point(Point(c)))
            }
            "LINESTRING" => {
                if self.empty() {
                    return Ok(Geometry::LineString(LineString::new(vec![])));
                }
                Ok(Geometry::LineString(LineString::new(self.coord_seq()?)))
            }
            "POLYGON" => {
                if self.empty() {
                    return Ok(Geometry::Polygon(Polygon::new(
                        LineString::new(vec![]),
                        vec![],
                    )));
                }
                let rings = self.ring_seq()?;
                Ok(Geometry::Polygon(make_polygon(rings)))
            }
            "MULTIPOINT" => {
                if self.empty() {
                    return Ok(Geometry::MultiPoint(MultiPoint::new(vec![])));
                }
                self.expect(b'(')?;
                let mut pts = Vec::new();
                loop {
                    // both `(x y)` and bare `x y` forms are allowed
                    if self.eat(b'(') {
                        pts.push(Point(self.coord()?));
                        self.expect(b')')?;
                    } else {
                        pts.push(Point(self.coord()?));
                    }
                    if !self.eat(b',') {
                        break;
                    }
                }
                self.expect(b')')?;
                Ok(Geometry::MultiPoint(MultiPoint::new(pts)))
            }
            "MULTILINESTRING" => {
                if self.empty() {
                    return Ok(Geometry::MultiLineString(MultiLineString::new(vec![])));
                }
                self.expect(b'(')?;
                let mut lines = vec![LineString::new(self.coord_seq()?)];
                while self.eat(b',') {
                    lines.push(LineString::new(self.coord_seq()?));
                }
                self.expect(b')')?;
                Ok(Geometry::MultiLineString(MultiLineString::new(lines)))
            }
            "MULTIPOLYGON" => {
                if self.empty() {
                    return Ok(Geometry::MultiPolygon(MultiPolygon::new(vec![])));
                }
                self.expect(b'(')?;
                let mut polys = vec![make_polygon(self.ring_seq()?)];
                while self.eat(b',') {
                    polys.push(make_polygon(self.ring_seq()?));
                }
                self.expect(b')')?;
                Ok(Geometry::MultiPolygon(MultiPolygon::new(polys)))
            }
            "GEOMETRYCOLLECTION" => {
                if self.empty() {
                    return Ok(Geometry::GeometryCollection(vec![]));
                }
                self.expect(b'(')?;
                let mut gs = vec![self.geometry()?];
                while self.eat(b',') {
                    gs.push(self.geometry()?);
                }
                self.expect(b')')?;
                Ok(Geometry::GeometryCollection(gs))
            }
            other => Err(WktError(format!("unknown geometry type {other:?}"))),
        }
    }

    /// `( (ring), (ring), ... )`
    fn ring_seq(&mut self) -> Result<Vec<Vec<Coord>>, WktError> {
        self.expect(b'(')?;
        let mut rings = vec![self.coord_seq()?];
        while self.eat(b',') {
            rings.push(self.coord_seq()?);
        }
        self.expect(b')')?;
        Ok(rings)
    }
}

fn make_polygon(mut rings: Vec<Vec<Coord>>) -> Polygon {
    if rings.is_empty() {
        return Polygon::new(LineString::new(vec![]), vec![]);
    }
    let exterior = LineString::new(rings.remove(0));
    let interiors = rings.into_iter().map(LineString::new).collect();
    Polygon::new(exterior, interiors)
}
