//! `wkb` — native Well-Known Binary reader/writer over [`crate::space::r2::Geometry`].
//!
//! A self-contained codec for 2-D OGC/ISO WKB (geometry types 1–7), tolerant of
//! either byte order on read and of the ISO `Z`/`M` type offsets and the EWKB
//! SRID flag (the SRID and any Z/M ordinates are skipped). Writes little-endian
//! 2-D WKB. No `geo_types`/`wkb`-crate round trip.
//!
//! `r2::Rect`/`r2::Triangle`/`r2::Line` have no WKB type, so they encode as the
//! equivalent POLYGON / LINESTRING.

use crate::space::r2::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};

/// A WKB decode error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WkbError(pub String);

impl std::fmt::Display for WkbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "invalid WKB: {}", self.0)
    }
}
impl std::error::Error for WkbError {}

// ---- writer (little-endian) -------------------------------------------------

/// Encode a geometry as little-endian 2-D WKB.
#[must_use]
pub fn to_wkb(geom: &Geometry) -> Vec<u8> {
    // Size the buffer exactly once up front so the byte pushes never reallocate.
    let mut out = Vec::with_capacity(wkb_len(geom));
    write_geom(&mut out, geom);
    out
}

const WKB_HDR: usize = 5; // 1-byte order + u32 geometry type
const WKB_COUNT: usize = 4; // a u32 element/ring/coord count
const WKB_COORD: usize = 16; // two f64 ordinates

/// Bytes a ring of `n` coordinates occupies: its count plus the coordinates.
fn ring_len(n: usize) -> usize {
    WKB_COUNT + WKB_COORD * n
}

/// Exact serialized length of `geom` in little-endian 2-D WKB — mirrors
/// [`write_geom`] byte-for-byte so [`to_wkb`] can size its buffer once.
fn wkb_len(geom: &Geometry) -> usize {
    match geom {
        Geometry::Point(_) => WKB_HDR + WKB_COORD,
        Geometry::Line(_) => WKB_HDR + ring_len(2),
        Geometry::LineString(ls) => WKB_HDR + ring_len(ls.0.len()),
        Geometry::Polygon(p) => polygon_len(p),
        Geometry::MultiPoint(mp) => WKB_HDR + WKB_COUNT + mp.0.len() * (WKB_HDR + WKB_COORD),
        Geometry::MultiLineString(mls) => {
            WKB_HDR + WKB_COUNT + mls.iter().map(|ls| WKB_HDR + ring_len(ls.0.len())).sum::<usize>()
        }
        Geometry::MultiPolygon(mp) => WKB_HDR + WKB_COUNT + mp.iter().map(polygon_len).sum::<usize>(),
        Geometry::Rect(r) => polygon_len(&r.to_polygon()),
        Geometry::Triangle(t) => polygon_len(&t.to_polygon()),
        Geometry::GeometryCollection(gs) => {
            WKB_HDR + WKB_COUNT + gs.iter().map(wkb_len).sum::<usize>()
        }
    }
}

fn polygon_len(p: &Polygon) -> usize {
    if p.exterior().0.is_empty() {
        return WKB_HDR + WKB_COUNT;
    }
    WKB_HDR
        + WKB_COUNT
        + ring_len(p.exterior().0.len())
        + p.interiors().iter().map(|h| ring_len(h.0.len())).sum::<usize>()
}

fn hdr(out: &mut Vec<u8>, type_code: u32) {
    out.push(1); // NDR / little-endian
    out.extend_from_slice(&type_code.to_le_bytes());
}

fn put_u32(out: &mut Vec<u8>, v: usize) {
    // WKB counts are u32; a geometry with >4B elements can't be represented.
    debug_assert!(v <= u32::MAX as usize, "WKB count {v} exceeds u32");
    out.extend_from_slice(&(v as u32).to_le_bytes());
}

fn put_coord(out: &mut Vec<u8>, c: Coord) {
    out.extend_from_slice(&c.x.to_le_bytes());
    out.extend_from_slice(&c.y.to_le_bytes());
}

fn put_ring(out: &mut Vec<u8>, coords: &[Coord]) {
    put_u32(out, coords.len());
    for c in coords {
        put_coord(out, *c);
    }
}

fn write_polygon(out: &mut Vec<u8>, p: &Polygon) {
    hdr(out, 3);
    if p.exterior().0.is_empty() {
        put_u32(out, 0);
        return;
    }
    put_u32(out, 1 + p.interiors().len());
    put_ring(out, &p.exterior().0);
    for hole in p.interiors() {
        put_ring(out, &hole.0);
    }
}

fn write_geom(out: &mut Vec<u8>, geom: &Geometry) {
    match geom {
        Geometry::Point(p) => {
            hdr(out, 1);
            put_coord(out, p.0);
        }
        Geometry::Line(l) => {
            hdr(out, 2);
            put_u32(out, 2);
            put_coord(out, l.start);
            put_coord(out, l.end);
        }
        Geometry::LineString(ls) => {
            hdr(out, 2);
            put_ring(out, &ls.0);
        }
        Geometry::Polygon(p) => write_polygon(out, p),
        Geometry::MultiPoint(mp) => {
            hdr(out, 4);
            put_u32(out, mp.0.len());
            for p in mp.iter() {
                hdr(out, 1);
                put_coord(out, p.0);
            }
        }
        Geometry::MultiLineString(mls) => {
            hdr(out, 5);
            put_u32(out, mls.0.len());
            for ls in mls.iter() {
                hdr(out, 2);
                put_ring(out, &ls.0);
            }
        }
        Geometry::MultiPolygon(mp) => {
            hdr(out, 6);
            put_u32(out, mp.0.len());
            for p in mp.iter() {
                write_polygon(out, p);
            }
        }
        Geometry::Rect(r) => write_polygon(out, &r.to_polygon()),
        Geometry::Triangle(t) => write_polygon(out, &t.to_polygon()),
        Geometry::GeometryCollection(gs) => {
            hdr(out, 7);
            put_u32(out, gs.len());
            for g in gs {
                write_geom(out, g);
            }
        }
    }
}

// ---- reader -----------------------------------------------------------------

/// Decode WKB bytes into a geometry. Requires the bytes to be **exactly** one
/// WKB geometry — trailing bytes after the geometry are rejected.
///
/// # Errors
/// Returns [`WkbError`] if the bytes are not valid 2-D WKB, or if there are
/// leftover bytes after the geometry (use [`from_wkb_prefix`] to allow them).
pub fn from_wkb(bytes: &[u8]) -> Result<Geometry, WkbError> {
    let (geom, consumed) = from_wkb_prefix(bytes)?;
    if consumed != bytes.len() {
        return Err(WkbError(format!(
            "trailing bytes: parsed {consumed} of {} bytes",
            bytes.len()
        )));
    }
    Ok(geom)
}

/// Decode one WKB geometry from the start of `bytes`, returning it along with
/// the number of bytes consumed. Unlike [`from_wkb`], extra trailing bytes are
/// left unread rather than rejected — useful when a geometry is one field of a
/// larger record.
///
/// # Errors
/// Returns [`WkbError`] if the leading bytes are not valid 2-D WKB.
pub fn from_wkb_prefix(bytes: &[u8]) -> Result<(Geometry, usize), WkbError> {
    let mut r = Reader { b: bytes, i: 0 };
    let geom = r.geometry()?;
    Ok((geom, r.i))
}

struct Reader<'a> {
    b: &'a [u8],
    i: usize,
}

impl Reader<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8], WkbError> {
        if self.i + n > self.b.len() {
            return Err(WkbError("truncated".into()));
        }
        let s = &self.b[self.i..self.i + n];
        self.i += n;
        Ok(s)
    }

    fn order(&mut self) -> Result<bool, WkbError> {
        match self.take(1)?[0] {
            0 => Ok(false),
            1 => Ok(true),
            other => Err(WkbError(format!("bad byte order {other}"))),
        }
    }

    fn u32(&mut self, le: bool) -> Result<u32, WkbError> {
        let s = self.take(4)?;
        let a = [s[0], s[1], s[2], s[3]];
        Ok(if le {
            u32::from_le_bytes(a)
        } else {
            u32::from_be_bytes(a)
        })
    }

    fn f64(&mut self, le: bool) -> Result<f64, WkbError> {
        let s = self.take(8)?;
        let mut a = [0u8; 8];
        a.copy_from_slice(s);
        Ok(if le {
            f64::from_le_bytes(a)
        } else {
            f64::from_be_bytes(a)
        })
    }

    /// Read a header, returning `(byte_order_le, base_type, extra_ordinates)`.
    fn header(&mut self) -> Result<(bool, u32, usize), WkbError> {
        let le = self.order()?;
        let raw = self.u32(le)?;
        if raw & 0x2000_0000 != 0 {
            self.u32(le)?; // skip EWKB SRID
        }
        // EWKB Z/M flags or ISO 1000/2000-offset both add ordinates per coord.
        let ewkb_z = raw & 0x8000_0000 != 0;
        let ewkb_m = raw & 0x4000_0000 != 0;
        let base = (raw & 0x0000_ffff) % 1000;
        let iso = (raw & 0x0000_ffff) / 1000; // 0=2D, 1=Z, 2=M, 3=ZM
        let extra = usize::from(ewkb_z) + usize::from(ewkb_m) + usize::from(iso == 1 || iso == 2) + 2 * usize::from(iso == 3);
        Ok((le, base, extra))
    }

    fn coord(&mut self, le: bool, extra: usize) -> Result<Coord, WkbError> {
        let x = self.f64(le)?;
        let y = self.f64(le)?;
        for _ in 0..extra {
            self.f64(le)?; // discard Z / M
        }
        Ok(Coord::new(x, y))
    }

    fn coords(&mut self, le: bool, extra: usize) -> Result<Vec<Coord>, WkbError> {
        let n = self.u32(le)? as usize;
        (0..n).map(|_| self.coord(le, extra)).collect()
    }

    fn geometry(&mut self) -> Result<Geometry, WkbError> {
        let (le, base, extra) = self.header()?;
        match base {
            1 => Ok(Geometry::Point(Point(self.coord(le, extra)?))),
            2 => Ok(Geometry::LineString(LineString::new(self.coords(le, extra)?))),
            3 => Ok(Geometry::Polygon(self.polygon(le, extra)?)),
            4 => {
                let n = self.u32(le)? as usize;
                let mut pts = Vec::with_capacity(n);
                for _ in 0..n {
                    match self.geometry()? {
                        Geometry::Point(p) => pts.push(p),
                        _ => return Err(WkbError("MultiPoint member not a Point".into())),
                    }
                }
                Ok(Geometry::MultiPoint(MultiPoint::new(pts)))
            }
            5 => {
                let n = self.u32(le)? as usize;
                let mut lines = Vec::with_capacity(n);
                for _ in 0..n {
                    match self.geometry()? {
                        Geometry::LineString(ls) => lines.push(ls),
                        _ => return Err(WkbError("MultiLineString member not a LineString".into())),
                    }
                }
                Ok(Geometry::MultiLineString(MultiLineString::new(lines)))
            }
            6 => {
                let n = self.u32(le)? as usize;
                let mut polys = Vec::with_capacity(n);
                for _ in 0..n {
                    match self.geometry()? {
                        Geometry::Polygon(p) => polys.push(p),
                        _ => return Err(WkbError("MultiPolygon member not a Polygon".into())),
                    }
                }
                Ok(Geometry::MultiPolygon(MultiPolygon::new(polys)))
            }
            7 => {
                let n = self.u32(le)? as usize;
                let mut gs = Vec::with_capacity(n);
                for _ in 0..n {
                    gs.push(self.geometry()?);
                }
                Ok(Geometry::GeometryCollection(gs))
            }
            other => Err(WkbError(format!("unsupported geometry type {other}"))),
        }
    }

    fn polygon(&mut self, le: bool, extra: usize) -> Result<Polygon, WkbError> {
        let num_rings = self.u32(le)? as usize;
        if num_rings == 0 {
            return Ok(Polygon::new(LineString::new(vec![]), vec![]));
        }
        let exterior = LineString::new(self.coords(le, extra)?);
        let mut interiors = Vec::with_capacity(num_rings - 1);
        for _ in 1..num_rings {
            interiors.push(LineString::new(self.coords(le, extra)?));
        }
        Ok(Polygon::new(exterior, interiors))
    }
}
