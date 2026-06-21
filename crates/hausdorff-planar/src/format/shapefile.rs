//! `format::shapefile` - native ESRI Shapefile `.shp` geometry reader.
//!
//! This module reads the geometry stream only. It intentionally does not read
//! `.dbf` attribute sidecars or `.shx` indexes; callers that need attributes can
//! pair record order with their own DBF reader. The reader maps supported shape
//! records into Hausdorff's 2-D [`Geometry`] model. Z and M shape variants keep
//! their XY coordinates and ignore Z/M ordinates.

use std::fmt;
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use crate::space::r2::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect,
};

const FILE_CODE: i32 = 9994;
const HEADER_SIZE_BYTES: usize = 100;
const RECORD_HEADER_SIZE_BYTES: usize = 8;

/// ESRI shapefile shape type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ShapeType {
    NullShape,
    Point,
    PolyLine,
    Polygon,
    MultiPoint,
    PointZ,
    PolyLineZ,
    PolygonZ,
    MultiPointZ,
    PointM,
    PolyLineM,
    PolygonM,
    MultiPointM,
    MultiPatch,
}

impl ShapeType {
    /// The integer shape-type code used on disk.
    #[must_use]
    pub const fn code(self) -> i32 {
        match self {
            ShapeType::NullShape => 0,
            ShapeType::Point => 1,
            ShapeType::PolyLine => 3,
            ShapeType::Polygon => 5,
            ShapeType::MultiPoint => 8,
            ShapeType::PointZ => 11,
            ShapeType::PolyLineZ => 13,
            ShapeType::PolygonZ => 15,
            ShapeType::MultiPointZ => 18,
            ShapeType::PointM => 21,
            ShapeType::PolyLineM => 23,
            ShapeType::PolygonM => 25,
            ShapeType::MultiPointM => 28,
            ShapeType::MultiPatch => 31,
        }
    }

    fn from_code(code: i32) -> Result<Self, ShapefileError> {
        match code {
            0 => Ok(ShapeType::NullShape),
            1 => Ok(ShapeType::Point),
            3 => Ok(ShapeType::PolyLine),
            5 => Ok(ShapeType::Polygon),
            8 => Ok(ShapeType::MultiPoint),
            11 => Ok(ShapeType::PointZ),
            13 => Ok(ShapeType::PolyLineZ),
            15 => Ok(ShapeType::PolygonZ),
            18 => Ok(ShapeType::MultiPointZ),
            21 => Ok(ShapeType::PointM),
            23 => Ok(ShapeType::PolyLineM),
            25 => Ok(ShapeType::PolygonM),
            28 => Ok(ShapeType::MultiPointM),
            31 => Ok(ShapeType::MultiPatch),
            other => Err(ShapefileError::UnsupportedShapeType(other)),
        }
    }
}

/// File-level `.shp` header.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShapefileHeader {
    /// Total file length in bytes, as declared by the file header.
    pub file_length_bytes: usize,
    /// Shapefile format version. The standard value is 1000.
    pub version: i32,
    /// File-level shape type. Non-null records must match this value.
    pub shape_type: ShapeType,
    /// File-level XY bounding box.
    pub bbox: Rect,
    /// File-level Z range, if meaningful for the file's shape type.
    pub z_range: [f64; 2],
    /// File-level M range, if meaningful for the file's shape type.
    pub m_range: [f64; 2],
}

/// One geometry record from a `.shp` file.
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeRecord {
    /// Big-endian record number from the record header.
    pub number: i32,
    /// Record-level shape type.
    pub shape_type: ShapeType,
    /// The decoded geometry. `None` represents a NullShape record.
    pub geometry: Option<Geometry>,
}

/// Native `.shp` geometry reader.
pub struct ShapefileReader<R> {
    source: R,
    header: ShapefileHeader,
    bytes_read: usize,
}

impl<R: Read> ShapefileReader<R> {
    /// Create a reader and consume the 100-byte shapefile header.
    ///
    /// # Errors
    /// Returns [`ShapefileError`] for invalid headers, unsupported shape types,
    /// or I/O failures.
    pub fn new(mut source: R) -> Result<Self, ShapefileError> {
        let header = read_header(&mut source)?;
        Ok(Self {
            source,
            header,
            bytes_read: HEADER_SIZE_BYTES,
        })
    }

    /// The parsed file header.
    #[must_use]
    pub fn header(&self) -> &ShapefileHeader {
        &self.header
    }

    /// Read the next shape record, or `Ok(None)` at the declared file end.
    ///
    /// # Errors
    /// Returns [`ShapefileError`] for truncated records, invalid record sizes,
    /// unsupported shape types, or malformed geometry payloads.
    pub fn next_record(&mut self) -> Result<Option<ShapeRecord>, ShapefileError> {
        if self.bytes_read >= self.header.file_length_bytes {
            return Ok(None);
        }

        let number = read_i32_be(&mut self.source)?;
        let content_words = read_i32_be(&mut self.source)?;
        let content_len = words_to_bytes(content_words, "record content length")?;
        let mut payload = vec![0; content_len];
        self.source.read_exact(&mut payload)?;
        self.bytes_read = self
            .bytes_read
            .checked_add(RECORD_HEADER_SIZE_BYTES + content_len)
            .ok_or(ShapefileError::InvalidRecord("file position overflow"))?;

        let record = parse_record(number, &payload, self.header.shape_type)?;
        Ok(Some(record))
    }

    /// Read all remaining shape records.
    ///
    /// # Errors
    /// Propagates [`next_record`](Self::next_record) errors.
    pub fn read_all(mut self) -> Result<Vec<ShapeRecord>, ShapefileError> {
        let mut records = Vec::new();
        while let Some(record) = self.next_record()? {
            records.push(record);
        }
        Ok(records)
    }

    /// Read all remaining shape records as geometries.
    ///
    /// # Errors
    /// Propagates [`next_record`](Self::next_record) errors.
    pub fn read_geometries(self) -> Result<Vec<Option<Geometry>>, ShapefileError> {
        self.read_all()
            .map(|records| records.into_iter().map(|record| record.geometry).collect())
    }
}

/// Read all geometry records from a `.shp` byte stream.
///
/// # Errors
/// Returns [`ShapefileError`] for invalid headers, records, or I/O failures.
pub fn read_shapes<R: Read>(source: R) -> Result<Vec<ShapeRecord>, ShapefileError> {
    ShapefileReader::new(source)?.read_all()
}

/// Read all geometries from a `.shp` byte stream.
///
/// # Errors
/// Returns [`ShapefileError`] for invalid headers, records, or I/O failures.
pub fn read_geometries<R: Read>(source: R) -> Result<Vec<Option<Geometry>>, ShapefileError> {
    ShapefileReader::new(source)?.read_geometries()
}

/// Read all geometry records from a `.shp` file path.
///
/// # Errors
/// Returns [`ShapefileError`] for file-open, header, record, or geometry errors.
pub fn read_shapes_path(path: impl AsRef<Path>) -> Result<Vec<ShapeRecord>, ShapefileError> {
    let file = File::open(path)?;
    read_shapes(file)
}

/// Read all geometries from a `.shp` file path.
///
/// # Errors
/// Returns [`ShapefileError`] for file-open, header, record, or geometry errors.
pub fn read_geometries_path(
    path: impl AsRef<Path>,
) -> Result<Vec<Option<Geometry>>, ShapefileError> {
    let file = File::open(path)?;
    read_geometries(file)
}

/// Shapefile decode error.
#[derive(Debug)]
pub enum ShapefileError {
    Io(io::Error),
    InvalidFileCode(i32),
    UnsupportedShapeType(i32),
    UnsupportedGeometry(ShapeType),
    InvalidHeader(&'static str),
    InvalidRecord(&'static str),
    NegativeCount {
        field: &'static str,
        value: i32,
    },
    RecordTooShort {
        needed: usize,
        remaining: usize,
    },
    MismatchedRecordShapeType {
        header: ShapeType,
        record: ShapeType,
    },
    OrphanedInnerRing,
}

impl fmt::Display for ShapefileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ShapefileError::Io(e) => write!(f, "shapefile I/O error: {e}"),
            ShapefileError::InvalidFileCode(code) => {
                write!(f, "invalid shapefile file code {code}")
            }
            ShapefileError::UnsupportedShapeType(code) => {
                write!(f, "unsupported shapefile shape type {code}")
            }
            ShapefileError::UnsupportedGeometry(shape_type) => {
                write!(f, "unsupported shapefile geometry {shape_type:?}")
            }
            ShapefileError::InvalidHeader(msg) => write!(f, "invalid shapefile header: {msg}"),
            ShapefileError::InvalidRecord(msg) => write!(f, "invalid shapefile record: {msg}"),
            ShapefileError::NegativeCount { field, value } => {
                write!(f, "negative shapefile count {field}={value}")
            }
            ShapefileError::RecordTooShort { needed, remaining } => write!(
                f,
                "shapefile record too short: need {needed} bytes, have {remaining}"
            ),
            ShapefileError::MismatchedRecordShapeType { header, record } => write!(
                f,
                "record shape type {record:?} does not match header shape type {header:?}"
            ),
            ShapefileError::OrphanedInnerRing => {
                write!(f, "polygon inner ring appeared before any outer ring")
            }
        }
    }
}

impl std::error::Error for ShapefileError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ShapefileError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for ShapefileError {
    fn from(value: io::Error) -> Self {
        ShapefileError::Io(value)
    }
}

fn read_header<R: Read>(source: &mut R) -> Result<ShapefileHeader, ShapefileError> {
    let file_code = read_i32_be(source)?;
    if file_code != FILE_CODE {
        return Err(ShapefileError::InvalidFileCode(file_code));
    }

    skip_exact(source, 20)?;
    let file_length_bytes = words_to_bytes(read_i32_be(source)?, "file length")?;
    if file_length_bytes < HEADER_SIZE_BYTES {
        return Err(ShapefileError::InvalidHeader(
            "declared file length is smaller than header",
        ));
    }

    let version = read_i32_le(source)?;
    let shape_type = ShapeType::from_code(read_i32_le(source)?)?;
    let min_x = read_f64_le(source)?;
    let min_y = read_f64_le(source)?;
    let max_x = read_f64_le(source)?;
    let max_y = read_f64_le(source)?;
    let min_z = read_f64_le(source)?;
    let max_z = read_f64_le(source)?;
    let min_m = read_f64_le(source)?;
    let max_m = read_f64_le(source)?;

    Ok(ShapefileHeader {
        file_length_bytes,
        version,
        shape_type,
        bbox: Rect::new(Coord::new(min_x, min_y), Coord::new(max_x, max_y)),
        z_range: [min_z, max_z],
        m_range: [min_m, max_m],
    })
}

fn parse_record(
    number: i32,
    payload: &[u8],
    header_shape_type: ShapeType,
) -> Result<ShapeRecord, ShapefileError> {
    let mut r = SliceReader::new(payload);
    let shape_type = ShapeType::from_code(r.i32_le()?)?;
    if shape_type != ShapeType::NullShape && shape_type != header_shape_type {
        return Err(ShapefileError::MismatchedRecordShapeType {
            header: header_shape_type,
            record: shape_type,
        });
    }

    let geometry = match shape_type {
        ShapeType::NullShape => None,
        ShapeType::Point | ShapeType::PointM | ShapeType::PointZ => Some(read_point(&mut r)?),
        ShapeType::MultiPoint | ShapeType::MultiPointM | ShapeType::MultiPointZ => {
            Some(read_multipoint(&mut r)?)
        }
        ShapeType::PolyLine | ShapeType::PolyLineM | ShapeType::PolyLineZ => {
            Some(read_polyline(&mut r)?)
        }
        ShapeType::Polygon | ShapeType::PolygonM | ShapeType::PolygonZ => {
            Some(read_polygon(&mut r)?)
        }
        ShapeType::MultiPatch => return Err(ShapefileError::UnsupportedGeometry(shape_type)),
    };

    Ok(ShapeRecord {
        number,
        shape_type,
        geometry,
    })
}

fn read_point(r: &mut SliceReader<'_>) -> Result<Geometry, ShapefileError> {
    let x = r.f64_le()?;
    let y = r.f64_le()?;
    Ok(Geometry::Point(Point::new(x, y)))
}

fn read_multipoint(r: &mut SliceReader<'_>) -> Result<Geometry, ShapefileError> {
    skip_bbox(r)?;
    let num_points = nonnegative_usize(r.i32_le()?, "num_points")?;
    let points = read_points(r, num_points)?;
    Ok(Geometry::MultiPoint(MultiPoint::new(
        points.into_iter().map(Point).collect(),
    )))
}

fn read_polyline(r: &mut SliceReader<'_>) -> Result<Geometry, ShapefileError> {
    let parts = read_parts_shape(r)?;
    let lines = parts.into_iter().map(LineString::new).collect::<Vec<_>>();
    if lines.len() == 1 {
        Ok(Geometry::LineString(lines.into_iter().next().unwrap()))
    } else {
        Ok(Geometry::MultiLineString(MultiLineString::new(lines)))
    }
}

fn read_polygon(r: &mut SliceReader<'_>) -> Result<Geometry, ShapefileError> {
    let rings = read_parts_shape(r)?;
    let mut polygons: Vec<Polygon> = Vec::new();
    let mut current: Option<(LineString, Vec<LineString>)> = None;

    for ring in rings {
        let line = LineString::new(ring);
        if ring_is_inner(&line) {
            let Some((_, holes)) = current.as_mut() else {
                return Err(ShapefileError::OrphanedInnerRing);
            };
            holes.push(line);
        } else {
            if let Some((exterior, holes)) = current.take() {
                polygons.push(Polygon::new(exterior, holes));
            }
            current = Some((line, Vec::new()));
        }
    }

    if let Some((exterior, holes)) = current.take() {
        polygons.push(Polygon::new(exterior, holes));
    }

    if polygons.len() == 1 {
        Ok(Geometry::Polygon(polygons.into_iter().next().unwrap()))
    } else {
        Ok(Geometry::MultiPolygon(MultiPolygon::new(polygons)))
    }
}

fn read_parts_shape(r: &mut SliceReader<'_>) -> Result<Vec<Vec<Coord>>, ShapefileError> {
    skip_bbox(r)?;
    let num_parts = nonnegative_usize(r.i32_le()?, "num_parts")?;
    let num_points = nonnegative_usize(r.i32_le()?, "num_points")?;
    let mut starts = Vec::with_capacity(num_parts);
    for _ in 0..num_parts {
        let start = nonnegative_usize(r.i32_le()?, "part_start")?;
        if start > num_points {
            return Err(ShapefileError::InvalidRecord(
                "part start is greater than point count",
            ));
        }
        starts.push(start);
    }

    let points = read_points(r, num_points)?;
    let mut parts = Vec::with_capacity(num_parts);
    for i in 0..starts.len() {
        let start = starts[i];
        let end = starts.get(i + 1).copied().unwrap_or(num_points);
        if end < start || end > num_points {
            return Err(ShapefileError::InvalidRecord("invalid part range"));
        }
        parts.push(points[start..end].to_vec());
    }
    Ok(parts)
}

fn read_points(r: &mut SliceReader<'_>, n: usize) -> Result<Vec<Coord>, ShapefileError> {
    let mut points = Vec::with_capacity(n);
    for _ in 0..n {
        points.push(Coord::new(r.f64_le()?, r.f64_le()?));
    }
    Ok(points)
}

fn skip_bbox(r: &mut SliceReader<'_>) -> Result<(), ShapefileError> {
    r.skip(4 * size_of::<f64>())
}

fn ring_is_inner(ring: &LineString) -> bool {
    ring.0
        .windows(2)
        .map(|pts| (pts[1].x - pts[0].x) * (pts[1].y + pts[0].y))
        .sum::<f64>()
        / 2.0
        < 0.0
}

fn nonnegative_usize(value: i32, field: &'static str) -> Result<usize, ShapefileError> {
    if value < 0 {
        return Err(ShapefileError::NegativeCount { field, value });
    }
    Ok(value as usize)
}

fn words_to_bytes(words: i32, field: &'static str) -> Result<usize, ShapefileError> {
    let words = nonnegative_usize(words, field)?;
    words
        .checked_mul(2)
        .ok_or(ShapefileError::InvalidRecord("byte length overflow"))
}

fn skip_exact<R: Read>(source: &mut R, n: usize) -> Result<(), ShapefileError> {
    let mut buf = [0; 32];
    let mut remaining = n;
    while remaining > 0 {
        let take = remaining.min(buf.len());
        source.read_exact(&mut buf[..take])?;
        remaining -= take;
    }
    Ok(())
}

fn read_i32_be<R: Read>(source: &mut R) -> Result<i32, ShapefileError> {
    let mut b = [0; 4];
    source.read_exact(&mut b)?;
    Ok(i32::from_be_bytes(b))
}

fn read_i32_le<R: Read>(source: &mut R) -> Result<i32, ShapefileError> {
    let mut b = [0; 4];
    source.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

fn read_f64_le<R: Read>(source: &mut R) -> Result<f64, ShapefileError> {
    let mut b = [0; 8];
    source.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

struct SliceReader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> SliceReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], ShapefileError> {
        let end = self
            .pos
            .checked_add(n)
            .ok_or(ShapefileError::InvalidRecord("record cursor overflow"))?;
        if end > self.bytes.len() {
            return Err(ShapefileError::RecordTooShort {
                needed: n,
                remaining: self.bytes.len().saturating_sub(self.pos),
            });
        }
        let out = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn skip(&mut self, n: usize) -> Result<(), ShapefileError> {
        self.take(n).map(|_| ())
    }

    fn i32_le(&mut self) -> Result<i32, ShapefileError> {
        let b = self.take(4)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn f64_le(&mut self) -> Result<f64, ShapefileError> {
        let b = self.take(8)?;
        Ok(f64::from_le_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }
}
