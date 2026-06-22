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

mod parse;
use parse::{parse_record, read_header, read_i32_be, words_to_bytes};
